function setup_cluster_workers(; n_local_workers = max(Sys.CPU_THREADS ÷ 2, 1))
    if haskey(ENV, "SLURM_JOB_ID") || haskey(ENV, "SLURM_JOBID")
        # one worker per SLURM task (pair with --ntasks-per-node=1 for one worker per node),
        # each worker's own thread pool sized to its allocated cores for local @threads use
        cpus_per_task = parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "1"))
        addprocs(SlurmManager(); exeflags = "--threads=$(cpus_per_task)")
    else
        addprocs(n_local_workers; exeflags = "--threads=$(max(Sys.CPU_THREADS ÷ n_local_workers, 1))")
    end
    # can't use the @everywhere macro here: it expands to a :toplevel expression, which is
    # only valid at file/module top level, not nested inside a function body
    Distributed.remotecall_eval(Main, Distributed.procs(), :(using End2EndThermalImg))
    nworkers()
end

# contiguous block partition of 1:n into nshards pieces; returns the shard_idx-th piece (1-indexed)
function shard_indices(n::Integer, nshards::Integer, shard_idx::Integer)
    @assert 1 <= shard_idx <= nshards
    base, rem = divrem(n, nshards)
    start = (shard_idx - 1) * base + min(shard_idx - 1, rem) + 1
    len = base + (shard_idx <= rem ? 1 : 0)
    start:(start + len - 1)
end

# shard_indices, but over the flattened CartesianIndices of an n-dimensional grid
# (e.g. the (freq, z) work grid), so a shard isn't tied to a single axis
function shard_grid_indices(dims::Tuple, nshards::Integer, shard_idx::Integer)
    inds = CartesianIndices(dims)
    range = shard_indices(length(inds), nshards, shard_idx)
    inds[range]
end

# pushes an already-computed `value` to every worker as `Main.name`. The assignment expression
# carries `value` directly (the same mechanism `@everywhere x = $value` uses to broadcast a
# runtime value, as opposed to broadcasting code) -- sent one worker at a time.
# Distributed.remotecall_eval fans a broadcast out to every target concurrently via @sync/@async,
# which for large payloads (n2f_kernels is ~GB-scale) was observed to corrupt messages under
# concurrent send to multiple local workers; sending sequentially avoids that.
function broadcast_value(name::Symbol, value)
    ex = Expr(:(=), name, value)
    for p in Distributed.procs()
        Distributed.remotecall_wait(Core.eval, p, Main, ex)
    end
    value
end

# runs f(args...) once (on the calling process) and broadcasts the result to every worker, so
# expensive freq-only artifacts (n2f_kernels, surrogates) are computed once and shared rather
# than recomputed per node.
compute_and_broadcast(name::Symbol, f, args...) = broadcast_value(name, f(args...))

# Computes plans_n2f/plans_PSF on process p (via a direct remotecall dispatch) and stores them
# there as Main.plans_n2f/Main.plans_PSF. `n` is passed in explicitly (queried by the caller via a
# plain, reliable Threads.nthreads() call) rather than left for get_fft_plans to determine
# internally: see get_fft_plans's own docs (src/forward.jl) for why that distinction matters when
# this ends up dispatched to a worker from within another compiled package method.
function _setup_fft_plans_on_process!(php, imghp, n)
    plans_n2f, plans_PSF = get_fft_plans(php, imghp; n)
    Core.eval(Main, :(plans_n2f = $plans_n2f))
    Core.eval(Main, :(plans_PSF = $plans_PSF))
    nothing
end

# Sets Main.plans_n2f/Main.plans_PSF on every WORKER (never on the driver -- see below), and
# returns the driver's own plans_n2f, plans_PSF as a plain return value for the caller to bind
# under whatever name it likes:
#
#   plans_n2f, plans_PSF = setup_fft_plans_distributed!(php, imghp)   # driver's own, local names
#
# Two distinct bugs had to be found and fixed to get here, both confirmed via isolated A/B tests
# while debugging a silently-wrong-but-perfectly-reproducible factor-of-2 error in geoms gradients
# (with accompanying intermittent segfaults) that traced back to this setup step:
#
# 1. NOT `@everywhere plans_n2f, plans_PSF = get_fft_plans(php, imghp)`, and NOT dispatched via an
#    inline `do...end` closure: get_fft_plans's internal `for t in 1:nthreads()` loop, when
#    get_fft_plans is dispatched to a worker via @everywhere (or even via remotecall from within
#    another compiled package method), was found to have nthreads() silently see 1 instead of the
#    worker's real thread count -- even though the very same Threads.nthreads() call is completely
#    reliable when made directly (not nested inside get_fft_plans) in that same context. That
#    silently creates length-1 plan arrays on multi-threaded workers, so any later
#    @threads-dispatched task with threadid() > 1 (recompute_fftPSFs_shard!, Δgeoms_at_pair)
#    indexes past the end of that array -- undefined behavior. Fixed by querying nthreads() at the
#    call site (reliable) and passing it in explicitly (_setup_fft_plans_on_process! above), and
#    by dispatching via a plain, already-named top-level function (no closure).
#
# 2. The driver must NEVER bind its own plans to Main.plans_n2f/Main.plans_PSF (the same names
#    used on workers) -- doing so, by ANY means (Core.eval, a plain top-level assignment, before
#    or after the worker dispatch loop, whether or not it's wrapped in a package function),
#    reproducibly corrupts every worker's ALREADY-CORRECTLY-SET plans_n2f/plans_PSF back into a
#    length-1 array, even though nothing about the workers' own state or process is touched by the
#    driver's assignment. The mechanism is unclear -- something in Julia/Distributed.jl treats a
#    same-named global bound on multiple processes at once specially -- but the effect is exact
#    and 100% reproducible. Avoided entirely by never giving the driver's own plans that name;
#    they're returned instead, for the caller to bind under any other name.
function setup_fft_plans_distributed!(php, imghp)
    for p in Distributed.workers()
        n = Distributed.remotecall_fetch(Threads.nthreads, p)
        Distributed.remotecall_wait(_setup_fft_plans_on_process!, p, php, imghp, n)
    end
    get_fft_plans(php, imghp)
end
