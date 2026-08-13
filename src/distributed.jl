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

# runs f(args...) once (on the calling process) and pushes the result to every worker as
# `Main.name`, so expensive freq-only artifacts (n2f_kernels, surrogates) are computed once
# and shared rather than recomputed per node. The assignment expression carries `value`
# directly (the same mechanism `@everywhere x = $value` uses to broadcast a runtime value,
# as opposed to broadcasting code) -- sent one worker at a time. Distributed.remotecall_eval
# fans a broadcast out to every target concurrently via @sync/@async, which for large
# payloads (n2f_kernels is ~GB-scale) was observed to corrupt messages under concurrent
# send to multiple local workers; sending sequentially avoids that.
function compute_and_broadcast(name::Symbol, f, args...)
    value = f(args...)
    ex = Expr(:(=), name, value)
    for p in Distributed.procs()
        Distributed.remotecall_wait(Core.eval, p, Main, ex)
    end
    value
end
