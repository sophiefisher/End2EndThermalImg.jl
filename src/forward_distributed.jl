# Resident, sharded fftPSFs + distributed image formation.
#
# fftPSFs is only ever consumed as a linear sum over (freq, z) pairs to form an image
# (see make_image_at_z/make_image_from_3D in forward.jl), so rather than materializing and
# gathering the full (huge) fftPSFs array every call, each worker keeps its own shard of
# fftPSFs resident and we distribute the cheap image-formation reduction instead: only a small
# `object` goes out and a small `image` comes back per call, not the huge intermediate array.
#
# Note: none of the functions below use @everywhere. This file is `include`d as part of the
# End2EndThermalImg module itself, and every process already loads that module (via
# `using End2EndThermalImg`, see setup_cluster_workers in distributed.jl) -- so these functions
# already exist identically on every process, exactly like get_incident_field/get_fftPSF do.
# @everywhere is only correct for broadcasting definitions from an already-running interactive
# session (a notebook/script) to an already-running cluster; using it inside package source
# would make every process re-broadcast these definitions to the whole cluster on every package
# load, which is redundant at best and broken at worst (e.g. during initial precompilation,
# before any cluster exists at all).
#
# Also note: surrogates/n2f_kernels/plans_n2f/plans_PSF are accessed below via `Main.`-qualified
# names, not as function arguments or module globals. They're set up as per-worker-resident
# globals by the calling notebook/script's existing @everywhere setup, which runs in Main --
# a different module than End2EndThermalImg, where this file's functions are compiled -- so a
# bare reference would incorrectly resolve to (nonexistent) End2EndThermalImg-module globals.
# This is not optional for plans_n2f/plans_PSF in particular: FFTW plans wrap process-local C
# state and are not meaningfully serializable, so they must stay wherever they were computed
# locally on each worker and can never be passed across a remotecall/pmap boundary.
#
# Lifecycle: assign_shard_indices!() and init_resident_incidents!(...) run once per job.
# update_geoms_and_fftPSFs!(geoms, ...) runs once per outer (geoms) optimization iteration.
# make_image_from_3D_distributed/make_noisy_image_from_3D_distributed/
# reconstruction_objective_distributed/reconstruct_Tmap_and_zmap_distributed can then be called
# cheaply many times (e.g. the inner reconstruction loop) against the resident shards.

function set_shard_index!(shard_idx, nshards)
    global _shard_idx = shard_idx
    global _nshards = nshards
    nothing
end

function assign_shard_indices!()
    ws = workers()
    for (i, p) in enumerate(ws)
        remotecall_wait(set_shard_index!, p, i, length(ws))
    end
end

function init_incidents_shard!(freqs, PSF_zcoords, php)
    grid_inds = shard_grid_indices((length(freqs), length(PSF_zcoords)), _nshards, _shard_idx)
    # compute into a plain Vector (each thread writes its own disjoint index i, so no
    # synchronization is needed) and only build the Dict afterward, sequentially -- a plain Dict
    # isn't thread-safe for concurrent setindex! from multiple threads (insertion can trigger
    # internal rehashing), unlike every other @threads loop in this codebase, which writes into a
    # pre-allocated Vector/Matrix at disjoint indices instead.
    results = Vector{Matrix{ComplexF64}}(undef, length(grid_inds))
    @threads for i in eachindex(grid_inds)
        idx = grid_inds[i]
        results[i] = get_incident_field(freqs[idx[1]], PSF_zcoords[idx[2]], php)
    end
    vals = Dict{Tuple{Int,Int}, Matrix{ComplexF64}}(Tuple(grid_inds[i]) => results[i] for i in eachindex(grid_inds))
    global _shard_grid_inds = grid_inds
    global _incidents_shard = vals
    nothing
end

function init_resident_incidents!(freqs, PSF_zcoords, php)
    @sync for p in workers()
        @async remotecall_wait(init_incidents_shard!, p, freqs, PSF_zcoords, php)
    end
end

function recompute_fftPSFs_shard!(geoms, freqs, php, imghp, PSF_scale)
    # see init_incidents_shard! for why this computes into a plain Vector (thread-safe: each
    # thread writes its own disjoint index i) and only builds the Dict afterward, sequentially,
    # rather than mutating a Dict concurrently from inside @threads.
    results = Vector{Matrix{ComplexF64}}(undef, length(_shard_grid_inds))
    @threads for i in eachindex(_shard_grid_inds)
        idx = _shard_grid_inds[i]
        iF, iZ = Tuple(idx)
        tid = threadid()
        results[i] = get_fftPSF(
            get_PSF_at_freq_and_z(
                freqs[iF], _incidents_shard[(iF, iZ)], Main.surrogates[iF], geoms,
                Main.n2f_kernels[iF], Main.plans_n2f[tid], php, imghp, PSF_scale),
            Main.plans_PSF[tid])
    end
    shard = Dict{Tuple{Int,Int}, Matrix{ComplexF64}}(
        Tuple(_shard_grid_inds[i]) => results[i] for i in eachindex(_shard_grid_inds))
    global _fftPSFs_shard = shard
    global _fftPSFs_generation = (isdefined(@__MODULE__, :_fftPSFs_generation) ? _fftPSFs_generation : 0) + 1
    nothing
end

# call once per outer geoms iteration; recomputes and keeps resident (never gathers) each
# worker's fftPSFs shard for the new geoms.
function update_geoms_and_fftPSFs!(geoms, freqs, php, imghp, PSF_scale)
    # dispatch strictly one worker at a time (remotecall_wait blocks until that call returns)
    # rather than queuing several remotecalls before fetching any of them: `geoms` is a large
    # payload (2048x2048), and a list of remotecall(...) futures fetched afterward only LOOKS
    # sequential -- remotecall itself is non-blocking, so it still queues the send to every
    # worker essentially concurrently, which was observed (after enough repeated calls, e.g. in
    # a finite-difference sweep) to corrupt a worker's ClusterSerializer message stream
    # (deserialize TypeErrors on later, unrelated messages to that worker). Same issue
    # broadcast_value (distributed.jl) and Δgeoms_distributed (backward_distributed.jl) already
    # work around the same way. (This also drops the redundant Main.geoms broadcast the previous
    # version sent first: recompute_fftPSFs_shard! takes geoms as a plain argument, and nothing
    # in the package reads Main.geoms.)
    for p in workers()
        remotecall_wait(recompute_fftPSFs_shard!, p, geoms, freqs, php, imghp, PSF_scale)
    end
    nothing
end

function make_image_from_3D_shard(object, weights, php, imghp)
    # ignore_derivatives matches make_image_from_3D (forward.jl): these depend only on imghp,
    # not on `object`, but Zygote can't differentiate through get_PSF_zcoords' LinRange
    # construction, so it must be told not to try.
    δ_Δz = ignore_derivatives(() -> get_discretized_δ_function(imghp.smoothness_order, imghp.PSF_Δz))
    PSF_zcoords = ignore_derivatives(() -> get_PSF_zcoords(imghp))
    C_interp_3D = get_C_interp_3D(object.zmap, PSF_zcoords, δ_Δz, imghp)
    B = get_black_body_spectrum(object.Tmap, php)
    plan_PSF = Main.plans_PSF[1]
    sum(_shard_grid_inds) do idx
        iF, iZ = Tuple(idx)
        spectrum = B[iF] .* C_interp_3D[:, :, iZ]
        weights[iF] .* convolve_with_fftPSF_at_freq_and_z(_fftPSFs_shard[(iF, iZ)], spectrum, plan_PSF)
    end
end

function make_image_from_3D_distributed(object, weights, php, imghp)
    parts = pmap(p -> remotecall_fetch(make_image_from_3D_shard, p, object, weights, php, imghp), workers())
    sum(parts)
end

make_noisy_image_from_3D_distributed(object, weights, noise, noise_type::AbstractNoiseType, php, imghp) =
    apply_noise(noise_type, make_image_from_3D_distributed(object, weights, php, imghp), noise)

function reconstruction_objective_distributed(object_flat, noisy_image, weights, α, β, T_background, z_middle, php, imghp)
    object = unflatten_object(object_flat)
    image = make_image_from_3D_distributed(object, weights, php, imghp)
    error_image = sum((image .- noisy_image) .^ 2)
    regularization_τ = α * sum((object.Tmap .- T_background) .^ 2)
    regularization_ζ = β * sum((object.zmap .- z_middle) .^ 2)
    error_image + regularization_τ + regularization_ζ
end

# same as reconstruct_Tmap_and_zmap (forward.jl), targeting the distributed objective/image
# functions above; fftPSFs/plan_PSF are dropped since they're implicit/resident on the workers.
# objective_scale multiplies the objective (and hence its gradient) as seen by NLopt only, to lift
# the gradient above LD_LBFGS's hard-coded absolute gradient-norm stop threshold (roughly
# 1e-9..1e-8, independent of xtol_rel/ftol_rel): the data-misfit gradient w.r.t. T is ~1e-14 at
# the uniform initial guess, so with objective_scale = 1 NLopt returns SUCCESS after one
# evaluation without moving. The minimizer is unchanged (alpha/beta terms scale with it), and
# objective_opt/objective_history are returned in the original, unscaled units.
# objective_scale = :auto sets the scale to 10/f0, where f0 is the (unscaled) objective at the
# initial guess, so NLopt sees an objective that starts at 10. (One extra forward evaluation, no
# gradient.) A positive number is used as-is. The scale actually used is returned.
# ftol_rel/ftol_abs (rechp fields, like xtol_rel/maxeval) default to 0.0 (NLopt's convention for
# "ignore this criterion") rather than being left unset: when unset, NLopt falls back to its own
# internal defaults for them, which this function never accounted for -- a plausible explanation
# (under investigation) for cases where NLopt reported :SUCCESS after only a handful of iterations
# with object_opt barely moved from object_init, despite objective_scale lifting the gradient well
# above the documented ~1e-8 absolute floor at that point. With these at 0.0, only xtol_rel, the
# gradient-norm floor, and maxeval can stop this function, identically across every call.
# object_init: optional warm start (an (; Tmap, zmap) object), overriding rechp's own
# initialize_object(imghp, rechp). Needed for implicit-differentiation-style finite-difference
# checks across a perturbed α/β (see jupyter_notebooks/2026-10-16-validate-grad-alpha-beta.ipynb):
# re-initializing from the same generic start at every perturbed α/β can converge to a different
# local optimum each time (this reconstruction is not convex), so the implicit function theorem's
# premise -- tracking a single, continuously-varying stationary point -- silently breaks. Warm-
# starting from a nearby known optimum (e.g. the unperturbed baseline's own object_opt) keeps every
# trial in the same basin.
function reconstruct_Tmap_and_zmap_distributed(noisy_image, weights, α, β, jhp::JobHyperParams;
    iteration_print = 50, verbose = false, objective_scale = 1.0, object_init = nothing)
    verbose && @info "Starting object reconstruction"
    @unpack php, imghp, rechp = jhp
    @unpack T_background, z_middle, algorithm, maxeval, xtol_rel, ftol_rel, ftol_abs = rechp

    objective_history = Float64[]
    object_init = object_init === nothing ? initialize_object(imghp, rechp) : object_init
    object_init_flat = flatten_object(object_init)

    if objective_scale === :auto
        f0 = reconstruction_objective_distributed(object_init_flat, noisy_image, weights, α, β, T_background, z_middle, php, imghp)
        if isfinite(f0) && f0 > 0
            objective_scale = 10 / f0
        else
            @warn "objective_scale = :auto: initial objective f0 = $f0 is not positive and finite; falling back to objective_scale = 1.0"
            objective_scale = 1.0
        end
        verbose && @info "objective_scale = :auto" f0 objective_scale
    end

    # T (~T_background) and z (~z_middle) live on very different absolute scales (K vs the depth
    # unit, typically ~10^2 vs ~10^5) -- LD_LBFGS operates on the raw flattened [T; z] vector with
    # no per-variable scaling, so a single xtol_rel can be satisfied for the whole vector while one
    # block (z) is still far from a stationary point (its absolute step just looks small relative
    # to its own huge scale). Optimize over a rescaled vector where both blocks are O(1) instead:
    # to_physical is a plain elementwise multiply, so Zygote differentiates straight through it and
    # nlopt_wrap_objective_autodiff's gradient comes out correctly chain-ruled for the scaled
    # variable with no extra work. object_opt (returned below) is still in physical units.
    T_scale = T_background
    z_scale = abs(z_middle)
    scale_vec = [fill(T_scale, imghp.objN^2); fill(z_scale, imghp.objN^2)]
    to_scaled(x_physical) = x_physical ./ scale_vec
    to_physical(x_scaled) = x_scaled .* scale_vec

    opt = Opt(algorithm, 2 * imghp.objN^2)
    lower_bounds!(opt, to_scaled([fill(T_LOWER_BOUND, imghp.objN^2); fill(imghp.PSF_zlb, imghp.objN^2)]))
    upper_bounds!(opt, to_scaled([fill(Inf, imghp.objN^2); fill(imghp.PSF_zub, imghp.objN^2)]))
    objective_lambda = x_scaled -> objective_scale * reconstruction_objective_distributed(to_physical(x_scaled), noisy_image, weights, α, β, T_background, z_middle, php, imghp)
    objective_wrapped_lambda = (x, grad) -> nlopt_wrap_objective_autodiff(x, grad, objective_lambda, objective_history; iteration_print = iteration_print, verbose = verbose)
    min_objective!(opt, objective_wrapped_lambda)
    xtol_rel!(opt, xtol_rel)
    ftol_rel!(opt, ftol_rel)
    ftol_abs!(opt, ftol_abs)
    maxeval!(opt, maxeval)

    (objective_opt, x_opt_scaled, return_value) = NLopt.optimize!(opt, to_scaled(object_init_flat))
    objective_opt /= objective_scale
    objective_history ./= objective_scale
    verbose && @info "Done object reconstruction"
    verbose && @info "Optimization results" objective_opt return_value
    object_opt = unflatten_object(to_physical(x_opt_scaled))
    (; objective_opt, object_opt, return_value, objective_history, objective_scale)
end
