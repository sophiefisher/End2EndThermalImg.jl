# Same per-candidate computation as one iteration of train_outer_optimization_distributed's loop
# body (src/optimize.jl), minus the geoms-gradient (not needed here -- geoms is fixed, whatever
# the outer optimization will start from) and minus any Adam update (resolve_initial_α_β below
# drives its own separate NLopt search over α/β instead). object/noise are fixed by the caller
# (sampled once, reused for every candidate (α, β) the search evaluates), unlike the outer
# optimization's own per-iteration resampling. geoms itself never changes across the sweep either
# (only α/β do), so unlike the outer loop, this does NOT call update_geoms_and_fftPSFs! itself --
# the caller (resolve_initial_α_β) does that once, up front, since fftPSFs don't depend on α/β at
# all and recomputing them on every candidate would be pure waste.
#
# No special-casing for T pinning at T_LOWER_BOUND here, unlike train_outer_optimization_distributed:
# confirmed numerically (see the conversation that led to this file) that both dB_dT and d2B_dT2
# stay finite at T_LOWER_BOUND (3.0 K) and well below it, so the Hessian/minres!/gradient chain
# can't produce NaN/Inf here the way it could with the old eps() bound. A candidate (α, β) bad
# enough to pin T just gets a correspondingly bad (but finite) L and gradient -- exactly the signal
# that should steer the search away from it, with no special handling needed.
function reconstruction_error_and_grad(object, noise, α, β, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp)
    @unpack php, imghp, rechp = jhp

    noisy_image = make_noisy_image_from_3D_distributed(object, weights, noise, imghp.noise_type, php, imghp)
    result = reconstruct_Tmap_and_zmap_distributed(noisy_image, weights, α, β, jhp; objective_scale = :auto)
    object_opt = result.object_opt

    fftPSFs = let
        fftPSFs = Matrix{Matrix{ComplexF64}}(undef, length(freqs), length(PSF_zcoords))
        for p in workers()
            grid_inds, shard = remotecall_fetch(p) do
                (_shard_grid_inds, _fftPSFs_shard)
            end
            for idx in grid_inds
                iF, iZ = Tuple(idx)
                fftPSFs[iF, iZ] = shard[(iF, iZ)]
            end
        end
        fftPSFs
    end

    _, active_linear_indices = get_active_variables(object_opt.zmap, imghp)
    free_linear_indices = setdiff(1:imghp.objN^2, active_linear_indices)
    b_T = vec(2 .* (object.Tmap .- object_opt.Tmap) ./ sum(abs2, object.Tmap))
    b_z = vec(2 .* (object.zmap .- object_opt.zmap) ./ sum(abs2, object.zmap))[free_linear_indices]
    b = [b_T; b_z]
    N_adjoint = 2 * imghp.objN^2 - length(active_linear_indices)

    A = build_hessian_operator_distributed(object, object_opt, fftPSFs, freqs, weights, driver_plans_PSF,
                                             noise, imghp.noise_type, α, β, php, imghp)
    Λ, minres_hist = IterativeSolvers.minres!(zeros(N_adjoint), A, b; log = true)

    Λ_T = Λ[1:imghp.objN^2]
    Λ_z = Λ[imghp.objN^2 + 1:end]
    dL_dα = 2 * dot(Λ_T, vec(object_opt.Tmap .- rechp.T_background))
    dL_dβ = 2 * dot(Λ_z, vec(object_opt.zmap .- rechp.z_middle)[free_linear_indices])

    L = L_reconstruction_error(object, object_opt)
    (; L, dL_dα, dL_dβ, object_opt, result, minres_hist)
end

# Resolves rechp.α_init/β_init: returns them unchanged if neither is :auto (no sweep at all), and
# otherwise runs an NLopt search in log-space (same reasoning as train_outer_optimization_distributed:
# a fixed learning step/tolerance should be relative, not absolute, across values that can span
# many orders of magnitude) over whichever of the two is :auto, holding the other fixed at its
# given literal value if only one is :auto. The objective is the reconstruction error of a single
# object/noise draw -- sampled once here, up front, and held fixed for every candidate the search
# evaluates (not resampled per iteration the way the outer optimization's own object/noise are).
#
# α_start/β_start are only the search's own starting point when the corresponding parameter is
# :auto (there's no literal rechp value to start from in that case) -- unrelated to whatever the
# literal value would otherwise be.
function resolve_initial_α_β(rechp, geoms, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp::JobHyperParams, PSF_scale;
    α_start = 1e-7, β_start = 1e-7, xtol_rel = 1e-6, maxeval = 100, verbose = true)
    @unpack imghp = jhp
    α_auto = rechp.α_init === :auto
    β_auto = rechp.β_init === :auto

    if !α_auto && !β_auto
        return (; α = rechp.α_init, β = rechp.β_init, history = NamedTuple[], object = nothing, noise = nothing, object_opt = nothing)
    end

    update_geoms_and_fftPSFs!(geoms, freqs, jhp.php, imghp, PSF_scale)
    object = get_object(imghp)
    noise = generate_noise(imghp)
    history = NamedTuple[]

    # x holds log(α) (if α_auto), then log(β) (if β_auto), in that order -- length 1 or 2
    # depending on which of the two are actually free.
    function αβ_from_x(x)
        i = 0
        α = α_auto ? exp(x[i += 1]) : rechp.α_init
        β = β_auto ? exp(x[i += 1]) : rechp.β_init
        (α, β)
    end

    n = count((α_auto, β_auto))
    opt = Opt(:LD_LBFGS, n)
    opt.min_objective = (x, grad) -> begin
        α, β = αβ_from_x(x)
        r = reconstruction_error_and_grad(object, noise, α, β, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp)
        if length(grad) > 0
            i = 0
            α_auto && (grad[i += 1] = r.dL_dα * α)
            β_auto && (grad[i += 1] = r.dL_dβ * β)
        end
        push!(history, (; iteration = length(history) + 1, α, β, L = r.L,
                         reconstruction_return_value = r.result.return_value,
                         minres_converged = r.minres_hist.isconverged))
        verbose && @info "resolve_initial_α_β iteration $(length(history))" α β L=r.L
        r.L
    end
    opt.xtol_rel = xtol_rel
    opt.maxeval = maxeval

    x0 = Float64[]
    α_auto && push!(x0, log(α_start))
    β_auto && push!(x0, log(β_start))

    (L_opt, x_opt, return_value) = NLopt.optimize(opt, x0)
    α_resolved, β_resolved = αβ_from_x(x_opt)
    verbose && @info "resolve_initial_α_β done" α=α_resolved β=β_resolved L_opt return_value

    # One more explicit evaluation at the final resolved (α, β) -- NLopt's own internal line
    # search may have evaluated a few points after the one it ultimately reports as optimal, so
    # object_opt from inside the loop above isn't guaranteed to correspond exactly to x_opt.
    r_final = reconstruction_error_and_grad(object, noise, α_resolved, β_resolved, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp)

    (; α = α_resolved, β = β_resolved, history, return_value, object, noise, object_opt = r_final.object_opt)
end
