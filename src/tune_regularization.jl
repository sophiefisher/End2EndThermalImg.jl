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
# T pinning at T_LOWER_BOUND is handled the same way as everywhere else now (get_free_and_active_
# indices, backward.jl), not specially -- a candidate (α, β) bad enough to pin T just gets a
# correspondingly bad (but finite and correct) L and gradient, exactly the signal that should
# steer the search away from it. The one genuinely degenerate case (every T and z pixel pinned,
# N_adjoint == 0) is guarded explicitly below, since that candidate is more plausible here than in
# train_outer_optimization_distributed: this search deliberately sweeps α/β across many orders of
# magnitude in log-space, including extreme values unlikely to arise from Adam's smaller steps.
function reconstruction_error_and_grad(object, noise, α, β, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp)
    @unpack php, imghp, rechp = jhp

    noisy_image = make_noisy_image_from_3D_distributed(object, weights, noise, imghp.noise_type, php, imghp)
    result = reconstruct_Tmap_and_zmap_distributed(noisy_image, weights, α, β, jhp; objective_scale = :auto)
    object_opt = result.object_opt

    idx = get_free_and_active_indices(object_opt, imghp)
    free_T, free_z = idx.free_T, idx.free_z
    n_active_T, n_active_z = length(idx.active_T), length(idx.active_z)
    N_adjoint = length(free_T) + length(free_z)
    L = L_reconstruction_error(object, object_opt)

    if N_adjoint == 0
        # Nothing is a free reconstruction variable left for α/β to influence -- no adjoint
        # equation to solve, so report zero sensitivity rather than attempting a degenerate 0x0
        # minres! solve (same reasoning as optimize.jl's own N_adjoint == 0 guard).
        @warn "Every T and z pixel pinned at its bound for (α, β) = ($α, $β) -- reporting zero dL_dα/dL_dβ" n_active_T n_active_z
        return (; L, dL_dα = 0.0, dL_dβ = 0.0, object_opt, result, minres_hist = nothing, n_active_T, n_active_z)
    end

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
        # If a worker died since assign_shard_indices! last partitioned the (freq, z) grid across
        # workers(), that worker's slice is permanently missing from every future gather (its
        # grid_inds are no longer reachable through any live worker) -- fail loudly here, with
        # which cells are missing, rather than a cryptic UndefRefError deep inside
        # build_hessian_operator_distributed/make_image_from_3D.
        missing_cells = [(iF, iZ) for iF in axes(fftPSFs, 1), iZ in axes(fftPSFs, 2) if !isassigned(fftPSFs, iF, iZ)]
        isempty(missing_cells) || error("fftPSFs gather is missing $(length(missing_cells)) (freq, z) cells " *
            "(e.g. $(first(missing_cells))) -- a worker likely died since assign_shard_indices! was last called; " *
            "restart cluster workers and re-run setup (assign_shard_indices!/init_resident_incidents!/" *
            "update_geoms_and_fftPSFs!) before retrying. workers() = $(workers())")
        fftPSFs
    end

    b_T = vec(2 .* (object.Tmap .- object_opt.Tmap) ./ sum(abs2, object.Tmap))[free_T]
    b_z = vec(2 .* (object.zmap .- object_opt.zmap) ./ sum(abs2, object.zmap))[free_z]
    b = [b_T; b_z]

    A = build_hessian_operator_distributed(object, object_opt, fftPSFs, freqs, weights, driver_plans_PSF,
                                             noise, imghp.noise_type, α, β, php, imghp)
    Λ, minres_hist = IterativeSolvers.minres!(zeros(N_adjoint), A, b; log = true)

    Λ_T = Λ[1:length(free_T)]
    Λ_z = Λ[length(free_T) + 1:end]
    dL_dα = 2 * dot(Λ_T, vec(object_opt.Tmap .- rechp.T_background)[free_T])
    dL_dβ = 2 * dot(Λ_z, vec(object_opt.zmap .- rechp.z_middle)[free_z])

    (; L, dL_dα, dL_dβ, object_opt, result, minres_hist, n_active_T, n_active_z)
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
                         minres_converged = r.minres_hist === nothing ? missing : r.minres_hist.isconverged,
                         n_active_T = r.n_active_T, n_active_z = r.n_active_z))
        verbose && @info "resolve_initial_α_β iteration $(length(history))" α β L=r.L n_active_T=r.n_active_T n_active_z=r.n_active_z
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

    (; α = α_resolved, β = β_resolved, history, return_value, object, noise, object_opt = r_final.object_opt,
       n_active_T = r_final.n_active_T, n_active_z = r_final.n_active_z)
end

# Saves resolve_initial_α_β's output into run_folder/tune_regularization/ -- a subfolder of the
# same optdata/opt_<datetime>/ run_and_save_outer_optimization_distributed (src/optimize.jl) uses
# for the main optimization, so everything from one cluster job lands together. Mirrors
# save_optimization_state's pattern (results.jld2 + history.csv + plots), plus two reconstruction
# plots: at the search's starting point (α_start, β_start) and at the resolved (α, β) --
# re-evaluating the former explicitly, the same way jupyter_notebooks/2026-10-23-tune-
# regularization.ipynb's two reconstruction sections do, since resolve_initial_α_β's own return
# value only carries the resolved reconstruction.
function save_tune_regularization_state(run_folder, tune_out, α_start, β_start, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp, PSF_scale)
    @unpack php, imghp, rechp = jhp
    @unpack α, β, history, return_value, object, noise, object_opt, n_active_T, n_active_z = tune_out

    tune_folder = joinpath(run_folder, "tune_regularization")
    mkdir(tune_folder)

    results_file = joinpath(tune_folder, "results.jld2")
    @save results_file php imghp rechp PSF_scale α β α_start β_start history return_value object noise object_opt n_active_T n_active_z
    CSV.write(joinpath(tune_folder, "history.csv"), DataFrame(history))

    plot_alpha_beta_history([h.α for h in history], [h.β for h in history];
        log_scale = true, savefile = "$(tune_folder)/alpha_beta_history.png")
    plot_objective_function_history([h.L for h in history]; savefile = "$(tune_folder)/loss_history.png")

    r_init = reconstruction_error_and_grad(object, noise, α_start, β_start, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp)
    plot_reconstruction_results(object, r_init.object_opt;
        extra_title = "α = $α_start, β = $β_start  (initial, before tuning).  n_active_T=$(r_init.n_active_T), n_active_z=$(r_init.n_active_z)",
        savefile = "$(tune_folder)/reconstruction_initial.png")

    plot_reconstruction_results(object, object_opt;
        extra_title = "α = $α, β = $β  (resolved).  n_active_T=$n_active_T, n_active_z=$n_active_z",
        savefile = "$(tune_folder)/reconstruction_resolved.png")

    results_file
end
