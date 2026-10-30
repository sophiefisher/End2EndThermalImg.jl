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
    object, noise = with_tuning_rng(jhp.php, imghp) do
        (get_object(imghp), generate_noise(imghp))
    end
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

# Top-level cache for resolve_initial_α_β results, analogous to inversedesign.jl's monochromatic-
# lens cache (find_or_design_monochromatic_lens/fields_match) -- a sibling of optdata/ and
# inversedesign/ at the project root, one subfolder per resolved (α, β), keyed by exact equality on
# the php/imghp/rechp/opthp fields that actually affect this search (TUNE_REGULARIZATION_*_FIELDS
# below), so repeated job submissions with the same setup skip re-running the sweep. This is
# distinct from run_folder/tune_regularization/ (see save_tune_regularization_state below) -- that
# is a per-job archival copy (with plots) of whichever result, fresh or cached, actually got used;
# this is the lookup table deciding whether a fresh sweep is needed at all.
#
# Requires seed_tuning_rng! (src/prepare.jl) to be called before every random draw this result
# depends on (resolve_initial_α_β's own object/noise draw, and PSF_scale's reference-object draw in
# scripts/endtoend_opt.jl) -- otherwise two runs with identical php/imghp/rechp/opthp would still
# resolve against different random realizations, and this cache would rarely hit in practice.
function ensure_tune_regularization_folder()
    root = "tune_regularization"
    isdir(root) || mkdir(root)
    root
end

# php/imghp: all raw (non-"Computed parameters") constructor fields -- the full forward/imaging
# model the sweep's reconstruction error is evaluated against (unlike the narrower monochromatic-
# lens design, nothing here is safe to drop). rechp: everything but α_init/β_init themselves (the
# thing being resolved) and the z_middle computed field. opthp: only geoms_init_type/
# geoms_init_path -- what determines the (otherwise unsaved) geoms, mirroring how inversedesign.jl
# keys on design inputs rather than the geoms array itself.
const TUNE_REGULARIZATION_PHP_FIELDS = (
    :λlb_μm, :λub_μm, :freq_order, :focal_length_μm, :num_unit_cells, :unit_cell_length_μm,
    :pillar_width_lb_μm, :pillar_width_ub_μm, :pillar_width_order, :pillar_height_μm,
    :pillar_material, :substrate_height_μm, :substrate_material, :nG
)
const TUNE_REGULARIZATION_IMGHP_FIELDS = (
    :objN, :imgN, :binN, :sampleN, :PSF_zlb_μm, :PSF_zub_μm, :PSF_zlen, :smoothness_order,
    :object_type, :noise_type
)
const TUNE_REGULARIZATION_RECHP_FIELDS = (
    :object_init_type, :T_background, :z_middle_μm, :algorithm, :maxeval, :xtol_rel, :ftol_rel, :ftol_abs
)
const TUNE_REGULARIZATION_OPTHP_FIELDS = (:geoms_init_type, :geoms_init_path)

# Scans tune_regularization/ most-recent-first (mirrors find_matching_monochromatic_design),
# returning the first subfolder whose saved php/imghp/rechp/opthp match on the relevant fields, or
# nothing if none do.
function find_matching_tune_regularization(php, imghp, rechp, opthp)
    root = ensure_tune_regularization_folder()
    for run_name in sort(readdir(root); rev = true)
        results_file = joinpath(root, run_name, "results.jld2")
        isfile(results_file) || continue
        saved = JLD2.load(results_file)
        all(haskey(saved, k) for k in ("php", "imghp", "rechp", "opthp")) || continue
        if fields_match(php, saved["php"], TUNE_REGULARIZATION_PHP_FIELDS) &&
           fields_match(imghp, saved["imghp"], TUNE_REGULARIZATION_IMGHP_FIELDS) &&
           fields_match(rechp, saved["rechp"], TUNE_REGULARIZATION_RECHP_FIELDS) &&
           fields_match(opthp, saved["opthp"], TUNE_REGULARIZATION_OPTHP_FIELDS)
            return joinpath(root, run_name)
        end
    end
    nothing
end

function save_tune_regularization_cache(php, imghp, rechp, opthp, tune_out)
    root = ensure_tune_regularization_folder()
    run_folder = joinpath(root, "resolved_$(now())")
    mkdir(run_folder)
    results_file = joinpath(run_folder, "results.jld2")
    @unpack α, β, history, return_value, object, noise, object_opt, n_active_T, n_active_z = tune_out
    @save results_file php imghp rechp opthp α β history return_value object noise object_opt n_active_T n_active_z
    run_folder
end

# Drop-in replacement for resolve_initial_α_β at call sites that want this cache -- same signature
# plus opthp up front (needed for the geoms_init_type/geoms_init_path cache key; not otherwise used
# by resolve_initial_α_β itself). On a cache hit, returns the saved result without touching geoms/
# freqs/weights/driver_plans_PSF/PSF_zcoords at all; on a miss, resolves fresh and saves it here
# before returning. Either way, the return value also carries cache_hit/cache_source (the matched
# cache folder, or nothing on a miss) -- save_tune_regularization_state below threads these through
# so a given job's own run_folder/tune_regularization/ records which case it was, not just the
# resolved α/β themselves (otherwise indistinguishable from the outside, since that folder's
# content looks the same either way).
function find_or_resolve_initial_α_β(opthp, rechp, geoms, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp::JobHyperParams, PSF_scale; kwargs...)
    @unpack php, imghp = jhp
    match = find_matching_tune_regularization(php, imghp, rechp, opthp)
    if match !== nothing
        @info "Reusing existing tune_regularization result" run_folder=match
        saved = JLD2.load(joinpath(match, "results.jld2"))
        return (; α = saved["α"], β = saved["β"], history = saved["history"], return_value = saved["return_value"],
                  object = saved["object"], noise = saved["noise"], object_opt = saved["object_opt"],
                  n_active_T = saved["n_active_T"], n_active_z = saved["n_active_z"],
                  cache_hit = true, cache_source = match)
    end
    @info "No matching tune_regularization result found; resolving a new one"
    tune_out = resolve_initial_α_β(rechp, geoms, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp, PSF_scale; kwargs...)
    save_tune_regularization_cache(php, imghp, rechp, opthp, tune_out)
    (; tune_out..., cache_hit = false, cache_source = nothing)
end

# Saves resolve_initial_α_β's output into run_folder/tune_regularization/ -- a subfolder of the
# same optdata/opt_<datetime>/ run_and_save_outer_optimization_distributed (src/optimize.jl) uses
# for the main optimization, so everything from one cluster job lands together. Mirrors
# save_optimization_state's pattern (results.jld2 + history.csv + plots), plus two reconstruction
# plots: at the search's starting point (α_start, β_start) and at the resolved (α, β) --
# re-evaluating the former explicitly, the same way jupyter_notebooks/2026-10-23-tune-
# regularization.ipynb's two reconstruction sections do, since resolve_initial_α_β's own return
# value only carries the resolved reconstruction.
#
# Called unconditionally regardless of whether tune_out came from a cache hit or a miss (same
# payload either way), so without cache_hit/cache_source this folder would look identical in both
# cases -- cache_status.txt and the reconstruction_resolved.png title make that distinction visible
# without needing to load results.jld2 to check.
function save_tune_regularization_state(run_folder, tune_out, α_start, β_start, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp, PSF_scale)
    @unpack php, imghp, rechp = jhp
    @unpack α, β, history, return_value, object, noise, object_opt, n_active_T, n_active_z, cache_hit, cache_source = tune_out

    tune_folder = joinpath(run_folder, "tune_regularization")
    mkdir(tune_folder)

    results_file = joinpath(tune_folder, "results.jld2")
    @save results_file php imghp rechp PSF_scale α β α_start β_start history return_value object noise object_opt n_active_T n_active_z cache_hit cache_source
    CSV.write(joinpath(tune_folder, "history.csv"), DataFrame(history))
    write(joinpath(tune_folder, "cache_status.txt"),
        cache_hit ? "cache hit\nsource: $cache_source\n" : "cache miss (freshly resolved)\n")

    plot_alpha_beta_history([h.α for h in history], [h.β for h in history];
        log_scale = true, savefile = "$(tune_folder)/alpha_beta_history.png")
    plot_objective_function_history([h.L for h in history]; savefile = "$(tune_folder)/loss_history.png")

    r_init = reconstruction_error_and_grad(object, noise, α_start, β_start, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp)
    plot_reconstruction_results(object, r_init.object_opt;
        extra_title = "α = $α_start, β = $β_start  (initial, before tuning).  n_active_T=$(r_init.n_active_T), n_active_z=$(r_init.n_active_z)",
        savefile = "$(tune_folder)/reconstruction_initial.png")

    cache_note = cache_hit ? "cache hit: $cache_source" : "freshly resolved"
    plot_reconstruction_results(object, object_opt;
        extra_title = "α = $α, β = $β  (resolved, $cache_note).  n_active_T=$n_active_T, n_active_z=$n_active_z",
        savefile = "$(tune_folder)/reconstruction_resolved.png")

    results_file
end
