function ensure_inversedesign_folder()
    folder = "inversedesign"
    if !isdir(folder)
        mkdir(folder)
        @info "Created folder: $folder"
    else
        @info "Folder already exists: $folder"
    end
    return folder
end

# Loads geoms saved by design_monochromatic_lens (opthp.geoms_init_type == "loaded", with
# opthp.geoms_init_path pointing at a run folder, e.g. "inversedesign/monochromatic_lens_<datetime>").
# The size check guards against loading a design made with a different num_unit_cells, which would
# otherwise only surface later as a confusing dimension mismatch somewhere downstream.
function load_geoms(php::PhysicsHyperParams, path::String)
    geoms_opt_flat = JLD2.load(joinpath(path, "results.jld2"), "geoms_opt_flat")
    geoms = unflatten_square_matrix(geoms_opt_flat)
    size(geoms) == (php.num_unit_cells, php.num_unit_cells) ||
        error("Loaded geoms size $(size(geoms)) doesn't match php.num_unit_cells = $(php.num_unit_cells)")
    geoms
end

# The php/imghp fields design_monochromatic_lens actually depends on (traced through its full call
# chain: surrogate loading/filename, frequency and PSF-depth grids, PSF pixel binning, NLopt
# bounds) -- everything else (e.g. imghp.noise_type/object_type/smoothness_order,
# php.substrate_height_μm, opthp entirely) can differ without invalidating a cached design. Only
# the raw/user-facing fields are listed; php/imghp's "computed" fields (wavcen, pillar_width_lb,
# PSF_zlb, etc.) are deterministic functions of these, so comparing the raw ones is sufficient.
const MONOCHROMATIC_DESIGN_PHP_FIELDS = (
    :λlb_μm, :λub_μm, :freq_order, :focal_length_μm, :num_unit_cells, :unit_cell_length_μm,
    :pillar_width_lb_μm, :pillar_width_ub_μm, :pillar_width_order, :pillar_height_μm,
    :pillar_material, :substrate_material, :nG
)
const MONOCHROMATIC_DESIGN_IMGHP_FIELDS = (:objN, :imgN, :binN, :sampleN, :PSF_zlb_μm, :PSF_zub_μm, :PSF_zlen)

fields_match(a, b, fields) = all(getfield(a, f) == getfield(b, f) for f in fields)

# Scans inversedesign/ for an existing design_monochromatic_lens run whose saved php/imghp match
# the current ones on the fields that function actually depends on. Folders are named
# monochromatic_lens_<datetime>, so sorting the names descending gives most-recent-first; returns
# the first (i.e. most recent) match, or nothing if there isn't one. Anything that isn't a
# results.jld2-containing folder with saved php/imghp (there's currently nothing else under
# inversedesign/, but just in case) is skipped rather than erroring.
function find_matching_monochromatic_design(php, imghp)
    root = ensure_inversedesign_folder()
    for run_name in sort(readdir(root); rev = true)
        run_folder = joinpath(root, run_name)
        results_file = joinpath(run_folder, "results.jld2")
        isfile(results_file) || continue
        saved = JLD2.load(results_file)
        haskey(saved, "php") && haskey(saved, "imghp") || continue
        if fields_match(php, saved["php"], MONOCHROMATIC_DESIGN_PHP_FIELDS) &&
           fields_match(imghp, saved["imghp"], MONOCHROMATIC_DESIGN_IMGHP_FIELDS)
            return run_folder
        end
    end
    nothing
end

# geoms_init_type == "monochromatic": reuses an existing design_monochromatic_lens run that matches
# (see find_matching_monochromatic_design), or designs and saves a fresh one (with
# design_monochromatic_lens's own defaults) if there isn't one yet -- a self-populating cache under
# inversedesign/, shared across any job that initializes this way.
function find_or_design_monochromatic_lens(php, imghp)
    match = find_matching_monochromatic_design(php, imghp)
    if match !== nothing
        @info "Reusing existing monochromatic lens design" run_folder=match
        return load_geoms(php, match)
    end
    @info "No matching monochromatic lens design found; designing a new one"
    out = design_monochromatic_lens(php, imghp)
    unflatten_square_matrix(out.geoms_opt_flat)
end

# TODO: add verbose option; apply to wrap objective function and @info statements
function design_monochromatic_lens(php, imghp; geoms_init_type = "uniform", xtol_rel = 1e-8, maxeval = 5000, iteration_print = 50, plot_PSFs = false, log_scale = true, verbose = true, target_objective = 1.0)
    datetime = now()
    verbose && @info "Designing monochromatic lens" datetime=datetime
    root = ensure_inversedesign_folder()
    run_name = "monochromatic_lens_$(datetime)"
    run_folder = joinpath(root, run_name)
    mkdir(run_folder)

    # assume center zcoord
    PSF_zcoords = End2EndThermalImg.get_PSF_zcoords(imghp)
    center_zcoord_idx = cld(length(PSF_zcoords), 2)
    center_zcoord = PSF_zcoords[center_zcoord_idx]
    
    # assume center frequency
    freqs = End2EndThermalImg.get_freq_chebpoints(php)
    center_freq_idx = cld(length(freqs), 2)
    center_freq = freqs[center_freq_idx]

    surrogate = load_surrogate_model_at_freq(php, center_freq_idx)
    incident = get_incident_field(center_freq, center_zcoord, php)
    plans_n2f, _ = get_fft_plans(php, imghp)
    n2f_kernel = get_n2f_kernel(center_freq, plans_n2f[1], php, imghp)

    function objective_at_scale(geoms_flat, PSF_scale)
        geoms = unflatten_square_matrix(geoms_flat)
        PSF = get_PSF_at_freq_and_z(center_freq, incident, surrogate, geoms, n2f_kernel, plans_n2f[1], php, imghp, PSF_scale)
        center_PSF_idx = cld(size(PSF, 1), 2)
        PSF[center_PSF_idx, center_PSF_idx]
    end

    objective_history = Float64[]
    geoms_init_flat = initialize_geoms(php::PhysicsHyperParams, geoms_init_type)[:]

    # PSF_scale enters the forward model purely multiplicatively (far_field_to_PSF's last step),
    # so objective_at_scale is exactly linear in it -- auto-choose it (same technique as
    # choose_PSF_scale) so the very first evaluation, at geoms_init_flat, equals target_objective,
    # rather than requiring the caller to guess a reasonable absolute scale by hand.
    reference_scale = 1.0
    objective_at_reference_scale = objective_at_scale(geoms_init_flat, reference_scale)
    PSF_scale = choose_PSF_scale(objective_at_reference_scale, reference_scale; target_mean = target_objective)
    verbose && @info "Auto-chose PSF_scale" PSF_scale objective_at_reference_scale
    objective(geoms_flat) = objective_at_scale(geoms_flat, PSF_scale)

    opt = Opt(:LD_MMA, php.num_unit_cells^2) # TODO: check algorithm choice
    opt.lower_bounds = fill(php.pillar_width_lb, php.num_unit_cells^2)
    opt.upper_bounds = fill(php.pillar_width_ub, php.num_unit_cells^2)
    opt.max_objective = (x, grad) -> nlopt_wrap_objective_autodiff(x, grad, objective, objective_history; iteration_print = iteration_print, verbose = verbose)
    opt.xtol_rel = xtol_rel
    opt.maxeval = maxeval

    (objective_opt, geoms_opt_flat, return_value) = NLopt.optimize(opt, geoms_init_flat)
    verbose && @info "Done designing monochromatic lens"
    verbose && @info "Optimization results" objective_opt return_value

    results_file = joinpath(run_folder, "results.jld2")
    @save results_file objective_opt geoms_opt_flat return_value objective_history center_freq center_zcoord datetime php imghp PSF_scale
    plot_geoms(geoms_opt_flat, php; savefile = "$(run_folder)/geoms.png")
    plot_objective_function_history(objective_history; savefile = "$(run_folder)/objective_history.png")
    if plot_PSFs
        surrogates = load_surrogate_models(php)
        incidents = get_incident_fields(freqs, PSF_zcoords, php)
        n2f_kernels = get_n2f_kernels(freqs, plans_n2f, php, imghp)
        PSFs = get_PSFs(freqs, incidents, surrogates, unflatten_square_matrix(geoms_opt_flat), n2f_kernels, plans_n2f, php, imghp, PSF_scale)
        vmin = minimum([minimum(PSF) for PSF in PSFs])
        vmax = maximum([maximum(PSF) for PSF in PSFs])
        for (z_idx, zcoord) in enumerate(PSF_zcoords)
            if z_idx == center_zcoord_idx
                extra_title = "[Design depth. index: $(z_idx)/$(length(PSF_zcoords))]"
                center_tag = "_center"
            else
                extra_title = "[index: $(z_idx)/$(length(PSF_zcoords))]"
                center_tag = ""
            end
            plot_PSFs_at_z(PSFs[:, z_idx], zcoord, freqs; savefile = "$(run_folder)/PSFs_z_$(zcoord)$(center_tag).png", vmin = vmin, vmax = vmax, extra_title = extra_title, log_scale = log_scale)
        end
    end

    verbose && @info "Saved results" file=results_file
    # TODO: return geoms unflattened
    (; results_file, objective_opt, geoms_opt_flat, return_value, objective_history, PSF_scale)
end

