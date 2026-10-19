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
            else
                extra_title = "[index: $(z_idx)/$(length(PSF_zcoords))]"
            end
            plot_PSFs_at_z(PSFs[:, z_idx], zcoord, freqs; savefile = "$(run_folder)/PSFs_z_$(zcoord).png", vmin = vmin, vmax = vmax, extra_title = extra_title, log_scale = log_scale)
        end
    end

    verbose && @info "Saved results" file=results_file
    # TODO: return geoms unflattened
    (; results_file, objective_opt, geoms_opt_flat, return_value, objective_history, PSF_scale)
end

