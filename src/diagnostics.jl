function finite_difference_gradient_central(f, x; ε = 1e-6)
    n = length(x)
    grad_fd = similar(x)
    for i in 1:n
        x_plus = copy(x);  x_plus[i] += ε
        x_minus = copy(x); x_minus[i] -= ε
        grad_fd[i] = (f(x_plus) - f(x_minus)) / (2ε)
    end
    return grad_fd
end

function test_reconstruction_gradients_by_component(noisy_image, fftPSFs, weights, plan_PSF, α, β, jhp;
    ε = 1e-6, figsize_x = DEFAULT_FIGSIZE_X, figsize_y = DEFAULT_FIGSIZE_Y)
    @unpack php, imghp, rechp = jhp
    @unpack T_background, z_middle = rechp
    
    object_init_flat = flatten_object(initialize_object(imghp, rechp))
    objective_lambda = object_flat -> reconstruction_objective(object_flat, noisy_image, fftPSFs, weights, plan_PSF, α, β, T_background, z_middle, php, imghp)

    grad_autodiff = Zygote.gradient(x -> objective_lambda(x), object_init_flat)[1]
    grad_fd = finite_difference_gradient_central(objective_lambda, object_init_flat; ε = ε)

    diff_norm = norm(grad_autodiff - grad_fd) / (norm(grad_fd))
    @info "Relative gradient error for ε=$(ε): $diff_norm"

    fig, ax = subplots(figsize=(figsize_x, figsize_y))
    ax.plot(grad_fd;       label = "Finite difference", linewidth = 1, alpha = 0.8)
    ax.plot(grad_autodiff; label = "Autodiff",          linewidth = 1, alpha = 0.8)
    ax.set_xlabel("Gradient component index")
    ax.set_ylabel("Gradient value")
    ax.set_title("Gradient comparison (rel. error = $(round(diff_norm, sigdigits=3)))")
    ax.legend()

    return diff_norm, fig
end

function test_reconstruction_gradients_random_perturbation(noisy_image, fftPSFs, weights, plan_PSF, α, β, jhp;
    n_ε = 40, figsize_x = DEFAULT_FIGSIZE_X, figsize_y = DEFAULT_FIGSIZE_Y, fontsize = 10)
    @unpack php, imghp, rechp = jhp
    @unpack T_background, z_middle = rechp
    Tmap_lb = 50; Tmap_ub = 1000
    Tmap_init = Tmap_lb .+ (Tmap_ub - Tmap_lb) .* rand(imghp.objN, imghp.objN)
    zmap_init = imghp.PSF_zlb .+ (imghp.PSF_zub - imghp.PSF_zlb) .* rand(imghp.objN, imghp.objN)
    object_init_flat = flatten_object((; Tmap = Tmap_init, zmap = zmap_init))
    objective_lambda = object_flat -> reconstruction_objective(object_flat, noisy_image, fftPSFs, weights, plan_PSF, α, β, T_background, z_middle, php, imghp)

    grad_autodiff = Zygote.gradient(x -> objective_lambda(x), object_init_flat)[1]

    δ = randn(length(object_init_flat))
    δ = δ / norm(δ)  # normalize so ε controls the step size

    f0 = objective_lambda(object_init_flat)
    ε_values = 10 .^ range(-12, 0; length = n_ε)
    rel_errors = similar(ε_values)
    for (k, ε) in enumerate(ε_values)
        f_perturbed = objective_lambda(object_init_flat + ε * δ)
        # finite difference directional derivative
        fd_directional = (f_perturbed - f0) / ε
        # AD directional derivative
        ad_directional = dot(grad_autodiff, δ)
        rel_errors[k] = abs(ad_directional - fd_directional) / abs(ad_directional)
    end

    matplotlib.rcParams["font.size"] = fontsize
    fig, ax = subplots(figsize=(figsize_x, figsize_y))
    ax.loglog(ε_values, rel_errors; marker = "o", linewidth = 1, markersize = 4)
    ax.set_xlabel("Perturbation size ε")
    ax.set_ylabel("Relative directional \n derivative error")
    ax.set_title("Gradient check: relative error vs. perturbation size")
    ax.grid(true; which = "both", linestyle = "--", alpha = 0.5)

    return ε_values, rel_errors, fig
end

function dC_dz_AD(δ_Δz, PSF_zcoord, z)
    diff = PSF_zcoord - z
    dδ_ddiff = Zygote.gradient(δ_Δz, diff)[1]
    return -dδ_ddiff
end

function test_dC_dz(imghp; figsize_x = DEFAULT_FIGSIZE_X, figsize_y = DEFAULT_FIGSIZE_Y)
    @unpack smoothness_order, PSF_Δz = imghp
    δ_Δz = get_discretized_δ_function(imghp)
    PSF_zcoords = get_PSF_zcoords(imghp)

    # test on a range of z values within the support of δ_Δz
    z_test = range(imghp.PSF_zlb, imghp.PSF_zub, length = 200)

    for PSF_zcoord in PSF_zcoords
        dC_hand = [dC_dz(z, PSF_zcoord, smoothness_order, PSF_Δz) for z in z_test]
        dC_ad   = [dC_dz_AD(δ_Δz, PSF_zcoord, z) for z in z_test]

        diff_norm = norm(dC_hand - dC_ad) / norm(dC_ad)
        @info "PSF_zcoord = $PSF_zcoord: relative difference = $diff_norm"
    end

    # plot comparison for the middle PSF_zcoord
    PSF_zcoord_mid = PSF_zcoords[length(PSF_zcoords) ÷ 2]
    dC_hand = [dC_dz(z, PSF_zcoord_mid, smoothness_order, PSF_Δz) for z in z_test]
    dC_ad   = [dC_dz_AD(δ_Δz, PSF_zcoord_mid, z) for z in z_test]

    fig, ax = subplots(figsize=(figsize_x, figsize_y))
    ax.semilogy(z_test, abs.(dC_hand); label = "Hand-derived", linewidth = 1, alpha = 0.8)
    ax.semilogy(z_test, abs.(dC_ad);   label = "AD",           linewidth = 1, alpha = 0.8, linestyle = "--")
    ax.set_xlabel("z")
    ax.set_ylabel("dC/dz")
    ax.set_title("dC_dz comparison at PSF_zcoord = $PSF_zcoord_mid")
    ax.legend()

    return fig
end

function get_fftPSFs_diffable(freqs, incidents, surrogates, geoms, n2f_kernels, plans_n2f, plans_PSF,
                               php::PhysicsHyperParams, imghp::ImagingHyperParams)
    PSF_zlen = imghp.PSF_zlen
    nF = length(freqs)
    [End2EndThermalImg.get_fftPSF(End2EndThermalImg.get_PSF_at_freq_and_z(
        freqs[iF],
        incidents[iF, iZ],
        surrogates[iF],
        geoms,
        n2f_kernels[iF],
        plans_n2f[1],
        php,
        imghp
    ), plans_PSF[1]) for iF in 1:nF, iZ in 1:PSF_zlen]
end

function test_get_fftPSFs_gradient(jhp; figsize_x = DEFAULT_FIGSIZE_X, figsize_y = DEFAULT_FIGSIZE_Y)
    @unpack php, imghp, opthp = jhp
    plans_n2f, plans_PSF = get_fft_plans(php, imghp)
    geoms = initialize_geoms(php, opthp)
    PSF_zcoords = get_PSF_zcoords(imghp)
    freqs = get_freq_chebpoints(php)
    incidents = get_incident_fields(freqs, PSF_zcoords, php)
    surrogates = load_surrogate_models(php)
    n2f_kernels = get_n2f_kernels(freqs, plans_n2f[1], php, imghp)

    fftPSFs = get_fftPSFs(freqs, incidents, surrogates, geoms, n2f_kernels, plans_n2f, plans_PSF, php, imghp)

    grad_diffable = Zygote.gradient(
        g -> sum(sum.(abs2, get_fftPSFs_diffable(freqs, incidents, surrogates, g, n2f_kernels, plans_n2f, plans_PSF, php, imghp))),
        geoms
    )[1]
    grad_precomputed = Zygote.gradient(
        g -> sum(sum.(abs2, get_fftPSFs_from_precomputed(g, fftPSFs, freqs, incidents, surrogates, n2f_kernels, plans_n2f, plans_PSF, php, imghp))),
        geoms
    )[1]

    diff_norm = norm(grad_diffable - grad_precomputed) / norm(grad_diffable)
    @info "Relative gradient error between diffable and precomputed: $diff_norm"

    rel_diff = abs.(grad_diffable - grad_precomputed) ./ (abs.(grad_diffable) .+ eps())

    fig, axes = subplots(1, 3; figsize=(3 * figsize_x, 1.5 * figsize_y))

    im0 = axes[0].imshow(grad_diffable; aspect="auto")
    axes[0].set_title("Diffable gradient")
    fig.colorbar(im0, ax=axes[0])

    im1 = axes[1].imshow(grad_precomputed; aspect="auto")
    axes[1].set_title("Precomputed gradient")
    fig.colorbar(im1, ax=axes[1])

    im2 = axes[2].imshow(rel_diff; aspect="auto")
    axes[2].set_title("Relative difference per component")
    fig.colorbar(im2, ax=axes[2])

    fig.suptitle("get_fftPSFs gradient comparison (rel. error = $(round(diff_norm, sigdigits=3)))")
    fig.tight_layout()

    @info "Timing diffable gradient:"
    @btime Zygote.gradient(
        g -> sum(sum.(abs2, get_fftPSFs_diffable($freqs, $incidents, $surrogates, g, $n2f_kernels, $plans_n2f, $plans_PSF, $php, $imghp))),
        $geoms
    )
    @info "Timing precomputed gradient:"
    @btime Zygote.gradient(
        g -> sum(sum.(abs2, get_fftPSFs_from_precomputed(g, $fftPSFs, $freqs, $incidents, $surrogates, $n2f_kernels, $plans_n2f, $plans_PSF, $php, $imghp))),
        $geoms
    )

    return diff_norm, fig
end