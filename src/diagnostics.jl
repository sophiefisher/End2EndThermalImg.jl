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