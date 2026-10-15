using End2EndThermalImg
using Distributed
using LinearAlgebra
using IterativeSolvers
using Zygote
using Statistics

nw = End2EndThermalImg.setup_cluster_workers()
println("nworkers = ", nworkers(), "  workers = ", workers())
println("driver threads = ", Threads.nthreads())
println("threads per worker = ", [remotecall_fetch(Threads.nthreads, p) for p in workers()])
println("driver hostname = ", gethostname())
println("hostname per worker = ", [remotecall_fetch(gethostname, p) for p in workers()])

# Same hyperparameters used throughout jupyter_notebooks/2026-09-23-develop-pipeling.ipynb and
# jupyter_notebooks/2026-10-14-add-PSF-scale.ipynb, except noise_type = RelativeNoise(0.2) (20%
# noise relative to the mean image), per this script's own request.
@everywhere begin
    php = PhysicsHyperParams(
        λlb_μm = 8.0,
        λub_μm = 12.0,
        freq_order = 20,
        focal_length_μm = 2.0e4,
        num_unit_cells = 2048,
        unit_cell_length_μm = 4.0,
        pillar_width_lb_μm = 1.8,
        pillar_width_ub_μm = 2.7,
        pillar_width_order = 160,
        pillar_height_μm = 10.0,
        pillar_material = "Si_no_absorption",
        substrate_height_μm = 300.0,
        substrate_material = "Si_no_absorption",
        nG = 1000
    )

    # no PSF_scale field (src/prepare.jl) -- chosen below via choose_PSF_scale and threaded
    # through explicitly, not stored here.
    imghp = ImagingHyperParams(
        objN = 32,
        imgN = 128,
        binN = 3,
        sampleN = 1,
        PSF_zlb_μm = -2.1e6,
        PSF_zub_μm = -0.9e6,
        PSF_zlen = 3,
        smoothness_order = 3.0,
        object_type = GaussianObject(Tlb = 136.15, Tub = 450.15, std_dev_T = 6.0, zlb_μm = -2.0e6, zub_μm = -1.0e6, std_dev_z = 6.0, php = php),
        noise_type = RelativeNoise(0.2),
        php = php
    )

    opthp = OptimizeHyperParams(
        geoms_init_type = "uniform"
    )

    rechp = ReconstructionHyperParams(
        object_init_type = "uniform_offset",
        T_background = 293.15,
        z_middle_μm = -1.5e6,
        algorithm = :LD_LBFGS,
        php = php
    )

    jhp = JobHyperParams(
        php = php,
        imghp = imghp,
        opthp = opthp,
        rechp = rechp
    )
end

α = 1e-7
β = 1e-7

# --- One-time setup, not included in any of the timings below ---

println("\n--- setup (not timed) ---")
@time End2EndThermalImg.compute_system_parameters(php, imghp)
@time driver_plans_n2f, driver_plans_PSF = End2EndThermalImg.setup_fft_plans_distributed!(php, imghp)
@time freqs = End2EndThermalImg.get_freq_chebpoints(php)
@time n2f_kernels = End2EndThermalImg.compute_and_broadcast(:n2f_kernels, End2EndThermalImg.get_n2f_kernels, freqs, driver_plans_n2f, php, imghp)

@everywhere begin
    surrogates = End2EndThermalImg.load_surrogate_models(php)
    geoms = End2EndThermalImg.initialize_geoms(php, opthp)
    PSF_zcoords = End2EndThermalImg.get_PSF_zcoords(imghp)
    weights = End2EndThermalImg.get_clenshaw_curtis_quadrature_weights(php)
end

@time End2EndThermalImg.assign_shard_indices!()
@time End2EndThermalImg.init_resident_incidents!(freqs, PSF_zcoords, php)

object = End2EndThermalImg.get_object(imghp)
noise = End2EndThermalImg.generate_noise(imghp)

# Choosing PSF_scale (so the reconstruction has a realistic scale to work with) is explicitly
# excluded from the forward/backward timing below, per request -- done once here instead. The
# second update_geoms_and_fftPSFs! call (at the chosen scale) primes the resident fftPSFs state
# before the timed section begins; the timed section's own first update_geoms_and_fftPSFs! call
# redoes this same recomputation, exactly as every outer-loop iteration would in practice.
println("\n--- choosing PSF_scale (not timed) ---")
reference_scale = 1.0
End2EndThermalImg.update_geoms_and_fftPSFs!(geoms, freqs, php, imghp, reference_scale)
image_reference = End2EndThermalImg.make_image_from_3D_distributed(object, weights, php, imghp)
PSF_scale = End2EndThermalImg.choose_PSF_scale(image_reference, reference_scale; target_mean = 1.0)
println("chosen PSF_scale = ", PSF_scale, "  (mean(image) at reference_scale=1.0 was ", mean(image_reference), ")")
End2EndThermalImg.update_geoms_and_fftPSFs!(geoms, freqs, php, imghp, PSF_scale)

# --- Timed forward + backward pass, run twice: the first run pays Julia's JIT compilation cost
# for every function in the chain; the second measures steady-state, post-compilation cost --
# i.e. what every later call inside an actual optimization loop would actually cost. Wrapped in a
# function (not left at top level) so the timed code benefits from normal Julia type
# specialization rather than the well-known slow-global-variable pitfall of top-level scripts.
function run_forward_and_backward(geoms, freqs, php, imghp, jhp, object, weights, noise, α, β,
                                   PSF_scale, driver_plans_PSF, PSF_zcoords, pass_label)
    println("\n=== $pass_label ===")

    t_update_geoms = @elapsed End2EndThermalImg.update_geoms_and_fftPSFs!(geoms, freqs, php, imghp, PSF_scale)
    println("[$pass_label] update_geoms_and_fftPSFs!            : $(t_update_geoms) s")

    t_noisy_image = @elapsed noisy_image = End2EndThermalImg.make_noisy_image_from_3D_distributed(
        object, weights, noise, imghp.noise_type, php, imghp)
    println("[$pass_label] make_noisy_image_from_3D_distributed : $(t_noisy_image) s")

    t_reconstruct = @elapsed result = End2EndThermalImg.reconstruct_Tmap_and_zmap_distributed(
        noisy_image, weights, α, β, jhp;
        maxeval = 2000, xtol_rel = 0, ftol_rel = 0.0, ftol_abs = 0.0, objective_scale = :auto)
    println("[$pass_label] reconstruct_Tmap_and_zmap_distributed: $(t_reconstruct) s",
            "  (return_value=", result.return_value, ", iterations=", length(result.objective_history), ")")
    object_opt = result.object_opt

    forward_total = t_update_geoms + t_noisy_image + t_reconstruct
    println("[$pass_label] forward model total                  : $(forward_total) s")

    t_gather = @elapsed fftPSFs = let
        fftPSFs = Matrix{Matrix{ComplexF64}}(undef, length(freqs), length(PSF_zcoords))
        for p in workers()
            grid_inds, shard = remotecall_fetch(p) do
                (End2EndThermalImg._shard_grid_inds, End2EndThermalImg._fftPSFs_shard)
            end
            for idx in grid_inds
                iF, iZ = Tuple(idx)
                fftPSFs[iF, iZ] = shard[(iF, iZ)]
            end
        end
        fftPSFs
    end
    println("[$pass_label] gather fftPSFs                       : $(t_gather) s")

    b = N_adjoint = nothing
    t_build_b = @elapsed begin
        _, active_linear_indices = End2EndThermalImg.get_active_variables(object_opt.zmap, imghp)
        free_linear_indices = setdiff(1:imghp.objN^2, active_linear_indices)
        b_T = vec(2 .* (object.Tmap .- object_opt.Tmap) ./ norm(object.Tmap)^2)
        b_z = vec(2 .* (object.zmap .- object_opt.zmap) ./ norm(object.zmap)^2)[free_linear_indices]
        b = [b_T; b_z]
        N_adjoint = 2 * imghp.objN^2 - length(active_linear_indices)
    end
    println("[$pass_label] build b                               : $(t_build_b) s")

    t_hessian = @elapsed A = End2EndThermalImg.build_hessian_operator_distributed(
        object, object_opt, fftPSFs, freqs, weights, driver_plans_PSF, noise, imghp.noise_type, α, β, php, imghp)
    println("[$pass_label] build_hessian_operator_distributed    : $(t_hessian) s")

    t_minres = @elapsed Λ, minres_hist = IterativeSolvers.minres!(zeros(N_adjoint), A, b; log = true)
    println("[$pass_label] minres! (adjoint solve)               : $(t_minres) s  (", minres_hist, ")")

    t_grad = @elapsed grad = Zygote.gradient(
        g -> End2EndThermalImg.Λ_dot_dFdo_distributed(
            Λ, object, object_opt, g, fftPSFs, freqs, weights, driver_plans_PSF, noise, imghp.noise_type, php, imghp, PSF_scale),
        geoms
    )[1]
    println("[$pass_label] Zygote.gradient(Λ_dot_dFdo_distributed): $(t_grad) s")

    backward_total = t_gather + t_build_b + t_hessian + t_minres + t_grad
    println("[$pass_label] backward model total                  : $(backward_total) s")

    total = forward_total + backward_total
    println("[$pass_label] TOTAL (forward + backward)            : $(total) s")

    (; forward_total, backward_total, total)
end

results_pass1 = run_forward_and_backward(geoms, freqs, php, imghp, jhp, object, weights, noise, α, β,
                                          PSF_scale, driver_plans_PSF, PSF_zcoords, "Pass 1 (includes compilation)")
results_pass2 = run_forward_and_backward(geoms, freqs, php, imghp, jhp, object, weights, noise, α, β,
                                          PSF_scale, driver_plans_PSF, PSF_zcoords, "Pass 2 (compilation-free)")

println("\n=== Summary (compilation-free, Pass 2) ===")
println("Forward model time : ", results_pass2.forward_total, " s")
println("Backward model time: ", results_pass2.backward_total, " s")
println("Total time         : ", results_pass2.total, " s")
