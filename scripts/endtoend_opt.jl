using End2EndThermalImg
using Distributed
using LinearAlgebra
using Statistics

nw = End2EndThermalImg.setup_cluster_workers()
println("nworkers = ", nworkers(), "  workers = ", workers())
println("driver threads = ", Threads.nthreads())
println("threads per worker = ", [remotecall_fetch(Threads.nthreads, p) for p in workers()])

# Same hyperparameters used throughout jupyter_notebooks/2026-10-19-run-optimization.ipynb and
# scripts/time_pipeline.jl.
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
        geoms_init_type = "uniform",
        maxeval = 10,          # outer (geoms, α, β) iterations -- small test run for now
        learning_rate = 1e-3
    )

    rechp = ReconstructionHyperParams(
        object_init_type = "uniform_offset",
        T_background = 293.15,
        z_middle_μm = -1.5e6,
        algorithm = :LD_LBFGS,
        α_init = 1e-7,
        β_init = 1e-7,
        maxeval = 2000,
        xtol_rel = 0.0,
        ftol_rel = 0.0,
        ftol_abs = 0.0,
        php = php
    )

    jhp = JobHyperParams(
        php = php,
        imghp = imghp,
        opthp = opthp,
        rechp = rechp
    )
end

println("\n--- setup ---")
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

# Choosing PSF_scale: form a true object's image at a reference scale, then rescale so
# mean(image) = 1.0 -- same approach as jupyter_notebooks/2026-10-14-add-PSF-scale.ipynb. Held
# fixed for the whole run (not re-chosen per outer iteration).
println("\n--- choosing PSF_scale ---")
reference_scale = 1.0
End2EndThermalImg.update_geoms_and_fftPSFs!(geoms, freqs, php, imghp, reference_scale)
object_reference = End2EndThermalImg.get_object(imghp)
image_reference = End2EndThermalImg.make_image_from_3D_distributed(object_reference, weights, php, imghp)
PSF_scale = End2EndThermalImg.choose_PSF_scale(image_reference, reference_scale; target_mean = 1.0)
println("chosen PSF_scale = ", PSF_scale, "  (mean(image) at reference_scale=1.0 was ", mean(image_reference), ")")
End2EndThermalImg.update_geoms_and_fftPSFs!(geoms, freqs, php, imghp, PSF_scale)

println("\n--- running outer optimization (", opthp.maxeval, " iterations) ---")
out = End2EndThermalImg.run_and_save_outer_optimization_distributed(
    geoms, rechp.α_init, rechp.β_init, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp, PSF_scale;
    verbose = true
)

println("\n=== Done ===")
println("run_folder   = ", out.run_folder)
println("results_file = ", out.results_file)
println("final α      = ", out.α)
println("final β      = ", out.β)
println("final L      = ", out.history[end].L)
