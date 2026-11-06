using End2EndThermalImg
using Distributed
using LinearAlgebra
using Statistics
using Random
using JLD2
using Dates

# Determines what freq_order (number of Chebyshev frequency points, both for the Clenshaw-Curtis
# quadrature over the wavelength band and for the surrogate model's own frequency grid) is needed
# for the forward-modeled image to converge. Moved out of a notebook (jupyter_notebooks/2026-11-05-
# test-frequency-order-convergence.ipynb) into this script because the full sweep (freq_order up
# to 80) needs far more memory than fits in a typical small interactive Jupyter allocation --
# resident incident fields alone run to several GB per worker at freq_order=80. Run via
# cluster_scripts/run_frequency_order_convergence.sh (sbatch), which gives this the same
# distributed resource shape as scripts/endtoend_opt.jl. Results are saved to
# convergence_studies/freq_order_<datetime>/results.jld2; the notebook just loads that file and
# plots -- no heavy computation happens there.
#
# To isolate only the effect of freq_order on image formation (not geometry, not object
# randomness, not depth sampling -- that's a separate, later convergence study): geometry and
# object are each fixed (computed/sampled once, reused unchanged across every freq_order);
# PSF_zlen is held fixed at the production default (3) throughout, since it's a separate numerical
# approximation (z-interpolation), not the frequency quadrature; every image is formed noiselessly
# (make_image_from_3D_distributed, not the noisy variant); PSF_scale is chosen once per geometry
# (at freq_order=20) and held fixed across the sweep for that geometry.

nw = End2EndThermalImg.setup_cluster_workers()
println("nworkers = ", nworkers(), "  workers = ", workers())
println("threads per worker = ", [remotecall_fetch(Threads.nthreads, p) for p in workers()])

# --- Fixed geometry, object, and PSF_scale (computed once at freq_order = 20) ---

@everywhere begin
    php_ref = PhysicsHyperParams(
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

    # noise_type is required by ImagingHyperParams but unused: every image below is formed with
    # make_image_from_3D_distributed (the noiseless path).
    imghp_ref = ImagingHyperParams(
        objN = 32,
        imgN = 128,
        binN = 3,
        sampleN = 1,
        PSF_zlb_μm = -2.1e6,
        PSF_zub_μm = -0.9e6,
        PSF_zlen = 3,
        smoothness_order = 3.0,
        object_type = UniformlyRandomTFixedDepthObject(Tlb = 136.15, Tub = 450.15, z_μm = -1.5e6, php = php_ref),
        noise_type = RelativeNoise(0.0),
        php = php_ref
    )
end

println("\n--- fixed geometries ---")
opthp_design = OptimizeHyperParams(geoms_init_type = "monochromatic")
@time geoms_mono = End2EndThermalImg.initialize_geoms(php_ref, opthp_design; imghp = imghp_ref)
println("extrema(geoms_mono) = ", extrema(geoms_mono), "  bounds = (", php_ref.pillar_width_lb, ", ", php_ref.pillar_width_ub, ")")

Random.seed!(2)
geoms_random = rand(php_ref.pillar_width_lb:eps():php_ref.pillar_width_ub, php_ref.num_unit_cells, php_ref.num_unit_cells)
println("extrema(geoms_random) = ", extrema(geoms_random))

println("\n--- fixed object ---")
Random.seed!(1)
object = End2EndThermalImg.get_object(imghp_ref)
println("extrema(object.Tmap) = ", extrema(object.Tmap), "  unique(object.zmap) = ", unique(object.zmap))

println("\n--- PSF_scale per geometry (at freq_order = $(php_ref.freq_order)) ---")
driver_plans_n2f_ref, driver_plans_PSF_ref = End2EndThermalImg.setup_fft_plans_distributed!(php_ref, imghp_ref)
freqs_ref = End2EndThermalImg.get_freq_chebpoints(php_ref)
n2f_kernels_ref = End2EndThermalImg.compute_and_broadcast(:n2f_kernels, End2EndThermalImg.get_n2f_kernels, freqs_ref, driver_plans_n2f_ref, php_ref, imghp_ref)

@everywhere begin
    surrogates = End2EndThermalImg.load_surrogate_models(php_ref)
    PSF_zcoords_ref = End2EndThermalImg.get_PSF_zcoords(imghp_ref)
    weights_ref = End2EndThermalImg.get_clenshaw_curtis_quadrature_weights(php_ref)
end

End2EndThermalImg.assign_shard_indices!()
End2EndThermalImg.init_resident_incidents!(freqs_ref, PSF_zcoords_ref, php_ref)

reference_scale = 1.0

End2EndThermalImg.update_geoms_and_fftPSFs!(geoms_mono, freqs_ref, php_ref, imghp_ref, reference_scale)
image_mono_at_reference = End2EndThermalImg.make_image_from_3D_distributed(object, weights_ref, php_ref, imghp_ref)
PSF_scale_mono = End2EndThermalImg.choose_PSF_scale(image_mono_at_reference, reference_scale; target_mean = 1.0)
println("PSF_scale_mono = ", PSF_scale_mono)

End2EndThermalImg.update_geoms_and_fftPSFs!(geoms_random, freqs_ref, php_ref, imghp_ref, reference_scale)
image_random_at_reference = End2EndThermalImg.make_image_from_3D_distributed(object, weights_ref, php_ref, imghp_ref)
PSF_scale_random = End2EndThermalImg.choose_PSF_scale(image_random_at_reference, reference_scale; target_mean = 1.0)
println("PSF_scale_random = ", PSF_scale_random)

# --- Sweep across freq_order ---

freq_orders = [5, 10, 20, 40, 80]
geometries = [("monochromatic", geoms_mono, PSF_scale_mono), ("random", geoms_random, PSF_scale_random)]

images = Dict{Tuple{String, Int}, Matrix{Float64}}()

for freq_order in freq_orders
    println("\n=== freq_order = $freq_order ===")

    @everywhere begin
        php = PhysicsHyperParams(
            λlb_μm = 8.0,
            λub_μm = 12.0,
            freq_order = $freq_order,
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

        imghp = ImagingHyperParams(
            objN = 32,
            imgN = 128,
            binN = 3,
            sampleN = 1,
            PSF_zlb_μm = -2.1e6,
            PSF_zub_μm = -0.9e6,
            PSF_zlen = 3,
            smoothness_order = 3.0,
            object_type = UniformlyRandomTFixedDepthObject(Tlb = 136.15, Tub = 450.15, z_μm = -1.5e6, php = php),
            noise_type = RelativeNoise(0.0),
            php = php
        )
    end

    t_setup = @elapsed begin
        driver_plans_n2f, driver_plans_PSF = End2EndThermalImg.setup_fft_plans_distributed!(php, imghp)
        freqs = End2EndThermalImg.get_freq_chebpoints(php)
        n2f_kernels = End2EndThermalImg.compute_and_broadcast(:n2f_kernels, End2EndThermalImg.get_n2f_kernels, freqs, driver_plans_n2f, php, imghp)

        @everywhere begin
            surrogates = End2EndThermalImg.load_surrogate_models(php)
            PSF_zcoords = End2EndThermalImg.get_PSF_zcoords(imghp)
            weights = End2EndThermalImg.get_clenshaw_curtis_quadrature_weights(php)
        end

        End2EndThermalImg.assign_shard_indices!()
        End2EndThermalImg.init_resident_incidents!(freqs, PSF_zcoords, php)
    end
    println("  setup: $(t_setup) s")

    for (label, geoms, PSF_scale) in geometries
        t_image = @elapsed begin
            End2EndThermalImg.update_geoms_and_fftPSFs!(geoms, freqs, php, imghp, PSF_scale)
            image = End2EndThermalImg.make_image_from_3D_distributed(object, weights, php, imghp)
        end
        images[(label, freq_order)] = image
        println("  $label: $(t_image) s  mean(image) = $(mean(image))")
    end
end

# --- Save results ---

root = "convergence_studies"
isdir(root) || mkdir(root)
datetime = now()
run_folder = joinpath(root, "freq_order_$(datetime)")
mkdir(run_folder)
results_file = joinpath(run_folder, "results.jld2")
@save results_file freq_orders geoms_mono geoms_random object PSF_scale_mono PSF_scale_random images php_ref imghp_ref

println("\n=== Done ===")
println("run_folder   = ", run_folder)
println("results_file = ", results_file)
