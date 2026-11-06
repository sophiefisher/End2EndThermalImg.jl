using End2EndThermalImg
using Distributed
using LinearAlgebra
using Statistics
using Random
using JLD2
using Dates

# Determines what imghp.PSF_zlen (number of z-nodes, evenly spaced across PSF_zlb_μm/PSF_zub_μm,
# at which exact PSFs are computed -- any object pixel's true depth between two nodes is
# interpolated, see get_C_interp_3D) is needed for the forward-modeled image to converge. Mirrors
# scripts/frequency_order_convergence.jl's structure/rationale (moved out of a notebook for the
# same memory reason -- see that script's own header comment), run via
# cluster_scripts/run_psf_zlen_convergence.sh (sbatch). Results are saved to
# convergence_studies/psf_zlen_<datetime>/results.jld2; the notebook just loads that file and
# plots.
#
# Unlike freq_order, PSF_zlen is unrelated to any object_type's own z-sampling (e.g.
# UniformlyRandomObject.zlen is just how finely the *synthetic test object's own ground truth* is
# drawn, not a forward-model numerical-accuracy parameter) -- what actually matters here is how
# rapidly the true PSF varies with z across the fixed PSF_zlb_μm/PSF_zub_μm range. So, unlike the
# freq_order study (which deliberately used a single fixed object depth to avoid confounding),
# this one deliberately uses an object whose pixels span a *range* of depths (UniformlyRandomObject,
# not the fixed-depth variant), so a single formed image exercises interpolation quality across
# many different z-values at once. Its own zlen is set large (1000) so its z-grid is effectively
# continuous and isn't itself a limiting factor.
#
# freq_order is held fixed at 40 throughout (the production value recommended by the freq_order
# convergence study), so frequency-quadrature error doesn't leak into this comparison. Since php
# (and therefore freqs/surrogates/n2f_kernels/FFT plans, none of which depend on PSF_zlen) stays
# identical across the whole sweep, those are only computed once -- only PSF_zcoords and the
# resident incident fields (which do depend on PSF_zlen) get rebuilt per iteration.

nw = End2EndThermalImg.setup_cluster_workers()
println("nworkers = ", nworkers(), "  workers = ", workers())
println("threads per worker = ", [remotecall_fetch(Threads.nthreads, p) for p in workers()])

# --- Fixed php (freq_order = 40 throughout) and reference imghp (PSF_zlen = 3, for the
# monochromatic-lens design cache / PSF_scale reference) ---

@everywhere begin
    php = PhysicsHyperParams(
        λlb_μm = 8.0,
        λub_μm = 12.0,
        freq_order = 40,
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
        object_type = UniformlyRandomObject(Tlb = 136.15, Tub = 450.15, zlb_μm = -2.0e6, zub_μm = -1.0e6, zlen = 1000, php = php),
        noise_type = RelativeNoise(0.0),
        php = php
    )
end

println("\n--- fixed geometries ---")
opthp_design = OptimizeHyperParams(geoms_init_type = "monochromatic")
@time geoms_mono = End2EndThermalImg.initialize_geoms(php, opthp_design; imghp = imghp_ref)
println("extrema(geoms_mono) = ", extrema(geoms_mono), "  bounds = (", php.pillar_width_lb, ", ", php.pillar_width_ub, ")")

Random.seed!(2)
geoms_random = rand(php.pillar_width_lb:eps():php.pillar_width_ub, php.num_unit_cells, php.num_unit_cells)
println("extrema(geoms_random) = ", extrema(geoms_random))

println("\n--- fixed object ---")
Random.seed!(1)
object = End2EndThermalImg.get_object(imghp_ref)
println("extrema(object.Tmap) = ", extrema(object.Tmap), "  extrema(object.zmap) = ", extrema(object.zmap))

println("\n--- setup (freq_order = $(php.freq_order), shared across every PSF_zlen below) ---")
driver_plans_n2f, driver_plans_PSF = End2EndThermalImg.setup_fft_plans_distributed!(php, imghp_ref)
freqs = End2EndThermalImg.get_freq_chebpoints(php)
n2f_kernels = End2EndThermalImg.compute_and_broadcast(:n2f_kernels, End2EndThermalImg.get_n2f_kernels, freqs, driver_plans_n2f, php, imghp_ref)

@everywhere begin
    surrogates = End2EndThermalImg.load_surrogate_models(php)
    weights = End2EndThermalImg.get_clenshaw_curtis_quadrature_weights(php)
end

println("\n--- PSF_scale per geometry (at PSF_zlen = $(imghp_ref.PSF_zlen)) ---")
@everywhere begin
    imghp = imghp_ref
    PSF_zcoords = End2EndThermalImg.get_PSF_zcoords(imghp)
end

End2EndThermalImg.assign_shard_indices!()
End2EndThermalImg.init_resident_incidents!(freqs, PSF_zcoords, php)

reference_scale = 1.0

End2EndThermalImg.update_geoms_and_fftPSFs!(geoms_mono, freqs, php, imghp_ref, reference_scale)
image_mono_at_reference = End2EndThermalImg.make_image_from_3D_distributed(object, weights, php, imghp_ref)
PSF_scale_mono = End2EndThermalImg.choose_PSF_scale(image_mono_at_reference, reference_scale; target_mean = 1.0)
println("PSF_scale_mono = ", PSF_scale_mono)

End2EndThermalImg.update_geoms_and_fftPSFs!(geoms_random, freqs, php, imghp_ref, reference_scale)
image_random_at_reference = End2EndThermalImg.make_image_from_3D_distributed(object, weights, php, imghp_ref)
PSF_scale_random = End2EndThermalImg.choose_PSF_scale(image_random_at_reference, reference_scale; target_mean = 1.0)
println("PSF_scale_random = ", PSF_scale_random)

# --- Sweep across PSF_zlen ---

psf_zlens = [3, 5, 9, 17, 33]
geometries = [("monochromatic", geoms_mono, PSF_scale_mono), ("random", geoms_random, PSF_scale_random)]

images = Dict{Tuple{String, Int}, Matrix{Float64}}()

for psf_zlen in psf_zlens
    println("\n=== PSF_zlen = $psf_zlen ===")

    @everywhere begin
        imghp = ImagingHyperParams(
            objN = 32,
            imgN = 128,
            binN = 3,
            sampleN = 1,
            PSF_zlb_μm = -2.1e6,
            PSF_zub_μm = -0.9e6,
            PSF_zlen = $psf_zlen,
            smoothness_order = 3.0,
            object_type = UniformlyRandomObject(Tlb = 136.15, Tub = 450.15, zlb_μm = -2.0e6, zub_μm = -1.0e6, zlen = 1000, php = php),
            noise_type = RelativeNoise(0.0),
            php = php
        )
        PSF_zcoords = End2EndThermalImg.get_PSF_zcoords(imghp)
    end

    t_setup = @elapsed begin
        End2EndThermalImg.assign_shard_indices!()
        End2EndThermalImg.init_resident_incidents!(freqs, PSF_zcoords, php)
    end
    println("  setup: $(t_setup) s")

    for (label, geoms, PSF_scale) in geometries
        t_image = @elapsed begin
            End2EndThermalImg.update_geoms_and_fftPSFs!(geoms, freqs, php, imghp, PSF_scale)
            image = End2EndThermalImg.make_image_from_3D_distributed(object, weights, php, imghp)
        end
        images[(label, psf_zlen)] = image
        println("  $label: $(t_image) s  mean(image) = $(mean(image))")
    end
end

# --- Save results ---

root = "convergence_studies"
isdir(root) || mkdir(root)
datetime = now()
run_folder = joinpath(root, "psf_zlen_$(datetime)")
mkdir(run_folder)
results_file = joinpath(run_folder, "results.jld2")
@save results_file psf_zlens geoms_mono geoms_random object PSF_scale_mono PSF_scale_random images php imghp_ref

println("\n=== Done ===")
println("run_folder   = ", run_folder)
println("results_file = ", results_file)
