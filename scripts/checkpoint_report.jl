using End2EndThermalImg
using Distributed
using Logging

# Route @info/@warn to stdout instead of the default stderr -- same reasoning as endtoend_opt.jl.
global_logger(ConsoleLogger(stdout))

nw = End2EndThermalImg.setup_cluster_workers()
println("nworkers = ", nworkers(), "  workers = ", workers())
println("driver threads = ", Threads.nthreads())
println("threads per worker = ", [remotecall_fetch(Threads.nthreads, p) for p in workers()])

# Regenerates the initial-vs-checkpoint comparison artifacts (geoms_comparison_checkpoint.png,
# PSFs_checkpoint.png, reconstruction_checkpoint*.png) for a run_endtoend_opt.sh job that was
# cancelled or timed out before reaching opthp.maxeval, so it never got through the real
# final-state block in run_and_save_outer_optimization_distributed (src/optimize.jl). Everything
# needed is either already saved in that run's own results.jld2, or cheaply/deterministically
# reconstructible from it (see save_checkpoint_report's docstring, src/optimize.jl) -- this script
# is meant to be run against any optdata/opt_<datetime>/ folder with a results.jld2, finished or
# not, via:
#   E2E_RUN_FOLDER=optdata/opt_<datetime> sbatch cluster_scripts/run_checkpoint_report.sh
run_folder = get(ENV, "E2E_RUN_FOLDER", "")
isempty(run_folder) && error("E2E_RUN_FOLDER must be set to the optdata/opt_<datetime> run folder to report on")
isdir(run_folder) || error("E2E_RUN_FOLDER does not exist: $run_folder")
isfile(joinpath(run_folder, "results.jld2")) || error("$run_folder has no results.jld2 -- nothing to report on")

println("\n--- loading checkpoint from ", run_folder, " ---")
state = End2EndThermalImg.load_checkpoint_state(run_folder)
php, imghp, opthp, rechp = state.php, state.imghp, state.opthp, state.rechp
PSF_scale = state.PSF_scale
geoms_init_flat, geoms_opt_flat = state.geoms_init_flat, state.geoms_opt_flat
α_opt, β_opt = state.α_opt, state.β_opt
iteration = state.history[end].iteration
println("checkpoint at iteration ", iteration, "  α = ", α_opt, "  β = ", β_opt)

jhp = JobHyperParams(php = php, imghp = imghp, opthp = opthp, rechp = rechp)

# Mirrors endtoend_opt.jl's @everywhere begin php = PhysicsHyperParams(...); ... end block, which
# there constructs php/imghp/opthp/rechp fresh identically on every worker from literal values.
# Here they're loaded at runtime on the driver only (from results.jld2), so they must be sent to
# every worker explicitly instead -- broadcast_value is the same primitive compute_and_broadcast
# (used below for n2f_kernels) is built on (src/distributed.jl).
for (name, value) in ((:php, php), (:imghp, imghp), (:opthp, opthp), (:rechp, rechp))
    End2EndThermalImg.broadcast_value(name, value)
end

println("\n--- setup ---")
@time driver_plans_n2f, driver_plans_PSF = End2EndThermalImg.setup_fft_plans_distributed!(php, imghp)
@time freqs = End2EndThermalImg.get_freq_chebpoints(php)
@time n2f_kernels = End2EndThermalImg.compute_and_broadcast(:n2f_kernels, End2EndThermalImg.get_n2f_kernels, freqs, driver_plans_n2f, php, imghp)

@everywhere begin
    surrogates = End2EndThermalImg.load_surrogate_models(php)
    PSF_zcoords = End2EndThermalImg.get_PSF_zcoords(imghp)
    weights = End2EndThermalImg.get_clenshaw_curtis_quadrature_weights(php)
end

@time End2EndThermalImg.assign_shard_indices!()
@time End2EndThermalImg.init_resident_incidents!(freqs, PSF_zcoords, php)

# Same two human-readable extra test objects as endtoend_opt.jl, rebuilt identically (deterministic
# given php/imghp, both loaded unchanged from the original run) rather than saved/reloaded.
Tlb, Tub = End2EndThermalImg.get_object_Textent(imghp.object_type)
zlb_μm, zub_μm = End2EndThermalImg.get_object_zextent(imghp.object_type)

mit_object_type = LoadedObject(
    Tmap_path = "objdata/MIT$(imghp.objN).csv",
    zmap_path = "objdata/MIT$(imghp.objN)_zmap.csv",
    Tlb = Tlb, Tub = Tub, zlb_μm = zlb_μm, zub_μm = zub_μm,
    php = php
)
gaussian_object_type = GaussianObject(
    Tlb = Tlb, Tub = Tub, std_dev_T = imghp.objN / 4,
    zlb_μm = zlb_μm, zub_μm = zub_μm, std_dev_z = imghp.objN / 4,
    php = php
)
extra_test_objects = ["MIT" => mit_object_type, "Gaussian" => gaussian_object_type]

println("\n--- computing checkpoint report ---")
End2EndThermalImg.save_checkpoint_report(run_folder, iteration, geoms_init_flat, geoms_opt_flat, α_opt, β_opt,
    freqs, weights, driver_plans_PSF, PSF_zcoords, jhp, PSF_scale;
    extra_test_objects = extra_test_objects, verbose = true)

println("\n=== Done ===")
println("run_folder = ", run_folder)
println("iteration  = ", iteration)
