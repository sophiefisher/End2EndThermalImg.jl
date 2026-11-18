using End2EndThermalImg
using Distributed
using Logging
using JLD2
using Dates

# Route @info/@warn to stdout -- same reasoning as endtoend_opt.jl/checkpoint_report.jl.
global_logger(ConsoleLogger(stdout))

nw = End2EndThermalImg.setup_cluster_workers()
println("nworkers = ", nworkers(), "  workers = ", workers())
println("driver threads = ", Threads.nthreads())
println("threads per worker = ", [remotecall_fetch(Threads.nthreads, p) for p in workers()])

# Tests whether rechp.xtol_rel can let the inner reconstruction (reconstruct_Tmap_and_zmap_distributed)
# stop well short of rechp.maxeval = 2000 without changing the reconstruction error -- currently
# rechp.xtol_rel = rechp.ftol_rel = rechp.ftol_abs = 0.0, so maxeval is the only thing that ever
# stops it (confirmed by MAXEVAL_REACHED on nearly every outer iteration of the real run), which is
# a large part of why each outer iteration costs ~30 minutes. See the planning discussion for why
# xtol_rel specifically (not ftol_rel/ftol_abs, which have a documented history of causing spurious
# premature :SUCCESS in this exact pipeline -- see the comment above
# reconstruct_Tmap_and_zmap_distributed, src/forward_distributed.jl).
#
# Moved out of a notebook for the same reason scripts/psf_zlen_convergence.jl/
# frequency_order_convergence.jl were: this needs real distributed workers (reconstruct_Tmap_and_zmap_distributed's
# dominant cost, make_image_from_3D_distributed, pmaps across workers() every inner NLopt
# iteration), which the standing interactive Jupyter session's own SLURM allocation is too small to
# provide (NumNodes=1, NumCPUs=2, NumTasks=1 -- SlurmManager has no spare task slot to launch a
# worker into, so setup_cluster_workers() hangs until launch_timeout). Run via:
#   sbatch cluster_scripts/run_test_reconstruction_xtol.sh
# Results are saved to convergence_studies/xtol_rel_<datetime>/results.jld2; the notebook
# (jupyter_notebooks/2026-11-18-test-reconstruction-maxeval.ipynb) just loads that file and plots.
#
# source_run_folder: the optdata/opt_<datetime>/ job this test's "initial" and "checkpoint"
# geoms/α/β come from -- override via E2E_SOURCE_RUN_FOLDER if testing against a different run.
source_run_folder = get(ENV, "E2E_SOURCE_RUN_FOLDER", "optdata/opt_2026-11-13T13:07:16.666")
isdir(source_run_folder) || error("E2E_SOURCE_RUN_FOLDER does not exist: $source_run_folder")

println("\n--- loading checkpoint from ", source_run_folder, " ---")
state = End2EndThermalImg.load_checkpoint_state(source_run_folder)
php, imghp, opthp, rechp = state.php, state.imghp, state.opthp, state.rechp
PSF_scale = state.PSF_scale
geoms_init_flat, geoms_opt_flat = state.geoms_init_flat, state.geoms_opt_flat
α_init, β_init = state.history[1].α, state.history[1].β   # resolved initial α/β, pre-training
α_opt, β_opt = state.α_opt, state.β_opt                    # α/β as of the checkpoint
checkpoint_iteration = state.history[end].iteration
println("checkpoint iteration = ", checkpoint_iteration)
println("α_init = ", α_init, "  β_init = ", β_init)
println("α_opt  = ", α_opt,  "  β_opt  = ", β_opt)

jhp = JobHyperParams(php = php, imghp = imghp, opthp = opthp, rechp = rechp)

# php/imghp/opthp/rechp are loaded at runtime on the driver only, so they need to be sent to every
# worker explicitly -- mirrors scripts/checkpoint_report.jl.
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

# Plain (unseeded) draw, not the deterministic with_tuning_rng one used for
# reconstruction_initial*.png. Reused as-is for every comparison below (both geoms states, every
# xtol_rel value), so xtol_rel (and separately, geoms state) is the only thing actually varying.
object = End2EndThermalImg.get_object(imghp)
noise = End2EndThermalImg.generate_noise(imghp)

# ReconstructionHyperParams (src/prepare.jl) is a plain struct with a hand-written keyword
# constructor, not @with_kw, so there's no built-in "clone with override" syntax -- this just
# rebuilds it, passing every field through unchanged except xtol_rel/maxeval.
with_xtol_rel(rechp, xtol_rel; maxeval = rechp.maxeval) = ReconstructionHyperParams(
    object_init_type = rechp.object_init_type,
    T_background = rechp.T_background,
    z_middle_μm = rechp.z_middle_μm,
    algorithm = rechp.algorithm,
    α_init = rechp.α_init,
    β_init = rechp.β_init,
    maxeval = maxeval,
    xtol_rel = xtol_rel,
    ftol_rel = rechp.ftol_rel,
    ftol_abs = rechp.ftol_abs,
    php = php
)

# ftol_rel/ftol_abs stay at 0.0 throughout -- only xtol_rel is being tested. maxeval_ceiling is a
# safety net, not the thing driving convergence here: if some xtol_rel value still hits it
# (MAXEVAL_REACHED instead of a tolerance-triggered return), that value was too tight to matter at
# this budget.
states = [
    ("initial",    geoms_init_flat, α_init, β_init),
    ("checkpoint", geoms_opt_flat,  α_opt,  β_opt),
]
xtol_rels = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
maxeval_ceiling = 2000

println("\n--- running sweep ---")
summary_rows = NamedTuple[]
object_opts = Dict{Tuple{String,Float64}, Any}()
image_opts = Dict{Tuple{String,Float64}, Any}()
noisy_images = Dict{String, Any}()

for (state_name, geoms_flat, α, β) in states
    geoms = End2EndThermalImg.unflatten_square_matrix(geoms_flat)
    End2EndThermalImg.gather_PSFs_distributed(geoms, freqs, php, imghp, PSF_scale, driver_plans_PSF, PSF_zcoords)
    noisy_image = End2EndThermalImg.make_noisy_image_from_3D_distributed(object, weights, noise, imghp.noise_type, php, imghp)
    noisy_images[state_name] = noisy_image

    for tol in xtol_rels
        rechp_test = with_xtol_rel(rechp, tol; maxeval = maxeval_ceiling)
        jhp_test = JobHyperParams(jhp; rechp = rechp_test)
        result = nothing
        t = @elapsed result = End2EndThermalImg.reconstruct_Tmap_and_zmap_distributed(noisy_image, weights, α, β, jhp_test; objective_scale = :auto)
        L = End2EndThermalImg.L_reconstruction_error(object, result.object_opt)
        n_iterations = length(result.objective_history)
        image_opt = End2EndThermalImg.make_image_from_3D_distributed(result.object_opt, weights, php, imghp)

        push!(summary_rows, (; state = state_name, xtol_rel = tol, return_value = result.return_value,
                              n_iterations = n_iterations, L = L, time_s = t))
        object_opts[(state_name, tol)] = result.object_opt
        image_opts[(state_name, tol)] = image_opt
        println("$state_name  xtol_rel=$tol  -> $(result.return_value) after $n_iterations iters, " *
                "L=$(round(L, sigdigits=4)), $(round(t, digits=1))s")
    end
end

# --- Save results ---

root = "convergence_studies"
isdir(root) || mkdir(root)
datetime = now()
run_folder = joinpath(root, "xtol_rel_$(datetime)")
mkdir(run_folder)
results_file = joinpath(run_folder, "results.jld2")
@save results_file source_run_folder checkpoint_iteration xtol_rels maxeval_ceiling object noise summary_rows object_opts image_opts noisy_images php imghp

println("\n=== Done ===")
println("run_folder   = ", run_folder)
println("results_file = ", results_file)
