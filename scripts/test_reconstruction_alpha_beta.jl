using End2EndThermalImg
using Distributed
using Logging
using JLD2
using Dates

# Route @info/@warn to stdout -- same reasoning as the other scripts/*.jl cluster jobs.
global_logger(ConsoleLogger(stdout))

nw = End2EndThermalImg.setup_cluster_workers()
println("nworkers = ", nworkers(), "  workers = ", workers())
println("driver threads = ", Threads.nthreads())
println("threads per worker = ", [remotecall_fetch(Threads.nthreads, p) for p in workers()])

# Diagnoses the real run's observation that the reconstructed depth map (zmap) converges to
# z_middle (the regularization target) both at the initial state and at the checkpoint, with β
# barely moving over 150 training iterations (~1.0e-7 throughout) -- consistent with
# resolve_initial_α_β (a local, unconstrained, single-starting-point LD_LBFGS search, see
# src/tune_regularization.jl) getting stuck near its own starting point (α_start=β_start=1e-7)
# rather than having genuinely explored the landscape.
#
# Sweeps a full (α, β) grid against the monochromatic-lens INITIAL geoms only (not the
# checkpoint's optimized geoms -- "for now", per the planning discussion), holding geoms and the
# test object/noise fixed throughout (only α/β vary), so gather_PSFs_distributed and the noisy
# image only need to be computed once. Moved out of a notebook for the same reason
# scripts/test_reconstruction_xtol.jl was: needs real distributed workers.
#
# Uses xtol_rel = 1e-5 (validated in convergence_studies/xtol_rel_2026-11-18T17:30:44.528 to match
# xtol_rel=0/maxeval=2000's L to 3-4 sig figs at ~3x fewer iterations) rather than the source run's
# own saved rechp (which predates that change), to keep this 64-reconstruction sweep affordable.
#
# Run via: sbatch cluster_scripts/run_test_reconstruction_alpha_beta.sh
source_run_folder = get(ENV, "E2E_SOURCE_RUN_FOLDER", "optdata/opt_2026-11-13T13:07:16.666")
isdir(source_run_folder) || error("E2E_SOURCE_RUN_FOLDER does not exist: $source_run_folder")

println("\n--- loading checkpoint from ", source_run_folder, " ---")
state = End2EndThermalImg.load_checkpoint_state(source_run_folder)
php, imghp, opthp, rechp = state.php, state.imghp, state.opthp, state.rechp
PSF_scale = state.PSF_scale
geoms_init_flat = state.geoms_init_flat

jhp = JobHyperParams(php = php, imghp = imghp, opthp = opthp, rechp = rechp)

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

println("\n--- forming noisy image (monochromatic-lens initial geoms, fixed throughout) ---")
geoms_init = End2EndThermalImg.unflatten_square_matrix(geoms_init_flat)
End2EndThermalImg.gather_PSFs_distributed(geoms_init, freqs, php, imghp, PSF_scale, driver_plans_PSF, PSF_zcoords)

# Same deterministic draw reconstruction_initial.png already used (with_tuning_rng(php, imghp),
# untagged) -- a direct re-examination of the exact case that showed the zmap collapse, not a
# fresh sample that might or might not reproduce it.
test_object, test_noise = End2EndThermalImg.with_tuning_rng(php, imghp) do
    (End2EndThermalImg.get_object(imghp), End2EndThermalImg.generate_noise(imghp))
end
noisy_image = End2EndThermalImg.make_noisy_image_from_3D_distributed(test_object, weights, test_noise, imghp.noise_type, php, imghp)

# ReconstructionHyperParams (src/prepare.jl) is a plain struct with a hand-written keyword
# constructor, not @with_kw -- see scripts/test_reconstruction_xtol.jl for the same helper.
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
rechp_test = with_xtol_rel(rechp, 1e-5; maxeval = 2000)
jhp_test = JobHyperParams(jhp; rechp = rechp_test)

# Overridable via E2E_ALPHAS/E2E_BETAS (comma-separated), e.g. for a finer/narrower follow-up
# sweep without duplicating this whole script -- same parsing style as endtoend_opt.jl's
# E2E_ALPHA_INIT/E2E_BETA_INIT. Defaults to the original 8x8 grid.
parse_float_list(s) = parse.(Float64, split(s, ","))
alphas = haskey(ENV, "E2E_ALPHAS") ? parse_float_list(ENV["E2E_ALPHAS"]) : [0.0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
betas = haskey(ENV, "E2E_BETAS") ? parse_float_list(ENV["E2E_BETAS"]) : [0.0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]

println("\n--- running (α, β) sweep: ", length(alphas), " x ", length(betas), " = ", length(alphas) * length(betas), " pairs ---")
summary_rows = NamedTuple[]
object_opts = Dict{Tuple{Float64,Float64}, Any}()

for α in alphas, β in betas
    local result, L, n_iterations, t
    try
        t = @elapsed result = End2EndThermalImg.reconstruct_Tmap_and_zmap_distributed(noisy_image, weights, α, β, jhp_test; objective_scale = :auto)
        L = End2EndThermalImg.L_reconstruction_error(test_object, result.object_opt)
        n_iterations = length(result.objective_history)
        object_opts[(α, β)] = result.object_opt
        push!(summary_rows, (; α, β, return_value = result.return_value, n_iterations, L, time_s = t))
        println("α=$α  β=$β  -> $(result.return_value) after $n_iterations iters, L=$(round(L, sigdigits=4)), $(round(t, digits=1))s")
    catch e
        @warn "Reconstruction failed" α β exception=(e, catch_backtrace())
        push!(summary_rows, (; α, β, return_value = :ERROR, n_iterations = 0, L = NaN, time_s = NaN))
    end
end

# --- Save results ---

root = "convergence_studies"
isdir(root) || mkdir(root)
datetime = now()
run_folder = joinpath(root, "alpha_beta_sweep_$(datetime)")
mkdir(run_folder)
results_file = joinpath(run_folder, "results.jld2")
@save results_file source_run_folder alphas betas test_object noisy_image summary_rows object_opts php imghp

println("\n=== Done ===")
println("run_folder   = ", run_folder)
println("results_file = ", results_file)
