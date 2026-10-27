function nlopt_wrap_objective_autodiff(params_flat::Vector, grad::Vector, objective, objective_history; iteration_print = 50, verbose = false)
    # withgradient, not a separate objective(params_flat) call followed by Zygote.gradient: the
    # latter runs the forward pass twice per iteration (once directly, once again inside
    # Zygote.gradient, which rebuilds the forward pass as part of differentiating it).
    # withgradient gets both from a single combined forward+backward pass; confirmed its (val,
    # grad) are bit-identical to the old two-call approach.
    if length(grad) > 0
        result = Zygote.withgradient(objective, params_flat)
        obj_val = result.val
        grad[:] = result.grad[1]
    else
        obj_val = objective(params_flat)
    end
    push!(objective_history, obj_val)
    iteration_count = length(objective_history)
    if verbose
        if iteration_count % iteration_print == 0
            @info "Iteration $iteration_count: objective = $obj_val"
        end
    end
    obj_val
end

L_reconstruction_error(object, object_opt) =
    sum(abs2, object.Tmap .- object_opt.Tmap) / sum(abs2, object.Tmap) +
    sum(abs2, object.zmap .- object_opt.zmap) / sum(abs2, object.zmap)

# Stochastic outer optimization over [geoms; α; β], minimizing the expected reconstruction error
# L_reconstruction_error(object, object_opt) across a fresh random (object, noise) draw every
# iteration (get_object/generate_noise) -- each iteration cold-starts the inner reconstruction
# (reconstruct_Tmap_and_zmap_distributed) from scratch, since the previous iteration's object_opt
# is not a meaningful warm start once the true object itself has changed. geoms is updated via
# Optimisers.jl's Adam and projected back into [pillar_width_lb, pillar_width_ub] (Adam itself is
# unconstrained). opthp.maxeval is the number of outer iterations (the only stopping criterion for
# now).
#
# α/β are optimized in log-space (logα = log(α), logβ = log(β), with dL/d(logα) = dL/dα · α by the
# chain rule, similarly for β) rather than directly. This was adopted after a direct run hit NaN:
# α, β start many orders of magnitude smaller than geoms (e.g. 1e-7 vs. an O(1) geoms scale), but
# Adam's first step for every parameter is ≈ ±learning_rate regardless of that parameter's own
# scale (its bias-corrected moment estimates reduce to m̂₁ = g, v̂₁ = g² exactly on iteration 1).
# With a shared learning_rate tuned for geoms, that first step overshot α/β straight past 0, got
# clamped there, and α = β = 0 exactly removes the inner reconstruction's only pull keeping T near
# T_background / z near z_middle -- T was then free to approach its lower bound (≈0 Kelvin), where
# dB_dT's exp(A·freq/T) overflows in both numerator and denominator, giving Inf/Inf = NaN, which
# then propagated through object_opt -> the Hessian -> minres! -> outer_gradient_distributed ->
# Adam's update -> geoms/α/β all at once. Log-space avoids this at the source: exp(anything) is
# always > 0 (no clamping needed for α/β at all), and a fixed learning_rate becomes an
# approximately *relative* (multiplicative) step in the original α/β space instead of an absolute
# one, so the first step is proportionate to α/β's own current scale rather than overshooting it.
#
# geoms/logα/logβ are carried as a NamedTuple of 1-element vectors for logα/logβ (not bare
# Float64s): Optimisers.jl's Functors.jl-based tree traversal silently skips bare Number leaves
# inside a NamedTuple (confirmed empirically -- a gradient supplied for a bare Float64 field is
# just dropped, not applied, with no error), so wrapping a scalar parameter in a length-1 Vector is
# the standard workaround.
# run_folder: if given, checkpoints the current state (same payload save_optimization_state writes
# for the final save) every opthp.save_every iterations, plus unconditionally on the last iteration
# regardless of whether it lands on a save_every multiple -- so a cancelled/killed long run never
# loses more than save_every iterations of progress, and the very latest state is always captured.
# Left as nothing (no I/O at all) by default, e.g. for quick in-memory use from a notebook.
function train_outer_optimization_distributed(geoms, α, β, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp::JobHyperParams, PSF_scale;
    run_folder = nothing, datetime = now(), verbose = true)
    @unpack php, imghp, opthp, rechp = jhp

    model = (; geoms = copy(geoms), logα = [log(α)], logβ = [log(β)])
    tree = Optimisers.setup(Optimisers.Adam(opthp.learning_rate), model)

    history = NamedTuple[]

    for iteration in 1:opthp.maxeval
        local object, object_opt, result, minres_hist, L, α_cur, β_cur, n_active_T, n_active_z
        iteration_time = @elapsed begin
            object = get_object(imghp)
            noise = generate_noise(imghp)
            α_cur = exp(model.logα[1])
            β_cur = exp(model.logβ[1])

            update_geoms_and_fftPSFs!(model.geoms, freqs, php, imghp, PSF_scale)
            noisy_image = make_noisy_image_from_3D_distributed(object, weights, noise, imghp.noise_type, php, imghp)

            result = reconstruct_Tmap_and_zmap_distributed(noisy_image, weights, α_cur, β_cur, jhp; objective_scale = :auto)
            object_opt = result.object_opt

            idx = get_free_and_active_indices(object_opt, imghp)
            free_T, free_z = idx.free_T, idx.free_z
            n_active_T, n_active_z = length(idx.active_T), length(idx.active_z)
            N_adjoint = length(free_T) + length(free_z)

            minres_hist = nothing
            if N_adjoint == 0
                # Every T and z pixel pinned at its bound -- nothing is left as a free
                # reconstruction variable, so there's no adjoint equation to solve and the
                # geoms/α/β sensitivity is undefined. Skip this iteration's update entirely
                # (geoms/α/β left unchanged) rather than solving a degenerate 0x0 system.
                @warn "Every T and z pixel pinned at its bound on iteration $iteration -- skipping this iteration's update (geoms/α/β left unchanged)"
            else
                fftPSFs = let
                    fftPSFs = Matrix{Matrix{ComplexF64}}(undef, length(freqs), length(PSF_zcoords))
                    for p in workers()
                        grid_inds, shard = remotecall_fetch(p) do
                            (_shard_grid_inds, _fftPSFs_shard)
                        end
                        for idx in grid_inds
                            iF, iZ = Tuple(idx)
                            fftPSFs[iF, iZ] = shard[(iF, iZ)]
                        end
                    end
                    fftPSFs
                end

                b_T = vec(2 .* (object.Tmap .- object_opt.Tmap) ./ sum(abs2, object.Tmap))[free_T]
                b_z = vec(2 .* (object.zmap .- object_opt.zmap) ./ sum(abs2, object.zmap))[free_z]
                b = [b_T; b_z]

                A = build_hessian_operator_distributed(object, object_opt, fftPSFs, freqs, weights, driver_plans_PSF,
                                                         noise, imghp.noise_type, α_cur, β_cur, php, imghp)
                Λ, minres_hist = IterativeSolvers.minres!(zeros(N_adjoint), A, b; log = true)

                grad_full = outer_gradient_distributed(Λ, object, object_opt, model.geoms, fftPSFs, freqs, weights, driver_plans_PSF,
                                                         noise, imghp.noise_type, rechp, php, imghp, PSF_scale)

                dL_dα, dL_dβ = grad_full[end - 1], grad_full[end]
                grads = (; geoms = reshape(grad_full[1:end - 2], size(model.geoms)), logα = [dL_dα * α_cur], logβ = [dL_dβ * β_cur])

                tree, model = Optimisers.update!(tree, model, grads)
                clamp!(model.geoms, php.pillar_width_lb, php.pillar_width_ub)
            end

            L = L_reconstruction_error(object, object_opt)
        end

        entry = (;
            iteration,
            L,
            α = α_cur,
            β = β_cur,
            iteration_time,
            reconstruction_return_value = result.return_value,
            reconstruction_iterations = length(result.objective_history),
            minres_converged = minres_hist === nothing ? missing : minres_hist.isconverged,
            minres_iterations = minres_hist === nothing ? missing : minres_hist.iters,
            n_active_T,
            n_active_z
        )
        push!(history, entry)
        verbose && @info "Iteration $iteration" entry.L entry.α entry.β entry.iteration_time entry.reconstruction_return_value entry.reconstruction_iterations entry.minres_converged entry.minres_iterations entry.n_active_T entry.n_active_z

        if run_folder !== nothing && (iteration % opthp.save_every == 0 || iteration == opthp.maxeval)
            save_optimization_state(run_folder, php, imghp, opthp, rechp, PSF_scale,
                                     model.geoms, exp(model.logα[1]), exp(model.logβ[1]), history, datetime)
            verbose && @info "Checkpointed at iteration $iteration" run_folder
        end
    end

    (; geoms = model.geoms, α = exp(model.logα[1]), β = exp(model.logβ[1]), history)
end

function ensure_optdata_folder()
    folder = "optdata"
    if !isdir(folder)
        mkdir(folder)
        @info "Created folder: $folder"
    else
        @info "Folder already exists: $folder"
    end
    return folder
end

# Writes the full current state to run_folder/results.jld2 (php/imghp/opthp/rechp/PSF_scale,
# geoms_opt_flat/α_opt/β_opt, history, datetime), run_folder/history.csv (the same history, for
# inspection outside Julia), and three plots (geoms.png, loss_history.png, alpha_beta_history.png).
# Called both periodically (checkpointing, from inside train_outer_optimization_distributed's loop)
# and for the final save -- same payload either way, just overwritten in place each time (not
# versioned per call), since the goal is "the latest state survives a cancelled job," not a full
# checkpoint history. geoms is saved flattened under the key geoms_opt_flat specifically so a saved
# run here is interchangeable with design_monochromatic_lens runs (inversedesign/*/results.jld2)
# via the existing geoms_init_type == "loaded" / load_geoms mechanism (src/inversedesign.jl) --
# either kind of saved run can seed either kind of future run.
function save_optimization_state(run_folder, php, imghp, opthp, rechp, PSF_scale, geoms, α, β, history, datetime)
    geoms_opt_flat = vec(geoms)
    α_opt = α
    β_opt = β

    results_file = joinpath(run_folder, "results.jld2")
    @save results_file php imghp opthp rechp PSF_scale geoms_opt_flat α_opt β_opt history datetime
    CSV.write(joinpath(run_folder, "history.csv"), DataFrame(history))

    plot_geoms(geoms_opt_flat, php; savefile = "$(run_folder)/geoms.png")
    plot_objective_function_history([e.L for e in history]; savefile = "$(run_folder)/loss_history.png")
    plot_alpha_beta_history([e.α for e in history], [e.β for e in history]; savefile = "$(run_folder)/alpha_beta_history.png")

    results_file
end

# Creates optdata/opt_<datetime>/ and runs train_outer_optimization_distributed with checkpointing
# into it (every opthp.save_every iterations, plus unconditionally on the last one -- see
# train_outer_optimization_distributed). Since the loop itself guarantees the final iteration is
# saved, nothing further needs saving here once it returns.
#
# α/β here are whatever the caller passed as rechp.α_init/β_init -- if either is :auto, it's
# resolved via resolve_initial_α_β (src/tune_regularization.jl) before the main loop starts, with
# its own artifacts saved into run_folder/tune_regularization/ (see save_tune_regularization_state)
# so everything from this job lands in one place. If neither is :auto, this is a no-op and α/β
# are used exactly as given, unchanged from before this resolution step existed.
function run_and_save_outer_optimization_distributed(geoms, α, β, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp::JobHyperParams, PSF_scale;
    verbose = true)
    @unpack rechp = jhp
    root = ensure_optdata_folder()
    datetime = now()
    run_name = "opt_$(datetime)"
    run_folder = joinpath(root, run_name)
    mkdir(run_folder)

    if α === :auto || β === :auto
        α_start, β_start = 1e-7, 1e-7
        verbose && @info "α_init/β_init includes :auto -- resolving initial α/β before the main optimization"
        tune_out = resolve_initial_α_β(rechp, geoms, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp, PSF_scale;
            α_start = α_start, β_start = β_start, verbose = verbose)
        save_tune_regularization_state(run_folder, tune_out, α_start, β_start, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp, PSF_scale)
        α, β = tune_out.α, tune_out.β
        verbose && @info "Resolved initial α/β" α β
    end

    out = train_outer_optimization_distributed(geoms, α, β, freqs, weights, driver_plans_PSF, PSF_zcoords, jhp, PSF_scale;
        run_folder = run_folder, datetime = datetime, verbose = verbose)

    results_file = joinpath(run_folder, "results.jld2")
    verbose && @info "Saved results" file=results_file

    (; results_file, run_folder, geoms = out.geoms, α = out.α, β = out.β, history = out.history)
end
