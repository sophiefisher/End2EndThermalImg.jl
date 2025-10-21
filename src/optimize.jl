function nlopt_wrap_objective_autodiff(params_flat::Vector, grad::Vector, objective, objective_history; iteration_print = 50, verbose = false)
    obj_val = objective(params_flat)
    push!(objective_history, obj_val)
    iteration_count = length(objective_history)
    if verbose
        if iteration_count % iteration_print == 0
            @info "Iteration $iteration_count: objective = $obj_val"
        end
    end

    if length(grad) > 0
        grad[:] = Zygote.gradient(x -> objective(x), params_flat)[1]
    end
    obj_val
end