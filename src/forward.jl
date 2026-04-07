#= @memoize function make_plan(size::Tuple)
    # TODO: do i want to thread the ffts?
    plan_fft(zeros(ComplexF64, size), flags=FFTW.MEASURE)
end

@memoize function make_plan!(size::Tuple)
    # TODO: do i want to thread the ffts?
    plan_fft!(zeros(ComplexF64, size), flags=FFTW.MEASURE)
end   

function planned_fft(x)
    # TODO: do i want to wrap this in ignore_derivatives?
    # if i write a rule for this function, i probably don't need it
    plan = ChainRulesCore.ignore_derivatives( ()-> make_plan(size(x)) )
    plan * x
end

function planned_ifft(x)
    # TODO: do i want to wrap this in ignore_derivatives?
    # if i write a rule for this function, i probably don't need it
    plan = ChainRulesCore.ignore_derivatives( ()-> make_plan(size(x)) )
    plan \ x
end

function planned_fft!(x)
    plan = make_plan!(size(x))
    plan * x
end

function planned_ifft!(x)
    plan = make_plan!(size(x))
    plan \ x
end =#

# TODO: think about whether these should be in-place
function get_fft_plans(php::PhysicsHyperParams, imghp::ImagingHyperParams)
    n2f_size = get_n2f_size(php, imghp)
    PSF_size = imghp.objN + imghp.imgN
    plans_n2f = Vector{FFTW.cFFTWPlan}(undef, nthreads())
    plans_PSF = Vector{FFTW.cFFTWPlan}(undef, nthreads())
    for t in 1:nthreads()
        plans_n2f[t] = plan_fft(zeros(ComplexF64, (n2f_size, n2f_size)), flags=FFTW.MEASURE)
        plans_PSF[t] = plan_fft(zeros(ComplexF64, (PSF_size, PSF_size)), flags=FFTW.MEASURE)
    end
    plans_n2f, plans_PSF
end

# Convolves an inpL x inpL array with the FFT of a centered kernel 
# of size kerL x kerL to produce an output of size (kerL - inpL) x (kerL - inpL).
function convolve(inp, kernel, plan)
    # inpL < kerL
    inpL = size(inp, 1)
    kerL = size(kernel, 1)
    outL = kerL - inpL

    arr_pad = [inp zeros(inpL, outL); zeros(outL, inpL) zeros(outL, outL)]
    out_pad = plan \ ((plan * arr_pad) .* kernel)

    out = out_pad[inpL+1:kerL, inpL+1:kerL]
    out
end

function efield_point_source(x, y, z, k)
    r = √(x^2 + y^2 + z^2) 
    ℯ ^ (k * r * im) / ( 4 * π * r)
end

# z is the distance between the object plane and the metasurface
function incident_field(freq, z, n, num_unit_cells, unit_cell_length)
    ω = 2 * π * freq
    k = n * ω
    grid = range(-num_unit_cells / 2 + 0.5, num_unit_cells / 2 - 0.5, length = num_unit_cells) .* unit_cell_length
    incident = [efield_point_source(x, y, z, k) for x in grid, y in grid]
    incident
end

function get_incident_field(freq, z, php::PhysicsHyperParams)
    get_substrate_ϵ = get_permittivity_function(php.substrate_material)
    λ_µm = convert_freq_unitless_to_λ_µm(freq, php)
    substrate_ϵ = get_substrate_ϵ(λ_µm)
    incident = incident_field(freq, z, √(substrate_ϵ), php.num_unit_cells, php.unit_cell_length)
    incident
end

function get_incident_fields(freqs, PSF_zcoords, php::PhysicsHyperParams)
    incidents = [get_incident_field(freq, z, php) for freq in freqs, z in PSF_zcoords]
    incidents
end

function get_near_field(incident_field, surrogate, geoms, sampleN)
    near = incident_field .* surrogate.(geoms)
    if sampleN == 1
        return near
    else
        return repeat(near, inner=(imghp.sampleN, imghp.sampleN))
    end
end

function get_near_field(incident_field, surrogate, geoms, imghp::ImagingHyperParams)
    get_near_field(incident_field, surrogate, geoms, imghp.sampleN)
end

# TODO: implement absolute scaling factor for the green's functions
# TODO: i think i can drop the -1 term
function efield_n2f_greens(x, y, z, k, ϵ, μ)
    r = √(x^2 + y^2 + z^2)
    z * (-1 + k * r * im) * ℯ ^ (k * r * im) / (4 * π * r^3) * (-μ / ϵ)
end

function n2f_kernel(freq, z, ϵ, μ, n2f_size, unit_cell_length, sampleN, plan_n2f)
    ω = 2 * π * freq
    n = √(ϵ*μ)
    k = n * ω
    gridout = range(-(n2f_size ÷ 2), (n2f_size ÷ 2) - 1, length = n2f_size  ) .* (unit_cell_length / sampleN)
    n2f_kernel = plan_n2f * [efield_n2f_greens(x, y, z, k, ϵ, μ) for x in gridout, y in gridout]
    n2f_kernel
end

function get_n2f_size(php::PhysicsHyperParams, imghp::ImagingHyperParams)
    @unpack num_unit_cells = php
    @unpack objN, imgN, binN, sampleN = imghp
    psfN = (objN + imgN)
    n2f_size = (num_unit_cells + binN*psfN)*sampleN
    n2f_size
end

function get_n2f_kernel(freq, plan_n2f, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    @unpack focal_length, unit_cell_length = php
    @unpack sampleN = imghp
    n2f_size = get_n2f_size(php, imghp)
    out = n2f_kernel(freq, focal_length, 1.0, 1.0, n2f_size, unit_cell_length, sampleN, plan_n2f)
    out
end

function get_n2f_kernels(freqs, plan_n2f, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    n2f_kernels = [get_n2f_kernel(freq, plan_n2f, php, imghp) for freq in freqs]
    n2f_kernels
end

function near_to_far_field(near_field, n2f_kernel, plan_n2f)
    far = convolve(near_field, n2f_kernel, plan_n2f)
    far
end

# TODO: should I really be dividing by freq to get photon count here?
function far_field_to_PSF(far_field, freq, unit_cell_length, binN, sampleN, PSF_scale)
    far_field_abs = abs.(far_field).^2
    psfN = size(far_field, 1) ÷ sampleN ÷ binN
    far_field_abs_integrated = reshape(far_field_abs, (sampleN * binN, psfN, sampleN * binN, psfN))
    far_field_abs_integrated = sum(far_field_abs_integrated, dims=(1, 3)) 
    # (unit_cell_length / sampleN) is the integration/sampling width for integrating over each subpixel
    # divide by freq to turn energy into photon count
    PSF = dropdims(far_field_abs_integrated, dims=(1, 3)) .* (unit_cell_length / sampleN) .* PSF_scale ./ freq 
    PSF
end

function far_field_to_PSF(far_field, freq, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    far_field_to_PSF(far_field, freq, php.unit_cell_length, imghp.binN, imghp.sampleN, imghp.PSF_scale)
end

function get_PSF_at_freq_and_z(freq, incident, surrogate, geoms, n2f_kernel, plan_n2f, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    near = get_near_field(incident, surrogate, geoms, imghp)
    far = near_to_far_field(near, n2f_kernel, plan_n2f)
    PSF = far_field_to_PSF(far, freq, php, imghp)
    PSF
end

# function get_PSFs_distributed(freqs, incidents, surrogates, geoms, n2f_kernels, php::PhysicsHyperParams, imghp::ImagingHyperParams)
#     PSF_zlen = imghp.PSF_zlen
#     #[get_PSF(freqs[iF], incidents[iF, iZ], surrogates[iF], geoms, n2f_kernels[iF], php, imghp) for iF in eachindex(freqs), iZ in 1:PSF_zlen]
#     pmap(CartesianIndices((eachindex(freqs),1:PSF_zlen))) do i 
#         iF = i[1]
#         iZ = i[2]
#         get_PSF_at_freq_and_z(freqs[iF], incidents[iF, iZ], surrogates[iF], geoms, n2f_kernels[iF], php, imghp)
#     end
# end

function get_PSFs(freqs, incidents, surrogates, geoms, n2f_kernels, plans_n2f,
                           php::PhysicsHyperParams, imghp::ImagingHyperParams)

    PSF_zlen = imghp.PSF_zlen
    nF = length(freqs)
    # TODO add size of PSF below? 
    PSFs = Matrix{Matrix{Float64}}(undef, nF, PSF_zlen)
    inds = CartesianIndices((1:nF, 1:PSF_zlen))

    @threads for idx in eachindex(inds)
        tid = threadid()
        I = inds[idx]
        iF, iZ = I[1], I[2]

        PSFs[iF, iZ] = get_PSF_at_freq_and_z(
            freqs[iF],
            incidents[iF, iZ],
            surrogates[iF],
            geoms,
            n2f_kernels[iF],
            plans_n2f[tid],
            php,
            imghp
        )
    end
    return PSFs
end

# TODO: make in-place?
get_fftPSF(PSF, plan_PSF) = plan_PSF * complex.(PSF)

# function get_fftPSFs_distributed(freqs, incidents, surrogates, geoms, n2f_kernels, php::PhysicsHyperParams, imghp::ImagingHyperParams)
#     PSF_zlen = imghp.PSF_zlen
#     pmap(CartesianIndices((eachindex(freqs),1:PSF_zlen))) do i 
#         iF = i[1]
#         iZ = i[2]
#         get_fftPSF(get_PSF_at_freq_and_z(freqs[iF], incidents[iF, iZ], surrogates[iF], geoms, n2f_kernels[iF], php, imghp))
#     end
# end

function get_fftPSFs(freqs, incidents, surrogates, geoms, n2f_kernels, plans_n2f, plans_PSF,
                           php::PhysicsHyperParams, imghp::ImagingHyperParams)

    PSF_zlen = imghp.PSF_zlen
    nF = length(freqs)
    fftPSFs = Matrix{Matrix{ComplexF64}}(undef, nF, PSF_zlen)
    inds = CartesianIndices((1:nF, 1:PSF_zlen))

    @threads for idx in eachindex(inds)
        tid = threadid()
        I = inds[idx]
        iF, iZ = I[1], I[2]

        fftPSFs[iF, iZ] = get_fftPSF(get_PSF_at_freq_and_z(
            freqs[iF],
            incidents[iF, iZ],
            surrogates[iF],
            geoms,
            n2f_kernels[iF],
            plans_n2f[tid],
            php,
            imghp
        ), plans_PSF[tid])
    end
    return fftPSFs
end

function f_δ(z, smoothness_order)
    if smoothness_order == Inf
        return z > 0.0 ? exp(-1 / z) : 0.0
    else
        return z > 0.0 ? z^(smoothness_order + 1) : 0.0
    end
end

# smoothness_order = 0 yields the triangle function
function get_discretized_δ_function(smoothness_order, Δz)
    g = z -> f_δ(z, smoothness_order) / ( f_δ(z, smoothness_order) + f_δ(1-z, smoothness_order) )
    δ = z -> (-g(z ./ Δz) - g(-z ./ Δz) + 1) .* (1/Δz)
    δ
end

function get_discretized_δ_function(imghp::ImagingHyperParams)
    δ = get_discretized_δ_function(imghp.smoothness_order, imghp.PSF_Δz)
    δ
end

function get_black_body_spectrum(Tmap_zslice, php::PhysicsHyperParams)
    freqs = ChainRulesCore.ignore_derivatives( ()-> get_freq_chebpoints(php))
    b = [(2 .* freq ^3 ) ./ (exp.(ħ .* (freq .* c .* 10^6 / php.wavcen) ./ (kB .* Tmap_zslice) ) .- 1) for freq in freqs]
    b
end

convolve_with_fftPSF_at_freq_and_z(fftPSF_at_freq_and_z, b_at_freq_and_z, plan_PSF) = real.(convolve(b_at_freq_and_z, fftPSF_at_freq_and_z, plan_PSF))

function make_image_at_z(spectrum_at_z, fftPSFs_at_z, weights, plan_PSF)
    image = sum(
        weights .* convolve_with_fftPSF_at_freq_and_z.(fftPSFs_at_z,
                                                       spectrum_at_z,
                                                       Ref(plan_PSF))
    )
    return image
end

# function make_image_from_3D!(image_buf, object, fftPSFs, weights, php::PhysicsHyperParams, imghp::ImagingHyperParams)
#     δ_Δz = get_discretized_δ_function(imghp)
#     PSF_zcoords = get_PSF_zcoords(imghp)

#     # TODO: turn repeated code into function (see make_image_from_3D)
#     C_interp_3D = zeros(imghp.objN, imghp.objN, imghp.PSF_zlen) # TODO: allocations
#     Tmap_indices = CartesianIndices((1:imghp.objN, 1:imghp.objN))
#     for Tmap_idx in Tmap_indices
#         z = object.zmap[Tmap_idx]
#         if z == PSF_zcoords[end]
#             PSF_zlower_idx = imghp.PSF_zlen - 1
#             PSF_zupper_idx = imghp.PSF_zlen
#         else
#             PSF_zlower_idx = searchsortedlast(PSF_zcoords, z)
#             PSF_zupper_idx = PSF_zlower_idx + 1
#         end
#         PSF_zlower = PSF_zcoords[PSF_zlower_idx]
#         PSF_zupper = PSF_zcoords[PSF_zupper_idx]

#         C_zlower = δ_Δz(PSF_zlower - z)
#         C_zupper = δ_Δz(PSF_zupper - z)
#         C_interp_3D[Tmap_idx, PSF_zlower_idx] = C_zlower
#         C_interp_3D[Tmap_idx, PSF_zupper_idx] = C_zupper
#     end
#     B = get_black_body_spectrum(object.Tmap, php)

#     fill!(image_buf, 0.0)
#     @views for iZ in eachindex(PSF_zcoords)
#         spectrum = [b .* C_interp_3D[:, :, iZ] for b in B]  # TODO: allocations
#         image_at_z = make_image_at_z(spectrum, fftPSFs[:, iZ], weights) # TODO: allocations
#         @. image_buf += image_at_z
#     end
#     return image_buf
# end

function get_C_interp_3D(zmap, PSF_zcoords, δ_Δz, imghp::ImagingHyperParams)
    return [let
        z = zmap[Tmap_idx1, Tmap_idx2]
        if z == PSF_zcoords[end]
            PSF_zlower_idx = imghp.PSF_zlen - 1
            PSF_zupper_idx = imghp.PSF_zlen
        else
            PSF_zlower_idx = searchsortedlast(PSF_zcoords, z)
            PSF_zupper_idx = PSF_zlower_idx + 1
        end

        PSF_zlower = PSF_zcoords[PSF_zlower_idx]
        PSF_zupper = PSF_zcoords[PSF_zupper_idx]

        C_zlower = δ_Δz(PSF_zlower - z)
        C_zupper = δ_Δz(PSF_zupper - z)
        i == PSF_zlower_idx ? C_zlower : i == PSF_zupper_idx ? C_zupper : 0.0
    end for  Tmap_idx1 in 1:imghp.objN, Tmap_idx2 in 1:imghp.objN, i in 1:imghp.PSF_zlen]
end

# TODO: pass the variables wrapped in ignore_derivatives to the function directly
function make_image_from_3D(object, fftPSFs, weights, plan_PSF, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    δ_Δz = ignore_derivatives( ()-> get_discretized_δ_function(imghp.smoothness_order, imghp.PSF_Δz))
    PSF_zcoords = ignore_derivatives( ()-> get_PSF_zcoords(imghp))

    C_interp_3D = get_C_interp_3D(object.zmap, PSF_zcoords, δ_Δz, imghp)
    B = get_black_body_spectrum(object.Tmap, php)

    image = sum(
        iZ -> begin
            spectrum = [b .* C_interp_3D[:, :, iZ] for b in B]
            image_at_z = make_image_at_z(spectrum, fftPSFs[:, iZ], weights, plan_PSF)
        end,
    eachindex(PSF_zcoords))
    image
end

# generate_noise!(noise_buf) = randn!(noise_buf)

# # generates noise inside function
# function make_noisy_image_from_3D!(
#     image_buf, noise_buf,
#     object, fftPSFs, weights,
#     php::PhysicsHyperParams, imghp::ImagingHyperParams)

#     make_image_from_3D!(image_buf, object, fftPSFs, weights, php, imghp)
#     noise_scale = mean(image_buf) * imghp.noise_level
#     generate_noise!(noise_buf)
#     @. image_buf += noise_scale * noise_buf
#     return image_buf
# end

generate_noise(imghp) = randn((imghp.imgN, imghp.imgN))

# pass noise to function
function make_noisy_image_from_3D(
    object, fftPSFs, weights, noise, noise_scale, plan_PSF, php::PhysicsHyperParams, imghp::ImagingHyperParams)

    image = make_image_from_3D(object, fftPSFs, weights, plan_PSF, php, imghp)
    noisy_image = image .+ noise_scale .* noise
    return noisy_image
end

function reconstruction_objective(object_flat, noisy_image, fftPSFs, weights, plan_PSF, α, β, T_background, z_middle, php, imghp)
    object = unflatten_object(object_flat)
    τmap = object.Tmap
    ζmap = object.zmap
    image = make_image_from_3D(object, fftPSFs, weights, plan_PSF, php, imghp)
    error_image = sum((image .- noisy_image).^2)
    regularization_τ = α * sum((τmap .- T_background).^2)
    regularization_ζ = β * sum((ζmap .- z_middle).^2)
    error_image + regularization_τ + regularization_ζ
end

# TODO: might want to set xtol_rel, maxeval, as rechp parameters
function reconstruct_Tmap_and_zmap(noisy_image, fftPSFs, weights, plan_PSF, α, β, jhp::JobHyperParams; xtol_rel = 1e-8, maxeval = 5000, iteration_print = 50, verbose = false)
    verbose && @info "Starting object reconstruction"
    @unpack php, imghp, rechp = jhp
    @unpack T_background, z_middle = rechp
    
    objective_history = Float64[]
    object_init = initialize_object(imghp, rechp)
    object_init_flat = flatten_object(object_init)
    
    opt = Opt(:LD_LBFGS, 2 * imghp.objN^2) # TODO: check algorithm choice
    lower_bounds!(opt, [fill(eps(), imghp.objN^2); fill(imghp.PSF_zlb, imghp.objN^2)])
    upper_bounds!(opt, [fill(Inf, imghp.objN^2); fill(imghp.PSF_zub, imghp.objN^2)])
    objective_lambda = object_flat -> reconstruction_objective(object_flat, noisy_image, fftPSFs, weights, plan_PSF, α, β, T_background, z_middle, php, imghp)
    objective_wrapped_lambda = (x, grad) -> nlopt_wrap_objective_autodiff(x, grad, objective_lambda, objective_history; iteration_print = iteration_print, verbose = verbose)
    min_objective!(opt, objective_wrapped_lambda)
    xtol_rel!(opt, xtol_rel)
    maxeval!(opt, maxeval)

    (objective_opt, object_opt_flat, return_value) = NLopt.optimize!(opt, object_init_flat)
    verbose && @info "Done object reconstruction"
    verbose && @info "Optimization results" objective_opt return_value
    object_opt = unflatten_object(object_opt_flat)
    (; objective_opt, object_opt, return_value, objective_history)
end

# function Tmap_reconstruction_objective(Tmap_flat, noisy_image, fftPSFs, weights, α, T_background, php, imghp)
#     τmap = unflatten_square_matrix(Tmap_flat)
#     ζmap = fill(imghp.object_type.z, imghp.objN, imghp.objN)
#     object = (Tmap = τmap, zmap = ζmap)
#     image = make_image_from_3D(object, fftPSFs, weights, php, imghp)
#     error_image = sum((image .- noisy_image).^2)
#     regularization_τ = α * sum((τmap .- T_background).^2)
#     error_image + regularization_τ 
# end

# # sets β = 0
# # only works for objects with a fixed depth
# function reconstruct_Tmap_for_fixed_depth(noisy_image, fftPSFs, weights, α, jhp::JobHyperParams; xtol_rel = 1e-8, maxeval = 5000, iteration_print = 50, verbose = false)
#     verbose && @info "Starting object reconstruction"
#     @unpack php, imghp, rechp = jhp
#     @unpack T_background = rechp
    
#     objective_history = Float64[]
#     object_init = initialize_object(imghp, rechp)
#     Tmap_init_flat = object_init.Tmap[:]
#     opt = Opt(:LD_LBFGS, imghp.objN^2) # TODO: check algorithm choice
#     lower_bounds!(opt, fill(eps(), imghp.objN^2))
#     upper_bounds!(opt, fill(Inf, imghp.objN^2))
#     objective_lambda = Tmap_flat -> Tmap_reconstruction_objective(Tmap_flat, noisy_image, fftPSFs, weights, α, T_background, php, imghp)
#     objective_wrapped_lambda = (x, grad) -> nlopt_wrap_objective_autodiff(x, grad, objective_lambda, objective_history; iteration_print = iteration_print, verbose = verbose)
#     min_objective!(opt, objective_wrapped_lambda)
#     xtol_rel!(opt, xtol_rel)
#     maxeval!(opt, maxeval)

#     (objective_opt, Tmap_opt_flat, return_value) = NLopt.optimize!(opt, Tmap_init_flat)
#     verbose && @info "Done Tmap reconstruction"
#     verbose && @info "Optimization results" objective_opt return_value
#     Tmap_opt = unflatten_square_matrix(Tmap_opt_flat)
#     zmap = fill(imghp.object_type.z, imghp.objN, imghp.objN)
#     object_opt = (Tmap = Tmap_opt, zmap = zmap)

#     (; objective_opt, object_opt, return_value, objective_history)
# end
