@memoize function make_plan(size::Tuple)
    plan_fft(zeros(ComplexF64, size), flags=FFTW.MEASURE)
end

@memoize function make_plan!(size::Tuple)
    plan_fft!(zeros(ComplexF64, size), flags=FFTW.MEASURE)
end   

function planned_fft(x)
    plan = ChainRulesCore.ignore_derivatives( ()-> make_plan(size(x)) )
    plan * x
end

function planned_ifft(x)
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
end

"""
    convolve(inp, kernel)

Convolves an inpL x inpL array with the FFT of a centered kernel 
of size kerL x kerL to produce an output of size (kerL - inpL) x (kerL - inpL).

"""
function convolve(inp, kernel)
    # inpL < kerL
    inpL = size(inp, 1)
    kerL = size(kernel, 1)
    outL = kerL - inpL

    arr_pad = [inp zeros(inpL, outL); zeros(outL, inpL) zeros(outL, outL)]
    out_pad = planned_ifft(planned_fft(arr_pad) .* kernel)
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
    near = repeat(near, inner=(sampleN, sampleN))
    near
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

function n2f_kernel(freq, z, ϵ, μ, n2f_size, unit_cell_length, sampleN)
    ω = 2 * π * freq
    n = √(ϵ*μ)
    k = n * ω
    gridout = range(-(n2f_size ÷ 2), (n2f_size ÷ 2) - 1, length = n2f_size  ) .* (unit_cell_length / sampleN)
    n2f_kernel = planned_fft!([efield_n2f_greens(x, y, z, k, ϵ, μ) for x in gridout, y in gridout])
    n2f_kernel
end

function get_n2f_size(php::PhysicsHyperParams, imghp::ImagingHyperParams)
    @unpack num_unit_cells = php
    @unpack objN, imgN, binN, sampleN = imghp
    psfN = (objN + imgN)
    n2f_size = (num_unit_cells + binN*psfN)*sampleN
    n2f_size
end

function get_n2f_kernel(freq, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    @unpack focal_length, unit_cell_length = php
    @unpack sampleN = imghp
    n2f_size = get_n2f_size(php, imghp)
    out = n2f_kernel(freq, focal_length, 1.0, 1.0, n2f_size, unit_cell_length, sampleN)
    out
end

function get_n2f_kernels(freqs, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    n2f_kernels = [get_n2f_kernel(freq, php, imghp) for freq in freqs]
    n2f_kernels
end

function near_to_far_field(near_field, n2f_kernel)
    far = convolve(near_field, n2f_kernel)
    far
end

# TODO: should I really be dividing by freq to get photon count here?
function far_field_to_PSF(far_field, freq, unit_cell_length, binN, sampleN)
    far_field_abs = abs.(far_field).^2
    psfN = size(far_field, 1) ÷ sampleN ÷ binN
    far_field_abs_integrated = reshape(far_field_abs, (sampleN * binN, psfN, sampleN * binN, psfN))
    far_field_abs_integrated = sum(far_field_abs_integrated, dims=(1, 3)) 
    # (unit_cell_length / sampleN) is the integration/sampling width for integrating over each subpixel
    # divide by freq to turn energy into photon count
    # TODO: to normalize correctly, also need to divide by factor of hbar here
    PSF = dropdims(far_field_abs_integrated, dims=(1, 3)) .* (unit_cell_length / sampleN) ./ freq 
    PSF
end

function far_field_to_PSF(far_field, freq, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    far_field_to_PSF(far_field, freq, php.unit_cell_length, imghp.binN, imghp.sampleN)
end

function get_PSF_at_freq_and_z(freq, incident, surrogate, geoms, n2f_kernel, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    near = get_near_field(incident, surrogate, geoms, imghp)
    far = near_to_far_field(near, n2f_kernel)
    PSF = far_field_to_PSF(far, freq, php, imghp)
    PSF
end

function get_PSFs_at_z(freqs, incidents_at_z, surrogates, geoms, n2f_kernels, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    [get_PSF_at_freq_and_z(freqs[iF], incidents_at_z[iF], surrogates[iF], geoms, n2f_kernels[iF], php, imghp) for iF in eachindex(freqs)]
end

function get_PSFs(freqs, incidents, surrogates, geoms, n2f_kernels, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    PSF_zlen = imghp.PSF_zlen
    #[get_PSF(freqs[iF], incidents[iF, iZ], surrogates[iF], geoms, n2f_kernels[iF], php, imghp) for iF in eachindex(freqs), iZ in 1:PSF_zlen]
    pmap(CartesianIndices((eachindex(freqs),1:PSF_zlen))) do i 
        iF = i[1]
        iZ = i[2]
        get_PSF_at_freq_and_z(freqs[iF], incidents[iF, iZ], surrogates[iF], geoms, n2f_kernels[iF], php, imghp)
    end
end

function f_δ(smoothness_order, z)
    if smoothness_order == Inf
        return z > 0.0 ? exp(-1 / z) : 0.0
    else
        return z > 0.0 ? z^(smoothness_order + 1) : 0.0
    end
end

# smoothness_order = 0 yields the triangle function
function get_discretized_δ_function(smoothness_order, Δz)
    g = z -> f_δ(smoothness_order, z) / ( f_δ(smoothness_order, z) + f_δ(smoothness_order, 1-z) )
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
    #b = reduce(hcat, b);
    #b = reshape(b, size(Tmap_zslice,1), size(Tmap_zslice,1), :)
    b
end

# fixed freq and fixed z
function convolve_with_PSF(PSF, b_freqslice)
    fftPSF = planned_fft(PSF)
    out = real.(convolve(b_freqslice, fftPSF))
    out
end

function make_image_at_z(spectrum, incidents_at_z, geoms, n2f_kernels, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    # b = get_black_body_spectrum(Tmap_zslice, php)
    freqs = ChainRulesCore.ignore_derivatives( ()-> get_freq_chebpoints(php))
    surrogates = ChainRulesCore.ignore_derivatives( ()-> load_surrogate_models(php))
    weights = ChainRulesCore.ignore_derivatives( ()-> get_clenshaw_curtis_quadrature_weights(php))
    PSFs = ChainRulesCore.ignore_derivatives( ()-> get_PSFs_at_z(freqs, incidents_at_z, surrogates, geoms, n2f_kernels, php, imghp))
    image = sum(weights .* map(convolve_with_PSF, PSFs, spectrum))
    image
end

function make_image_at_z(spectrum, freqs, incidents_at_z, surrogates, geoms, n2f_kernels, weights, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    PSFs = ChainRulesCore.ignore_derivatives( ()-> get_PSFs_at_z(freqs, incidents_at_z, surrogates, geoms, n2f_kernels, php, imghp))
    image = sum(weights .* map(convolve_with_PSF, PSFs, spectrum))
    image
end

function make_image_at_z(spectrum, PSFs, weights)
    image = sum(weights .* map(convolve_with_PSF, PSFs, spectrum))
    image
end

function make_image_from_3D_inplace(object, incidents, geoms, n2f_kernels, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    δ_Δz = get_discretized_δ_function(imghp)
    PSF_zcoords = get_PSF_zcoords(imghp)

    C_interp_3D = zeros(imghp.objN, imghp.objN, imghp.PSF_zlen) 
    Tmap_indices =  CartesianIndices((1:imghp.objN, 1:imghp.objN))
    for Tmap_idx in Tmap_indices
        z = object.zmap[Tmap_idx]
        PSF_zlower_idx = searchsortedlast(PSF_zcoords,z)
        PSF_zlower = PSF_zcoords[PSF_zlower_idx]
        PSF_zupper_idx = PSF_zlower_idx + 1
        PSF_zupper = PSF_zcoords[PSF_zupper_idx]

        C_zlower = δ_Δz(PSF_zlower - z)
        C_zupper = δ_Δz(PSF_zupper - z)
        C_interp_3D[Tmap_idx, PSF_zlower_idx] = C_zlower
        C_interp_3D[Tmap_idx, PSF_zupper_idx] = C_zupper
    end
    B = get_black_body_spectrum(object.Tmap, php)

    image = zeros(imghp.imgN, imghp.imgN)
    for iZ = eachindex(PSF_zcoords)
        spectrum = [b .* C_interp_3D[:, :, iZ] for b in B]
        image_at_z = make_image_at_z(spectrum, incidents[:, iZ], geoms, n2f_kernels, php, imghp)
        image = image + image_at_z 
    end
    # TODO: need to add noise # actually I think noise should be added in a separate function
    # TODO: rewrite in-place operations
    image
end

# TODO: pass the variables wrapped in ignore_derivatives to the function directly
function make_image_from_3D_oop(object, incidents, geoms, n2f_kernels, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    δ_Δz = ChainRulesCore.ignore_derivatives( ()-> get_discretized_δ_function(imghp.smoothness_order, imghp.PSF_Δz))
    PSF_zcoords = ChainRulesCore.ignore_derivatives( ()-> get_PSF_zcoords(imghp))
    freqs = ChainRulesCore.ignore_derivatives( ()-> get_freq_chebpoints(php))
    surrogates = ChainRulesCore.ignore_derivatives( ()-> load_surrogate_models(php))
    weights = ChainRulesCore.ignore_derivatives( ()-> get_clenshaw_curtis_quadrature_weights(php))

    C_interp_3D = [let
        z = object.zmap[Tmap_idx1, Tmap_idx2]
        PSF_zlower_idx = searchsortedlast(PSF_zcoords, z)
        PSF_zupper_idx = PSF_zlower_idx + 1

        PSF_zlower = PSF_zcoords[PSF_zlower_idx]
        PSF_zupper = PSF_zcoords[PSF_zupper_idx]

        C_zlower = δ_Δz(PSF_zlower - z)
        C_zupper = δ_Δz(PSF_zupper - z)
        i == PSF_zlower_idx ? C_zlower : i == PSF_zupper_idx ? C_zupper : 0.0
    end for  Tmap_idx1 in 1:imghp.objN, Tmap_idx2 in 1:imghp.objN, i in 1:imghp.PSF_zlen]
    B = get_black_body_spectrum(object.Tmap, php)

    image = sum(
        iZ -> begin
            spectrum = [b .* C_interp_3D[:, :, iZ] for b in B]
            image_at_z = make_image_at_z(spectrum, freqs, incidents[:, iZ], surrogates, geoms, n2f_kernels, weights, php, imghp)
        end,
    eachindex(PSF_zcoords))
    # TODO: need to add noise # actually I think noise should be added in a separate function
    image
end