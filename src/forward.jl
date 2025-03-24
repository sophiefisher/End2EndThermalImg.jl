"""
    convolve(inp, kernel)

Convolves an inpL x inpL array with the FFT of a centered kernel 
of size kerL x kerL to produce an output of size (kerL - inpL) x (kerL - inpL).

"""
# TODO: why does this assume the kernel is already fourier transformed?
function convolve(inp, kernel)
    # inpL < kerL
    inpL = size(inp, 1)
    kerL = size(kernel, 1)
    outL = kerL - inpL

    arr_pad = [inp zeros(inpL, outL); zeros(outL, inpL) zeros(outL, outL)]
    out_pad = ifft(fft(arr_pad) .* kernel)
    out = out_pad[inpL+1:kerL, inpL+1:kerL]
    out
end

# z is the distance between the object plane and the metasurface
function incident_field(freq, z, n, num_unit_cells, unit_cell_length)
    ω = 2 * π * freq
    k = n * ω

    function efield(x, y)
        r = √(x^2 + y^2 + z^2) 
        ℯ ^ (k * r * im) / ( 4 * π * r)
    end

    grid = range(-num_unit_cells / 2 + 0.5, num_unit_cells / 2 - 0.5, length = num_unit_cells) .* unit_cell_length
    incident = [efield(x, y) for x in grid, y in grid]
    incident
end

function get_incident_field(freq, z, php::PhysicsHyperParams)
    get_substrate_ϵ = get_permittivity_function(php.substrate_material)
    λ_µm = convert_freq_unitless_to_λ_µm(freq, php)
    substrate_ϵ = get_substrate_ϵ(λ_µm)
    incident = incident_field(freq, z, √(substrate_ϵ), php.num_unit_cells, php.unit_cell_length)
    incident
end

function get_near_field(incident_field, surrogate, geoms, sampleN)
    near = incident_field .* surrogate.(geoms)
    near = repeat(near, inner=(sampleN, sampleN))
end

function get_near_field(incident_field, surrogate, geoms, imghp::ImagingHyperParams)
    get_near_field(incident_field, surrogate, geoms, imghp.sampleN)
end

# TODO: implement absolute scaling factor for the green's functions
function n2f_kernel(freq, z, ϵ, μ, n2f_size, unit_cell_length, sampleN)
    ω = 2 * π * freq
    n = √(ϵ*μ)
    k = n * ω

    function efield(x, y)
        r = √(x^2 + y^2 + z^2)
        z * (-1 + k * r * im) * ℯ ^ (k * r * im) / (4 * π * r^3)
    end

    gridout = range(-(n2f_size ÷ 2), (n2f_size ÷ 2) - 1, length = n2f_size  ) .* (unit_cell_length / sampleN)
    fft([efield(x, y) * -μ / ϵ for x in gridout, y in gridout])
end

function get_n2f_kernel(freq, focal_length, num_unit_cells, unit_cell_length, psfN, binN, sampleN)
    n2f_size = (num_unit_cells + binN*psfN)*sampleN
    n2f_kernel(freq, focal_length, 1.0, 1.0, n2f_size, unit_cell_length, sampleN)
end

function get_n2f_kernel(freq, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    @unpack focal_length, num_unit_cells, unit_cell_length = php
    @unpack objN, imgN, binN, sampleN = imghp
    psfN = (objN + imgN)
    get_n2f_kernel(freq, focal_length, num_unit_cells, unit_cell_length, psfN, binN, sampleN)
end

function near_to_far_field(near_field, n2f_kernel)
    far = convolve(near_field, n2f_kernel)
end

function far_field_to_PSF(far_field, freq, unit_cell_length, binN, sampleN)
    far_field_abs = abs.(far_field).^2
    psfN = size(far_field, 1) ÷ sampleN ÷ binN
    far_field_abs_integrated = reshape(far_field_abs, (sampleN * binN, psfN, sampleN * binN, psfN))
    far_field_abs_integrated = sum(far_field_abs_integrated, dims=(1, 3)) 
    # (unit_cell_length / sampleN) is the integration/sampling width for integrating over each subpixel
    # divide by freq to turn energy into photon constructor
    # TODO: to normalize correctly, also need to divide by factor of hbar here
    PSF = dropdims(far_field_abs_integrated, dims=(1, 3)) .* (unit_cell_length / sampleN) ./ freq 
end

function far_field_to_PSF(far_field, freq, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    far_field_to_PSF(far_field, freq, php.unit_cell_length, imghp.binN, imghp.sampleN)
end

function get_PSF(freq, z, surrogate, geoms, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    incident = get_incident_field(freq, z, php)
    near = get_near_field(incident, surrogate, geoms, imghp)
    n2f_kernel = get_n2f_kernel(freq, php, imghp)
    far = near_to_far_field(near, n2f_kernel)
    PSF = far_field_to_PSF(far, freq, php, imghp)
end

# smoothness_order = 0 yields the triangle function
function get_discretized_δ_function(smoothness_order, Δz)
    if smoothness_order == Inf
        f(z) = z > 0 ? exp(-1 / z) : 0
    else
        f(z) = z > 0 ? z^(smoothness_order + 1) : 0
    end
    g(z) = f(z) / ( f(z) + f(1-z) )
    δ(z) = (-g(z ./ Δz) - g(-z ./ Δz) + 1) .* (1/Δz)
end

function get_discretized_δ_function(imghp::ImagingHyperParams)
    get_discretized_δ_function(imghp.smoothness_order, imghp.PSF_Δz)
end

function get_black_body_spectrum(Tmap, freqs, php::PhysicsHyperParams)
    b = [(2 .* freq ^3 ) ./ (exp.(ħ .* (freq .* c .* 10^6 / php.wavcen) ./ (kB .* Tmap) ) .- 1) for freq in freqs]
    #b = reduce(hcat, b)
    #b = reshape(b, size(Tmap, 1), size(Tmap, 1), :)
    b
end

function convolve_PSF_with_black_body_spectrum(PSF, b)
    fftPSF = fft(PSF)
    convolve(b, fftPSF)
end