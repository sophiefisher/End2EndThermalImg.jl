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

function get_n2f_kernel(freq, jhp::JobHyperParams)
    php, imghp = jhp.php, jhp.imghp
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