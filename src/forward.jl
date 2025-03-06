"""
    convolve(inp, kernel)

Convolves an inpL x inpL array with the FFT of a centered kernel 
of size kerL x kerL to produce an output of size (kerL - inpL) x (kerL - inpL).

"""
# TODO: why does this assume the kernel is already fourier transformed?
function convolve(inp, kernel)
    # inpL < kerL
    inpL = size(inp)[1]
    kerL = size(kernel)[1]
    outL = kerL - inpL

    arr_pad = [inp zeros(inpL, outL); zeros(outL, inpL) zeros(outL, outL)]
    out_pad = ifft(fft(arr_pad) .* kernel)
    out = out_pad[inpL+1:kerL, inpL+1:kerL]
    out
end

function incident_field(z, freq, n, num_unit_cells, unit_cell_length)
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

function get_near_field(incident_field, surrogate, geoms)
    near = incident_field .* surrogate.(geoms)
    near
end

# TODO: implement absolute scaling factor for the green's functions
# TODO: this might change when i implement image sampling
function n2f_kernel(z, freq, ϵ, μ, n2f_size, unit_cell_length)
    ω = 2 * π * freq
    n = √(ϵ*μ)
    k = n * ω

    function efield(x, y)
        r = √(x^2 + y^2 + z^2)
        z * (-1 + k * r * im) * ℯ ^ (k * r * im) / (4 * π * r^3)
    end

    gridout = range(-(n2f_size ÷ 2), (n2f_size ÷ 2) - 1, length = n2f_size  ) .* unit_cell_length
    fft([efield(x, y) * -μ / ϵ for x in gridout, y in gridout])
end

function get_n2f_kernel(num_unit_cells, unit_cell_length, psfN, binN, freq, z)
    # TODO: this might change when i implement image sampling
    n2f_size = num_unit_cells + binN*psfN
    n2f_kernel(z, freq, 1.0, 1.0, n2f_size, unit_cell_length)
end

function get_n2f_kernel(jhp::JobHyperParams, freq, z)
    php, imghp = jhp.php, jhp.imghp
    @unpack num_unit_cells, unit_cell_length = php
    @unpack objN, imgN, binN = imghp
    psfN = (objN + imgN)
    get_n2f_kernel(num_unit_cells, unit_cell_length, psfN, binN, freq, z)
end

function near_to_far_field(near_field, n2f_kernel)
    far = convolve(near_field, n2f_kernel)
end

function far_field_to_PSF(far_field, binN, freq)
    psfN = size(far_field)[1] ÷ binN
    far_field_reshaped = reshape(far_field, (binN, psfN, binN, psfN))
    far_field_reshaped_abs = (abs.(far_field_reshaped)).^2
    far_field_binned = sum(far_field_reshaped_abs, dims=(1, 3))
    PSF = dropdims(far_field_binned, dims=(1,3)) ./ freq 
end

function far_field_to_PSF(far_field, imghp::ImagingHyperParams, freq)
    far_field_to_PSF(far_field, imghp.binN, freq)
end