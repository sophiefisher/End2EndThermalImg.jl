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

# TODO: implement absolute scaling factor for the green's functions
function greens(z, freq, ϵ, μ, n2f_size, unit_cell_length)
    ω = 2 * π * freq
    n = √(ϵ*μ)
    k = n * ω

    function efield(x, y)
        r = √(x^2 + y^2 + z^2)
        z * (-1 + k * r * im) * ℯ ^ (k * r * im) / (4 * π * r^3)
    end

    gridout = range(-(n2f_size ÷ 2), (n2f_size ÷ 2) - 1, length = n2f_size  ) .* unit_cell_length
    [efield(x, y) * -μ / ϵ for x in gridout, y in gridout]
end

# assumes incidence from the infinite substrate
function get_incident_field(php::PhysicsHyperParams, freq::AbstractFloat, z::AbstractFloat)
    get_substrate_ϵ = get_permittivity_function(php.substrate_material)
    λ_µm = convert_freq_unitless_to_λ_µm(freq, php)
    substrate_ϵ = get_substrate_ϵ(λ_µm)
    incident = incident_field(z, freq, √(substrate_ϵ), php.num_unit_cells, php.unit_cell_length)
    incident
end