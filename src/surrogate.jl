@memoize function get_permittivity_function(material)
    if material == "Si_no_absorption"
        filepath = "materials/Si_n_Shkondin.csv"
        data = CSV.read(filepath, DataFrame)
        λ_µm = data[:,:wl]
        n = data[:,:n]
        itp = linear_interpolation(λ_µm, n)
        return get_permittivity(λ_µm) = itp(λ_µm) # wavelength must be in microns
    end
end

convert_freq_unitless_to_λ_µm(freq, wavcen) = wavcen / freq

function convert_freq_unitless_to_λ_µm(freq, php::PhysicsHyperParams)
    convert_freq_unitless_to_λ_µm(freq, php.wavcen)
end

function get_freq_chebpoints(php::PhysicsHyperParams)
    # first element in the list is fub, last is flb
    chebpoints(php.freq_order, php.freqlb, php.frequb) 
end

function get_width_chebpoints(php::PhysicsHyperParams)
    # first element in the list is ub, last is lb
    chebpoints(php.pillar_width_order, php.pillar_width_lb, php.pillar_width_ub)
end

function get_transmission(freq, pillar_width, pillar_height, pillar_epsilon, unit_cell_length, substrate_epsilon, nG)
    # nG is number of fourier components in RCWA (truncation order)
    L1 = pylist([unit_cell_length,0])
    L2 = pylist([0,unit_cell_length])
    N_grid_points = pyint(10000) # number of grid points for the patterned layers
    obj = grcwa[].obj(pyint(nG), L1, L2, freq, 0.0, 0.0, verbose = 1 )
    obj.Add_LayerUniform(3.0, substrate_epsilon)
    obj.Add_LayerGrid(pillar_height, N_grid_points, N_grid_points)
    obj.Add_LayerUniform(3.0, 1.0)  
    obj.Init_Setup()
    epgrid = numpy[].zeros((N_grid_points, N_grid_points))
    boundary = unit_cell_length / 2.0
    xarr = numpy[].linspace(-boundary, boundary, N_grid_points)
    yarr = numpy[].linspace(-boundary, boundary, N_grid_points)
    x, y = numpy[].meshgrid(xarr,yarr,indexing=pystr("ij") )
    ind = numpy[].logical_and(pyabs(x) < pillar_width / 2.0, pyabs(y) < pillar_width / 2.0)
    epgrid[ind] = pyfloat(1.0)
    epgrid = (pillar_epsilon - 1.0) * epgrid + 1.0
    obj.GridLayer_geteps(epgrid.flatten())
    obj.a0 = numpy[].zeros(2*obj.nG, dtype = pytype(pycomplex(1.0 )))
    obj.a0[0] = pycomplex(1.0)
    obj.bN = numpy[].zeros(2*obj.nG, dtype = pytype(pycomplex(1.0 )))
    aN, b0 = obj.GetAmplitudes(2, 0.1)
    transmission = aN[0]/obj.a0[0] * numpy[].sqrt(numpy[].sqrt(substrate_epsilon))
    pyconvert(Complex{Float64}, transmission)
end

function compute_surrogate_transmission_matrix(php::PhysicsHyperParams)
    @unpack pillar_height, unit_cell_length = php
    get_pillar_ϵ = get_permittivity_function(php.pillar_material)
    get_substrate_ϵ = get_permittivity_function(php.substrate_material)
    freq_chebpoints = get_freq_chebpoints(php)
    width_chebpoints = get_width_chebpoints(php)

    transmission_matrix = Matrix{Complex{Float64}}(undef, length(freq_chebpoints), length(width_chebpoints))
    PythonCall.GIL.@unlock Threads.@threads for j in eachindex(width_chebpoints)
    width = width_chebpoints[j]
        PythonCall.GIL.@lock for (i, freq) in enumerate(freq_chebpoints)
            λ_µm = convert_freq_unitless_to_λ_µm(freq, php)
            pillar_ϵ = get_pillar_ϵ(λ_µm)
            substrate_ϵ = get_substrate_ϵ(λ_µm)
            transmission = get_transmission(freq, width, pillar_height, pillar_ϵ, unit_cell_length, substrate_ϵ, php.nG)
            transmission_matrix[i, j] = transmission
            GC.gc()
        end
    end
    transmission_matrix
end

# surrogate is labeled with unitful quantities, but computes models for unitless quantities
function get_surrogate_label(php::PhysicsHyperParams)
    @unpack λlb_μm, λub_μm, freq_order = php
    @unpack unit_cell_length_μm, pillar_width_lb_μm, pillar_width_ub_μm, pillar_width_order, pillar_height_μm  = php
    @unpack pillar_material, substrate_material, nG = php
    "surrogate_$(λlb_μm)_$(λub_μm)_$(freq_order)_$(unit_cell_length_μm)_$(pillar_width_lb_μm)_$(pillar_width_ub_μm)_$(pillar_width_order)_$(pillar_height_μm)_$(pillar_material)_$(substrate_material)_$(nG)"
end

get_surrogate_filename(php::PhysicsHyperParams) = "surdata/$(get_surrogate_label(php)).csv"

function compute_and_save_surrogate_transmission_matrix(php::PhysicsHyperParams)
    transmission_matrix = compute_surrogate_transmission_matrix(php)
    filename = get_surrogate_filename(php)
    CSV.write(filename, Tables.table(transmission_matrix) )
end

# one for each frequency
function load_surrogate_models(php::PhysicsHyperParams)
    filepath = get_surrogate_filename(php)
    if isfile(filepath)
        transmission_matrix = CSV.read(filepath, DataFrame, types=Complex{Float64})
        return [chebinterp(collect(transmission_matrix[i,:]), php.pillar_width_lb, php.pillar_width_ub) for i in 1:php.freq_order + 1]
    else
        error("file does not exist: surrogate models have not yet been computed and saved")
    end
end

function plot_surrogate_models(php::PhysicsHyperParams)
    datetime = now()
    surrogates = load_surrogate_models(php)
    surrogate_label = get_surrogate_label(php)
    freq_chebpoints = get_freq_chebpoints(php)
    width_chebpoints = get_width_chebpoints(php)
    widths_linear = LinRange(php.pillar_width_lb, php.pillar_width_ub, 1000)
    f = figure(figsize=(8, 2.5*(php.freq_order+1) ))
    for i in eachindex(freq_chebpoints)
        freq = freq_chebpoints[i]
        surrogate_vals_linear = surrogates[i].(widths_linear)
        surrogate_vals_chebpoints = surrogates[i].(width_chebpoints)
        subplot(php.freq_order+1, 2, i*2 - 1)
        plot(widths_linear, mod.(angle.(surrogate_vals_linear),2*pi)/(2*pi), color = "tab:blue")
        plot(width_chebpoints, mod.(angle.(surrogate_vals_chebpoints),2*pi)/(2*pi), color = "tab:blue", ".")
        xlabel("Pillar width (unitless)")
        ylim([0,1])
        ylabel(L"\angle t / 2\pi")
        title(L"f = %$freq")

        subplot(php.freq_order+1, 2, i*2 )
        plot(widths_linear, abs.(surrogate_vals_linear).^2, color = "tab:orange")
        plot(width_chebpoints, abs.(surrogate_vals_chebpoints).^2, color = "tab:orange", ".")
        xlabel("Pillar width (unitless)")
        ylim([0,1])
        ylabel(L"|t^2|")
        title(L"f = %$freq")
    end
    tight_layout()
    savefig("plots/$(surrogate_label)_$(datetime).png")
    plotclose()
    return nothing
end