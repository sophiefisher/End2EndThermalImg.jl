function get_permittivity_function(material)
    if material == "Si_no_absorption"
        filepath = "materials/Si_n_Shkondin.csv"
        data = CSV.read(filepath, DataFrame)
        λ_µm = data[:,:wl]
        n = data[:,:n]
        itp = linear_interpolation(λ_µm, n)
        return get_permittivity(λ_µm) = itp(λ_µm) # wavelength must be in microns
    end
end

function get_freq_chebpoints(php::PhysicsHyperParams)
    flb, fub = get_freq_bounds(php)
    # first element in the list is fub, last is flb
    chebpoints(php.f_order, flb, fub) 
end

function get_width_chebpoints(php::PhysicsHyperParams)
    php_ul = get_php_unitless(php)
    # first element in the list is ub, last is lb
    chebpoints(php.pillar_width_order, php_ul.pillar_width_lb, php_ul.pillar_width_ub)
end

function get_transmission(freq, pillar_width, pillar_height, pillar_epsilon, unit_cell_length, substrate_epsilon, nG)
    grcwa = pyimport("grcwa")
    numpy = pyimport("numpy")
    # nG is number of fourier components in RCWA (truncation order)
    L1 = pylist([unit_cell_length,0])
    L2 = pylist([0,unit_cell_length])
    N_grid_points = pyint(10000) # number of grid points for the patterned layers
    
    obj = grcwa.obj(pyint(nG), L1, L2, freq, 0.0, 0.0, verbose = 1 )
    obj.Add_LayerUniform(3.0, substrate_epsilon)
    obj.Add_LayerGrid(pillar_height, N_grid_points, N_grid_points)
    obj.Add_LayerUniform(3.0, 1.0)  
    obj.Init_Setup()

    epgrid = numpy.zeros((N_grid_points, N_grid_points))
    boundary = unit_cell_length / 2.0
    xarr = numpy.linspace(-boundary, boundary, N_grid_points)
    yarr = numpy.linspace(-boundary, boundary, N_grid_points)
    x, y = numpy.meshgrid(xarr,yarr,indexing=pystr("ij") )
    ind = numpy.logical_and(pyabs(x) < pillar_width / 2.0, pyabs(y) < pillar_width / 2.0)
    epgrid[ind] = pyfloat(1.0)

    epgrid = (pillar_epsilon - 1.0) * epgrid + 1.0
    obj.GridLayer_geteps(epgrid.flatten())

    obj.a0 = numpy.zeros(2*obj.nG, dtype = pytype(pycomplex(1.0 )))
    obj.a0[0] = pycomplex(1.0)
    obj.bN = numpy.zeros(2*obj.nG, dtype = pytype(pycomplex(1.0 )))

    aN,b0 = obj.GetAmplitudes(2, 0.1)
    transmission = aN[0]/obj.a0[0] * numpy.sqrt(numpy.sqrt(substrate_epsilon))
    pyconvert(Complex{Float64}, transmission)
end