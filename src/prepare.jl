struct PhysicsHyperParams{FloatType <: AbstractFloat, IntType <: Integer}
    # bandwidth hyper params
    λlb::FloatType # wavelength lower bound (in units of µm)
    λub::FloatType # wavelength upper bound (in units of µm)
    f_order::IntType # order of Chebyshev polynomial in the frequency; number of points is f_order + 1
    # geometry and material hyper params
    focal_length::FloatType # focal length (in units of µm)
    num_unit_cells::IntType # number of unit cells on one side of the square metasurface (total number of unit cells is num_unit_cells^2)
    unit_cell_length::FloatType # side length of each square unit cell (in units of µm)
    pillar_width_lb::FloatType # square pillar width lower bound (in units of µm)
    pillar_width_ub::FloatType # square pillar width upper bound (in units of µm)
    pillar_width_order::IntType # order of Chebyshev polynomial in the pillar width; number of points is pillar_width_order + 1
    pillar_height::FloatType # square pillar height in z (in units of µm)
    pillar_material::String # material of the pillars
    substrate_height::FloatType # substrate height in z (in units of µm) # TODO: do I still need this anywhere if I assume infinite substrate?
    substrate_material::String # material of the substrate
end

function get_wavcen(php::PhysicsHyperParams)
    freq_center = ((1/php.λlb) + (1/php.λub))/2
    round(1 / freq_center,sigdigits=3)
end

function get_freq_bounds(php::PhysicsHyperParams)
    wavcen = get_wavcen(php)
    flb = wavcen / php.λub
    fub = wavcen / php.λlb
    (; flb, fub)
end

function get_php_unitless(php::PhysicsHyperParams)
    wavcen = get_wavcen(php)
    flb, fub = get_freq_bounds(php)
    focal_length = php.focal_length / wavcen
    unit_cell_length = php.unit_cell_length / wavcen
    pillar_width_lb = php.pillar_width_lb / wavcen
    pillar_width_ub = php.pillar_width_ub / wavcen
    pillar_height = php.pillar_height / wavcen
    substrate_height = php.substrate_height / wavcen
    (; flb, fub, focal_length, unit_cell_length, pillar_width_lb, 
    pillar_width_ub, pillar_height, substrate_height)
end
