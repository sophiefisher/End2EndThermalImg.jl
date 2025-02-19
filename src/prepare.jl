struct PhysicsHyperParams{FloatType <: AbstractFloat, IntType <: Integer}
    # bandwidth hyper params
    λlb_μm::FloatType # wavelength lower bound (in units of µm)
    λub_μm::FloatType # wavelength upper bound (in units of µm)
    f_order::IntType # order of Chebyshev polynomial in the frequency; number of points is f_order + 1
    # geometry and material hyper params
    focal_length_μm::FloatType # focal length (in units of µm)
    num_unit_cells::IntType # number of unit cells on one side of the square metasurface (total number of unit cells is num_unit_cells^2)
    unit_cell_length_μm::FloatType # side length of each square unit cell (in units of µm)
    pillar_width_lb_μm::FloatType # square pillar width lower bound (in units of µm)
    pillar_width_ub_μm::FloatType # square pillar width upper bound (in units of µm)
    pillar_width_order::IntType # order of Chebyshev polynomial in the pillar width; number of points is pillar_width_order + 1
    pillar_height_μm::FloatType # square pillar height in z (in units of µm)
    pillar_material::String # material of the pillars
    substrate_height_μm::FloatType # substrate height in z (in units of µm) # TODO: do I still need this anywhere if I assume infinite substrate?
    substrate_material::String # material of the substrate
    nG::IntType # number of fourier components in RCWA (truncation order)

    # Computed unitless parameters
    λlb::FloatType # wavelength lower bound (unitless, normalized by λ correponding to center freq)
    λub::FloatType # wavelength upper bound (unitless, normalized by λ correponding to center freq)
    focal_length::FloatType # focal length (unitless, normalized by λ correponding to center freq)
    unit_cell_length::FloatType # side length of each square unit cell (unitless, normalized by λ correponding to center freq)
    pillar_width_lb::FloatType # square pillar width lower bound (unitless, normalized by λ correponding to center freq)
    pillar_width_ub::FloatType # square pillar width upper bound (unitless, normalized by λ correponding to center freq)
    pillar_height::FloatType # square pillar height in z (unitless, normalized by λ correponding to center freq)
    substrate_height::FloatType # substrate height in z (unitless, normalized by λ correponding to center freq)
end

function PhysicsHyperParams(; 
    λlb_μm::FloatType,
    λub_μm::FloatType,
    f_order::IntType,
    focal_length_μm::FloatType,
    num_unit_cells::IntType,
    unit_cell_length_μm::FloatType,
    pillar_width_lb_μm::FloatType,
    pillar_width_ub_μm::FloatType,
    pillar_width_order::IntType,
    pillar_height_μm::FloatType,
    pillar_material::String,
    substrate_height_μm::FloatType,
    substrate_material::String,
    nG::IntType
) where {FloatType <: AbstractFloat, IntType <: Integer}
    wavcen = get_wavcen(λlb_μm, λub_μm)
    λlb = λlb_μm / wavcen
    λub = λub_μm / wavcen
    focal_length = focal_length_μm / wavcen
    unit_cell_length = unit_cell_length_μm / wavcen
    pillar_width_lb = pillar_width_lb_μm / wavcen
    pillar_width_ub = pillar_width_ub_μm / wavcen
    pillar_height = pillar_height_μm / wavcen
    substrate_height = substrate_height_μm / wavcen
    return PhysicsHyperParams{FloatType, IntType}(
        λlb_μm,
        λub_μm,
        f_order,
        focal_length_μm,
        num_unit_cells,
        unit_cell_length_μm,
        pillar_width_lb_μm,
        pillar_width_ub_μm,
        pillar_width_order,
        pillar_height_μm,
        pillar_material,
        substrate_height_μm,
        substrate_material,
        nG,
        λlb,
        λub,
        focal_length,
        unit_cell_length,
        pillar_width_lb,
        pillar_width_ub,
        pillar_height,
        substrate_height
    )
end

@with_kw struct ImagingHyperParams{IntType <: Integer}
    objL::IntType # number of object pixels in the x and y directions 
    imgL::IntType # number of image pixels in the x and y directions
    binL::IntType # how much to bin each sensor pixel (binL x binL subpixels)
end

@with_kw struct OptimizeHyperParams
    geoms_init_type::String # how to initialize the metasurface for the end-to-end
end

@with_kw struct ReconstructionHyperParams
    T_init_type::String # how to initialize the temperature map for reconstruction
end

@with_kw struct JobHyperParams
    php::PhysicsHyperParams
    imghp::ImagingHyperParams
    opthp::OptimizeHyperParams
    rechp::ReconstructionHyperParams
end

function get_wavcen(λlb, λub)
    freq_center = ((1/λlb) + (1/λub)) / 2
    wavcen = round(1 / freq_center, sigdigits = 3)
    wavcen
end

get_wavcen(php::PhysicsHyperParams) = get_wavcen(php.λlb, php.λub)

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

function prepare_geoms(jhp::JobHyperParams)
    geoms_init_type = jhp.opthp.geoms_init_type
    @unpack num_unit_cells = jhp.php
    php_unitless = get_php_unitless(jhp.php)
    @unpack pillar_width_lb, pillar_width_ub = php_unitless
    
    if geoms_init_type == "uniform"
        return fill((pillar_width_lb + pillar_width_ub)/2, num_unit_cells, num_unit_cells)
    end
end