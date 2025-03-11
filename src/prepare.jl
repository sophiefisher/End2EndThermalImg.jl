# TODO: should wavcen be a parameter of PhysicsHyperParams?
struct PhysicsHyperParams{FloatType <: AbstractFloat, IntType <: Integer}
    # bandwidth hyper params
    λlb_μm::FloatType # wavelength lower bound (in units of µm)
    λub_μm::FloatType # wavelength upper bound (in units of µm)
    freq_order::IntType # order of Chebyshev polynomial in the frequency; number of points is freq_order + 1
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

    # Computed unitless parameters, normalized by λ correponding to center freq
    λlb::FloatType # wavelength lower bound 
    λub::FloatType # wavelength upper bound 
    freqlb::FloatType # frequency lower bound 
    frequb::FloatType # frequency upper bound
    focal_length::FloatType # focal length 
    unit_cell_length::FloatType # side length of each square unit cell 
    pillar_width_lb::FloatType # square pillar width lower bound 
    pillar_width_ub::FloatType # square pillar width upper bound 
    pillar_height::FloatType # square pillar height in z 
    substrate_height::FloatType # substrate height in z
end

function PhysicsHyperParams(; 
    λlb_μm::FloatType,
    λub_μm::FloatType,
    freq_order::IntType,
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
    freqlb = 1 / λub
    frequb = 1 / λlb
    focal_length = focal_length_μm / wavcen
    unit_cell_length = unit_cell_length_μm / wavcen
    pillar_width_lb = pillar_width_lb_μm / wavcen
    pillar_width_ub = pillar_width_ub_μm / wavcen
    pillar_height = pillar_height_μm / wavcen
    substrate_height = substrate_height_μm / wavcen
    return PhysicsHyperParams{FloatType, IntType}(
        λlb_μm,
        λub_μm,
        freq_order,
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
        freqlb,
        frequb,
        focal_length,
        unit_cell_length,
        pillar_width_lb,
        pillar_width_ub,
        pillar_height,
        substrate_height
    )
end

@with_kw struct ImagingHyperParams{IntType <: Integer}
    objN::IntType # number of object pixels in the x and y directions 
    imgN::IntType # number of image pixels in the x and y directions
    binN::IntType # how much to bin each sensor pixel (binN x binN subpixels)
    sampleN::IntType # how many points to sample per subpixel (which has length unit_cell_length, i.e. the metasurface unit cell length)
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

function get_wavcen(λlb_μm, λub_μm)
    freq_center = ((1/λlb_μm) + (1/λub_μm)) / 2
    wavcen = round(1 / freq_center, sigdigits = 3)
    wavcen
end

get_wavcen(php::PhysicsHyperParams) = get_wavcen(php.λlb_μm, php.λub_μm)

function prepare_geoms(jhp::JobHyperParams)
    geoms_init_type = jhp.opthp.geoms_init_type
    @unpack num_unit_cells, pillar_width_lb, pillar_width_ub = jhp.php
    
    if geoms_init_type == "uniform"
        return fill((pillar_width_lb + pillar_width_ub)/2, num_unit_cells, num_unit_cells)
    end
end