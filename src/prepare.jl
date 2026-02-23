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

    # Computed parameters
    wavcen::FloatType # wavelength corresponding to center frequency
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
        wavcen,
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

abstract type AbstractObjectType end

struct ImagingHyperParams{FloatType <: AbstractFloat, IntType <: Integer}
    objN::IntType # number of object pixels in the x and y directions 
    imgN::IntType # number of image pixels in the x and y directions
    binN::IntType # how much to bin each sensor pixel (binN x binN subpixels)
    sampleN::IntType # how many points to sample per subpixel (which has length unit_cell_length, i.e. the metasurface unit cell length)
    PSF_zlb_μm::FloatType
    PSF_zub_μm::FloatType
    PSF_zlen::IntType
    PSF_scale::FloatType
    smoothness_order::FloatType # smoothness order of the discretized δ function
    object_type::AbstractObjectType # type of object to generate
    noise_level # noise percentage of the mean image

    # Computed parameters
    PSF_zlb::FloatType
    PSF_zub::FloatType
    PSF_Δz_μm::FloatType
    PSF_Δz::FloatType
end

function ImagingHyperParams(; 
    objN::IntType,
    imgN::IntType,
    binN::IntType,
    sampleN::IntType,
    PSF_zlb_μm::FloatType,
    PSF_zub_μm::FloatType,
    PSF_zlen::IntType,
    PSF_scale::FloatType,
    smoothness_order::FloatType,
    object_type::AbstractObjectType,
    noise_level::FloatType,
    php::PhysicsHyperParams
) where {FloatType <: AbstractFloat, IntType <: Integer}
    wavcen = php.wavcen
    PSF_zlb = PSF_zlb_μm / wavcen
    PSF_zub = PSF_zub_μm / wavcen
    PSF_Δz_μm = (PSF_zub_μm - PSF_zlb_μm) / (PSF_zlen - 1)
    PSF_Δz = PSF_Δz_μm / wavcen
    
    return ImagingHyperParams{FloatType, IntType}(
        objN,
        imgN,
        binN,
        sampleN,
        PSF_zlb_μm,
        PSF_zub_μm,
        PSF_zlen,
        PSF_scale,
        smoothness_order,
        object_type,
        noise_level,
        PSF_zlb,
        PSF_zub,
        PSF_Δz_μm,
        PSF_Δz
    )
end

struct GaussianObject{FloatType <: AbstractFloat} <: AbstractObjectType
    Tlb::FloatType # offset of the gaussian in temperature (in units of Kelvin)
    Tub::FloatType # peak of the gaussian in temperature (in units of Kelvin)
    std_dev_T::FloatType # standard deviation of the gaussian in temperature (in units of Kelvin)
    zlb_μm::FloatType # offset of the gaussian in z (in units of μm)
    zub_μm::FloatType # peak of the gaussian in z (in units of μm)
    std_dev_z::FloatType # standard deviation of the gaussian in z (in units of μm)

    # Computed parameters
    zlb::FloatType # offset of the gaussian in z (unitless)
    zub::FloatType # peak of the gaussian in z (unitless)
end

function GaussianObject(; 
    Tlb::FloatType, 
    Tub::FloatType,
    std_dev_T::FloatType,
    zlb_μm::FloatType,
    zub_μm::FloatType,
    std_dev_z::FloatType,
    php::PhysicsHyperParams
) where {FloatType <: AbstractFloat}
    wavcen = php.wavcen
    zlb = zlb_μm / wavcen
    zub = zub_μm / wavcen

    return GaussianObject{FloatType}(
        Tlb,
        Tub,
        std_dev_T,
        zlb_μm,
        zub_μm,
        std_dev_z,
        zlb,
        zub
    )
end

# uniformly random Tmap and uniformly random depth map (within bounds)
struct UniformlyRandomObject{FloatType <: AbstractFloat, IntType <: Integer} <: AbstractObjectType
    Tlb::FloatType # lower bound of the temperature (in units of Kelvin)
    Tub::FloatType # upper bound of the temperature (in units of Kelvin)
    zlb_μm::FloatType # lower bound z coordinate of the object (assumes the metasurface is at z = 0, so this should be negative) (in units of μm)
    zub_μm::FloatType # upper bound z coordinate of the object (assumes the metasurface is at z = 0, so this should be negative) (in units of μm)
    zlen::IntType # discretization width of the z range (unlike for temperature, this is not eps() since the range is too large)

    # Computed parameters
    zlb::FloatType
    zub::FloatType
    Δz_μm::FloatType # TODO: is this actually used anywhere? # Δz of the z range (in units of µm)
    Δz::FloatType # TODO: is this actually used anywhere? # Δz of the z range (unitless)
end

function UniformlyRandomObject(; 
    Tlb::FloatType, 
    Tub::FloatType,
    zlb_μm::FloatType,
    zub_μm::FloatType,
    zlen::IntType, 
    php::PhysicsHyperParams
) where {FloatType <: AbstractFloat, IntType <: Integer}
    wavcen = php.wavcen
    zlb = zlb_μm / wavcen
    zub = zub_μm / wavcen
    Δz_μm = (zub_μm - zlb_μm) / (zlen - 1)
    Δz = Δz_μm / wavcen
    
    return UniformlyRandomObject{FloatType, IntType}(
        Tlb,
        Tub,
        zlb_μm,
        zub_μm,
        zlen,
        zlb,
        zub,
        zub_μm,
        Δz
    )
end

# uniformly random Tmap (within bounds) and at a fixed depth
struct UniformlyRandomTFixedDepthObject{FloatType <: AbstractFloat} <: AbstractObjectType
    Tlb::FloatType # lower bound of the temperature (in units of Kelvin)
    Tub::FloatType # upper bound of the temperature (in units of Kelvin)
    z_μm::FloatType # z coordinate of the object (assumes the metasurface is at z = 0, so this should be negative) (in units of μm)

    # Computed parameters
    z::FloatType
end

function UniformlyRandomTFixedDepthObject(; 
    Tlb::FloatType, 
    Tub::FloatType,
    z_μm::FloatType,
    php::PhysicsHyperParams
) where {FloatType <: AbstractFloat}
    wavcen = php.wavcen
    z = z_μm / wavcen
    
    return UniformlyRandomTFixedDepthObject{FloatType}(
        Tlb,
        Tub,
        z_μm,
        z
    )
end

@with_kw struct OptimizeHyperParams
    geoms_init_type::String # how to initialize the metasurface for the end-to-end
end

struct ReconstructionHyperParams{FloatType <: AbstractFloat}
    object_init_type::String # how to initialize the object for reconstruction
    T_background::FloatType
    z_middle_μm::FloatType

    # Computed parameters
    z_middle::FloatType
end

function ReconstructionHyperParams(; 
    object_init_type::String,
    T_background::FloatType,
    z_middle_μm::FloatType,
    php::PhysicsHyperParams
) where {FloatType <: AbstractFloat}
    wavcen = php.wavcen
    z_middle = z_middle_μm / wavcen
    
    return ReconstructionHyperParams{FloatType}(
        object_init_type,
        T_background,
        z_middle_μm,
        z_middle
    )
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

function initialize_geoms(php::PhysicsHyperParams, opthp::OptimizeHyperParams)
    geoms_init_type = opthp.geoms_init_type
    initialize_geoms(php, geoms_init_type)
end

function initialize_geoms(php::PhysicsHyperParams, geoms_init_type)
    @unpack num_unit_cells, pillar_width_lb, pillar_width_ub = php
    
    if geoms_init_type == "uniform"
        return fill((pillar_width_lb + pillar_width_ub)/2, num_unit_cells, num_unit_cells)
    end
end

function initialize_object(imghp::ImagingHyperParams, rechp::ReconstructionHyperParams)
    @unpack objN = imghp
    @unpack object_init_type, T_background = rechp
    PSF_zcoords = get_PSF_zcoords(imghp)

    if object_init_type == "uniform"
        Tmap = fill(T_background, objN, objN)
        center_zcoord_idx = cld(length(PSF_zcoords), 2)
        center_zcoord = PSF_zcoords[center_zcoord_idx]
        zmap = fill(center_zcoord, objN, objN)
        return (; Tmap, zmap)
    end
end

gaussian_2D(x, y, offset, amplitude, std_dev) = offset + amplitude*exp(-(x^2 + y^2) / (2*std_dev^2))

function get_object(object_type::GaussianObject, imghp::ImagingHyperParams)
    N = imghp.objN
    coords = collect(-div(N,2):div(N-1,2))
    X = repeat(coords, 1, N)
    Y = repeat(coords', N, 1)
    Tmap = gaussian_2D.(X, Y, object_type.Tlb, (object_type.Tub - object_type.Tlb), object_type.std_dev_T)
    zmap = gaussian_2D.(X, Y, object_type.zlb, (object_type.zub - object_type.zlb),  object_type.std_dev_z)
    (; Tmap, zmap)
end

get_object_zrange(object_type::UniformlyRandomObject) = LinRange(object_type.zlb, object_type.zub, object_type.zlen)

function get_object(object_type::UniformlyRandomObject, imghp::ImagingHyperParams)
    Tmap = rand(object_type.Tlb:eps():object_type.Tub, imghp.objN, imghp.objN)
    object_zrange = get_object_zrange(object_type)
    zmap = rand(object_zrange, imghp.objN, imghp.objN)
    (; Tmap, zmap)
end

function get_object(object_type::UniformlyRandomTFixedDepthObject, imghp::ImagingHyperParams)
    Tmap = rand(object_type.Tlb:eps():object_type.Tub, imghp.objN, imghp.objN)
    zmap = fill(object_type.z, imghp.objN, imghp.objN)
    (; Tmap, zmap)
end

function get_object(imghp::ImagingHyperParams)
    get_object(imghp.object_type, imghp)
end

function get_PSF_zcoords(imghp::ImagingHyperParams)
    LinRange(imghp.PSF_zlb, imghp.PSF_zub, imghp.PSF_zlen)
end

# weights are symmetric, so don't need to reverse them
function get_clenshaw_curtis_quadrature_weights(php::PhysicsHyperParams)
    weights = ClenshawCurtisQuadrature(php.freq_order + 1).weights .* (php.frequb .- php.freqlb)
    weights
end

unflatten_square_matrix(matrix_flat) = reshape(matrix_flat, round(Int, sqrt(length(matrix_flat))), round(Int, sqrt(length(matrix_flat))))

flatten_object(object) = [object.Tmap[:]; object.zmap[:]]

function unflatten_object(object_flat)
    N = length(object_flat) ÷ 2  
    Tmap_flat = object_flat[1:N]
    Tmap = unflatten_square_matrix(Tmap_flat)
    zmap_flat = object_flat[N+1:end]
    zmap = unflatten_square_matrix(zmap_flat)
    (; Tmap, zmap)
end

prepare_noise_buf(imghp) = Array{Float64}(undef, imghp.imgN, imghp.imgN)

prepare_image_buf(imghp) = Array{Float64}(undef, imghp.imgN, imghp.imgN)

# TODO: add object pixel sizes for different depths
function compute_system_parameters(php::PhysicsHyperParams, imghp::ImagingHyperParams)
    @info "Computing system parameters"

    image_pixel_size = imghp.binN * php.unit_cell_length_μm
    @info "Image pixel size = $(round(image_pixel_size,digits=4)) μm"

    metasurface_size = php.unit_cell_length_μm * php.num_unit_cells
    @info "Metasurface size = $(round(metasurface_size,digits=4)) μm [$( round(metasurface_size / 1e4,digits=4)) cm]"

    NA = sin(atan( php.unit_cell_length_μm * php.num_unit_cells / (2 * php.focal_length_μm) ))
    @info "NA: $( round(NA,digits=4 ))"
end