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

# Noise model applied in image formation (make_noisy_image_from_3D[_distributed]): FixedNoise
# adds noise at a fixed absolute level (image .+ noise_scale .* noise, the original behavior);
# RelativeNoise scales it by the mean image value instead (image .+ mean(image) .* noise_scale
# .* noise), so noise_scale parameterizes a percentage of the mean rather than an absolute level.
abstract type AbstractNoiseType end

struct FixedNoise{FloatType <: AbstractFloat} <: AbstractNoiseType
    noise_scale::FloatType
end

struct RelativeNoise{FloatType <: AbstractFloat} <: AbstractNoiseType
    noise_scale::FloatType
end

apply_noise(nt::FixedNoise, image, noise) = image .+ nt.noise_scale .* noise
apply_noise(nt::RelativeNoise, image, noise) = image .+ mean(image) .* nt.noise_scale .* noise

struct ImagingHyperParams{FloatType <: AbstractFloat, IntType <: Integer}
    objN::IntType # number of object pixels in the x and y directions
    imgN::IntType # number of image pixels in the x and y directions
    binN::IntType # how much to bin each sensor pixel (binN x binN subpixels)
    sampleN::IntType # how many points to sample per subpixel (which has length unit_cell_length, i.e. the metasurface unit cell length)
    PSF_zlb_μm::FloatType
    PSF_zub_μm::FloatType
    PSF_zlen::IntType
    smoothness_order::FloatType # smoothness order of the discretized δ function
    object_type::AbstractObjectType # type of object to generate
    noise_type::AbstractNoiseType # noise model (and scale) applied in image formation

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
    smoothness_order::FloatType,
    object_type::AbstractObjectType,
    noise_type::AbstractNoiseType,
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
        smoothness_order,
        object_type,
        noise_type,
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
    geoms_init_path::String = "" # only used when geoms_init_type == "loaded" (see load_geoms, src/inversedesign.jl)
    maxeval::Int = 1000 # number of outer (geoms, α, β) optimization iterations; the only stopping criterion for now
    learning_rate::Float64 = 1e-3 # Adam learning rate used by train_outer_optimization_distributed (src/optimize.jl)
    save_every::Int = 10 # checkpoint interval (in outer iterations) for run_and_save_outer_optimization_distributed;
                         # the last iteration is always saved too, regardless of this value
end

struct ReconstructionHyperParams{FloatType <: AbstractFloat}
    object_init_type::String # how to initialize the object for reconstruction
    T_background::FloatType
    z_middle_μm::FloatType
    algorithm::Symbol # NLopt algorithm used by reconstruct_Tmap_and_zmap[_distributed], e.g. :LD_LBFGS, :LD_MMA
    α_init::Union{FloatType, Symbol} # initial regularization weights for the outer optimization
    β_init::Union{FloatType, Symbol} # (eventually optimized jointly with geoms). Either a literal
                                      # value, or :auto, which resolve_initial_α_β
                                      # (src/tune_regularization.jl) resolves via its own sweep.
    maxeval::Int # NLopt stopping criteria used by reconstruct_Tmap_and_zmap[_distributed]
    xtol_rel::FloatType
    ftol_rel::FloatType
    ftol_abs::FloatType

    # Computed parameters
    z_middle::FloatType
end

function ReconstructionHyperParams(;
    object_init_type::String,
    T_background::FloatType,
    z_middle_μm::FloatType,
    algorithm::Symbol = :LD_LBFGS,
    α_init::Union{FloatType, Symbol},
    β_init::Union{FloatType, Symbol},
    maxeval::Int = 5000,
    xtol_rel::FloatType = 1e-8,
    ftol_rel::FloatType = 0.0,
    ftol_abs::FloatType = 0.0,
    php::PhysicsHyperParams
) where {FloatType <: AbstractFloat}
    wavcen = php.wavcen
    z_middle = z_middle_μm / wavcen

    return ReconstructionHyperParams{FloatType}(
        object_init_type,
        T_background,
        z_middle_μm,
        algorithm,
        α_init,
        β_init,
        maxeval,
        xtol_rel,
        ftol_rel,
        ftol_abs,
        z_middle
    )
end


@with_kw struct JobHyperParams
    php::PhysicsHyperParams
    imghp::ImagingHyperParams
    opthp::OptimizeHyperParams
    rechp::ReconstructionHyperParams
end

# Recursively writes "label = value" for each field of x, one per line, indenting into nested
# struct fields (e.g. imghp.object_type/imghp.noise_type) rather than falling back to their
# unreadable single-line default show. Generic over any struct via fieldnames/getfield, so it
# stays correct automatically as fields are added/removed from the hyperparameter types, rather
# than needing a hand-written formatter per type.
function write_param_lines(io, label, x; indent = "")
    if isstructtype(typeof(x)) && !(x isa Union{AbstractString, Symbol, Number, Bool})
        println(io, indent, label, " (", typeof(x), "):")
        for field in fieldnames(typeof(x))
            write_param_lines(io, field, getfield(x, field); indent = indent * "  ")
        end
    else
        println(io, indent, label, " = ", x)
    end
end

# Writes a human-readable summary of this run's hyperparameters to filepath, for visual
# inspection in a run folder without loading results.jld2 in Julia -- the same role history.csv
# plays for the optimization history (src/optimize.jl). Not a replacement for results.jld2 (the
# programmatic source of truth, needed for reloading/analysis), just a convenience view.
function write_params_summary(filepath, php, imghp, opthp, rechp, system_params, PSF_scale)
    open(filepath, "w") do io
        write_param_lines(io, :php, php)
        println(io)
        write_param_lines(io, :imghp, imghp)
        println(io)
        write_param_lines(io, :opthp, opthp)
        println(io)
        write_param_lines(io, :rechp, rechp)
        println(io)
        write_param_lines(io, :system_params, system_params)
        for w in system_param_warnings(system_params)
            println(io, "WARNING: ", w)
        end
        println(io)
        write_param_lines(io, :PSF_scale, PSF_scale)
    end
end

function get_wavcen(λlb_μm, λub_μm)
    freq_center = ((1/λlb_μm) + (1/λub_μm)) / 2
    wavcen = round(1 / freq_center, sigdigits = 3)
    wavcen
end

# imghp: only required when opthp.geoms_init_type == "monochromatic" -- kept as an optional
# keyword (not a required positional argument) so every existing call site (which never uses that
# mode) keeps working unchanged.
function initialize_geoms(php::PhysicsHyperParams, opthp::OptimizeHyperParams; imghp = nothing)
    opthp.geoms_init_type == "loaded" && return load_geoms(php, opthp.geoms_init_path)
    if opthp.geoms_init_type == "monochromatic"
        imghp === nothing && error("initialize_geoms: geoms_init_type == \"monochromatic\" requires imghp (pass it as a keyword argument)")
        return find_or_design_monochromatic_lens(php, imghp)
    end
    initialize_geoms(php, opthp.geoms_init_type)
end

function initialize_geoms(php::PhysicsHyperParams, geoms_init_type)
    @unpack num_unit_cells, pillar_width_lb, pillar_width_ub = php

    if geoms_init_type == "uniform"
        return fill((pillar_width_lb + pillar_width_ub)/2, num_unit_cells, num_unit_cells)
    else
        error("Unrecognized geoms_init_type: $geoms_init_type")
    end
end

function initialize_object(imghp::ImagingHyperParams, rechp::ReconstructionHyperParams)
    @unpack objN = imghp
    @unpack object_init_type, T_background = rechp
    PSF_zcoords = get_PSF_zcoords(imghp)

    # TODO: might want to base center_zcoord off of the object z limits instead of the PSF limits
    # but then might have to modify code based on object_type
    if object_init_type == "uniform"
        Tmap = fill(T_background, objN, objN)
        center_zcoord_idx = cld(length(PSF_zcoords), 2)
        center_zcoord = PSF_zcoords[center_zcoord_idx]
        zmap = fill(center_zcoord, objN, objN)
        return (; Tmap, zmap)
    elseif object_init_type == "uniform_offset"
        # same as "uniform", except zmap starts halfway between the center PSF_zcoords node and
        # its upper neighbor rather than exactly on the node -- get_C_interp_3D's depth-interpolation
        # weights have zero derivative wrt z exactly at a node (the discretized delta function is
        # flat at both ends of its support), so a reconstruction started there gets an identically
        # zero z-gradient and z never moves; the midpoint between nodes is where that derivative is
        # largest, giving z the strongest starting gradient
        Tmap = fill(T_background, objN, objN)
        center_zcoord_idx = cld(length(PSF_zcoords), 2)
        zmap = fill(PSF_zcoords[center_zcoord_idx] + imghp.PSF_Δz / 2, objN, objN)
        return (; Tmap, zmap)
    else
        error("unknown object_init_type: $object_init_type")
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

# (zlb_μm, zub_μm) spanning where this object_type actually places z values, used by
# compute_system_parameters to check that this lies within imghp.PSF_zlb_μm/PSF_zub_μm. For
# GaussianObject, despite the field names (an offset and a peak, not literal bounds), the sampled
# zmap = zlb + (zub-zlb)*exp(-...) is a convex combination of the two (exp(...) ∈ (0,1]), so it
# never leaves [zlb_μm, zub_μm] either -- same bound as the uniform case, for a different reason.
get_object_zextent(object_type::GaussianObject) = (object_type.zlb_μm, object_type.zub_μm)
get_object_zextent(object_type::UniformlyRandomObject) = (object_type.zlb_μm, object_type.zub_μm)
get_object_zextent(object_type::UniformlyRandomTFixedDepthObject) = (object_type.z_μm, object_type.z_μm)

# (Tlb, Tub) spanning where this object_type actually places T values -- one generic method (not
# dispatched per subtype like get_object_zextent) since every AbstractObjectType subtype stores
# this the same way, as a literal (Tlb, Tub) pair with the same convex-combination argument as
# get_object_zextent's GaussianObject case applying here too.
get_object_Textent(object_type::AbstractObjectType) = (object_type.Tlb, object_type.Tub)

# Deterministic seed derived from (php, imghp) -- the same values always hash to the same seed
# (confirmed: Julia's default struct hash is fieldwise/value-based, not objectid-based, so this is
# stable across separate Julia processes/job submissions, not just within one).
tuning_seed(php::PhysicsHyperParams, imghp::ImagingHyperParams) = hash((php, imghp))

# Runs f() with the global RNG temporarily seeded to tuning_seed(php, imghp), then restores the RNG
# to whatever state it was in before this call -- so f()'s draws become a reproducible function of
# (php, imghp), while everything that runs AFTER this call (in particular
# train_outer_optimization_distributed's own per-iteration object/noise draws, src/optimize.jl) is
# left exactly as if this call had never happened, not merely "seeded here and left seeded from
# then on" (which would otherwise make the outer loop's randomness depend on incidental details
# like whether find_or_resolve_initial_α_β's cache was a hit or a miss). Confirmed empirically: a
# draw taken right after this call is bit-for-bit identical to the same draw taken with no call to
# this function in between at all.
#
# Used to make choosing PSF_scale's reference object (scripts/endtoend_opt.jl) and
# resolve_initial_α_β's object/noise draw (src/tune_regularization.jl) reproducible -- which is
# also what lets find_or_resolve_initial_α_β's cache (src/tune_regularization.jl) actually hit
# across separate job submissions with the same hyperparameters, rather than comparing against
# results from an unreproducible random draw.
function with_tuning_rng(f, php::PhysicsHyperParams, imghp::ImagingHyperParams)
    saved_rng = copy(Random.default_rng())
    Random.seed!(tuning_seed(php, imghp))
    try
        f()
    finally
        copy!(Random.default_rng(), saved_rng)
    end
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
# Returns a NamedTuple of system quantities derived from php/imghp (not stored on either struct
# itself), and prints it via write_param_lines in the same "label (Type):\n  field = value" style
# as the job hyperparameters themselves, for a consistent log format. The returned NamedTuple is
# also threaded through to write_params_summary by callers that want it recorded in params.txt
# (e.g. run_and_save_outer_optimization_distributed, src/optimize.jl).
function compute_system_parameters(php::PhysicsHyperParams, imghp::ImagingHyperParams, rechp::ReconstructionHyperParams)
    image_pixel_size_μm = imghp.binN * php.unit_cell_length_μm
    image_size_μm = image_pixel_size_μm * imghp.imgN
    metasurface_size_μm = php.unit_cell_length_μm * php.num_unit_cells
    metasurface_size_cm = metasurface_size_μm / 1e4
    NA = sin(atan( php.unit_cell_length_μm * php.num_unit_cells / (2 * php.focal_length_μm) ))
    f_number = php.focal_length_μm / metasurface_size_μm
    FOV_degrees = 2 * atan(image_size_μm / (2 * php.focal_length_μm)) * 180 / π

    # Diffraction limit (λ / (2·NA)) at the band edges and center -- object-side quantities like
    # object_pixel_size/object_size/object_angle (computed in an older version of this package)
    # are deliberately NOT included here: they depend on a single scalar depth, which this package
    # no longer has (imghp has a z-range for PSF sampling instead, generally distinct from
    # imghp.object_type's own z-extent), so porting them isn't a direct translation.
    λmid_μm = (php.λlb_μm + php.λub_μm) / 2
    difflim_λlb_μm = php.λlb_μm / (2 * NA)
    difflim_λmid_μm = λmid_μm / (2 * NA)
    difflim_λub_μm = php.λub_μm / (2 * NA)

    # Overdetermination check: object unknowns are Tmap + zmap, each objN×objN; measurements are
    # the single imgN×imgN sensor image (sampleN/PSF_zlen are internal integration/sampling
    # parameters, not separate acquisitions, so they don't factor in here).
    n_object_pixels = 2 * imghp.objN^2
    n_image_pixels = imghp.imgN^2
    overdetermination_ratio = n_image_pixels / n_object_pixels
    is_overdetermined = overdetermination_ratio > 1

    # Plane-wave / far-field check: each object point illuminates the metasurface as an exact
    # spherical wave (see efield_point_source in forward.jl); the worst case (largest deviation
    # from a plane wave) is the closest point in z and the shortest wavelength. max_phase_error_rad
    # is the paraxial/Fresnel quadratic phase term, valid only when paraxial_ratio (= tan of the
    # half-angle the aperture subtends from the object point) is itself small -- is_paraxial gates
    # whether max_phase_error_rad/is_farfield can be trusted at all.
    z_min_μm = min(abs(imghp.PSF_zlb_μm), abs(imghp.PSF_zub_μm))
    paraxial_ratio = (metasurface_size_μm / 2) / z_min_μm
    is_paraxial = paraxial_ratio < 0.1
    max_phase_error_rad = π * metasurface_size_μm^2 / (4 * php.λlb_μm * z_min_μm)
    is_farfield = max_phase_error_rad < π / 8

    # Object-extent check: imghp.object_type's own z-extent (get_object_zextent, dispatched per
    # object_type) should lie entirely within imghp.PSF_zlb_μm/PSF_zub_μm -- otherwise part of the
    # object is placed outside the z-range the PSFs are actually sampled over.
    object_zlb_μm, object_zub_μm = get_object_zextent(imghp.object_type)
    is_object_within_PSF_bounds = object_zlb_μm >= imghp.PSF_zlb_μm && object_zub_μm <= imghp.PSF_zub_μm

    # Reconstruction-initialization checks: rechp.T_background/z_middle_μm should lie within (and,
    # ideally but not strictly necessarily, at the midpoint of) the object's own T/z extent -- the
    # midpoint is the natural unbiased starting guess, so being off-center is worth flagging even
    # though it isn't necessarily wrong.
    object_Tlb, object_Tub = get_object_Textent(imghp.object_type)
    is_T_background_within_object_Textent = object_Tlb <= rechp.T_background <= object_Tub
    is_T_background_at_midpoint = isapprox(rechp.T_background, (object_Tlb + object_Tub) / 2)
    is_z_middle_within_object_zextent = object_zlb_μm <= rechp.z_middle_μm <= object_zub_μm
    is_z_middle_at_midpoint = isapprox(rechp.z_middle_μm, (object_zlb_μm + object_zub_μm) / 2)

    system_params = (; image_pixel_size_μm, image_size_μm, metasurface_size_μm, metasurface_size_cm,
                        NA, f_number, FOV_degrees, λmid_μm, difflim_λlb_μm, difflim_λmid_μm, difflim_λub_μm,
                        n_object_pixels, n_image_pixels, overdetermination_ratio, is_overdetermined,
                        z_min_μm, paraxial_ratio, is_paraxial, max_phase_error_rad, is_farfield,
                        object_zlb_μm, object_zub_μm, is_object_within_PSF_bounds,
                        object_Tlb, object_Tub, is_T_background_within_object_Textent, is_T_background_at_midpoint,
                        is_z_middle_within_object_zextent, is_z_middle_at_midpoint)
    write_param_lines(stdout, :system_params, system_params)
    for w in system_param_warnings(system_params)
        @warn w
    end
    system_params
end

# Human-readable messages for any system_params condition (computed by compute_system_parameters)
# that doesn't hold -- shared between the log (compute_system_parameters, via @warn) and params.txt
# (write_params_summary, via plain println) so the conditions/wording can't drift between the two.
function system_param_warnings(system_params)
    warnings = String[]
    system_params.is_overdetermined || push!(warnings,
        "system is not overdetermined (overdetermination_ratio = $(system_params.overdetermination_ratio) ≤ 1): more object unknowns than sensor measurements")
    system_params.is_paraxial || push!(warnings,
        "paraxial approximation does not hold (paraxial_ratio = $(system_params.paraxial_ratio) ≥ 0.1): max_phase_error_rad/is_farfield may not be reliable")
    system_params.is_farfield || push!(warnings,
        "object illumination is not sufficiently far-field/plane-wave-like (max_phase_error_rad = $(system_params.max_phase_error_rad) ≥ π/8)")
    system_params.is_object_within_PSF_bounds || push!(warnings,
        "object z-extent ($(system_params.object_zlb_μm), $(system_params.object_zub_μm)) μm is not entirely within imghp.PSF_zlb_μm/PSF_zub_μm")
    system_params.is_T_background_within_object_Textent || push!(warnings,
        "rechp.T_background is not within the object's T-extent ($(system_params.object_Tlb), $(system_params.object_Tub))")
    system_params.is_T_background_at_midpoint || push!(warnings,
        "rechp.T_background is not at the midpoint of the object's T-extent ($(system_params.object_Tlb), $(system_params.object_Tub)) -- not strictly required, but is the natural unbiased starting guess")
    system_params.is_z_middle_within_object_zextent || push!(warnings,
        "rechp.z_middle_μm is not within the object's z-extent ($(system_params.object_zlb_μm), $(system_params.object_zub_μm)) μm")
    system_params.is_z_middle_at_midpoint || push!(warnings,
        "rechp.z_middle_μm is not at the midpoint of the object's z-extent ($(system_params.object_zlb_μm), $(system_params.object_zub_μm)) μm -- not strictly required, but is the natural unbiased starting guess")
    warnings
end