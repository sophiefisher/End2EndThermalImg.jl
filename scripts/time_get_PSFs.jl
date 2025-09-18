using End2EndThermalImg
using BenchmarkTools

php = PhysicsHyperParams(
    λlb_μm = 8.0, 
    λub_μm = 12.0, 
    freq_order = 20, 
    focal_length_μm = 20000.0,
    num_unit_cells = 2048,
    unit_cell_length_μm = 4.0,
    pillar_width_lb_μm = 1.8,
    pillar_width_ub_μm = 2.7,
    pillar_width_order = 80,
    pillar_height_μm = 10.0,
    pillar_material = "Si_no_absorption",
    substrate_height_μm = 300.0,
    substrate_material = "Si_no_absorption",
    nG = 1000
)

imghp = ImagingHyperParams(
    objN = 32,
    imgN = 64,
    binN = 1,
    sampleN = 1,
    PSF_zlb_μm = -1.0e9,
    PSF_zub_μm = -9.5e8,
    PSF_zlen = 10,
    smoothness_order = 3.0,
    object_type = UniformlyRandomObject(Tlb = 136.15, Tub = 450.15, zlb_μm = -1.0e9, zub_μm =  -9.51e8, zlen = 100, php = php),
    php = php
)

opthp = OptimizeHyperParams(
    geoms_init_type = "uniform"
)

geoms = End2EndThermalImg.initialize_geoms(php, opthp)
PSF_zcoords = End2EndThermalImg.get_PSF_zcoords(imghp)
freqs = End2EndThermalImg.get_freq_chebpoints(php)
surrogates = End2EndThermalImg.load_surrogate_models(php)
incidents = End2EndThermalImg.get_incident_fields(freqs, PSF_zcoords, php);
n2f_kernels = End2EndThermalImg.get_n2f_kernels(freqs, php, imghp);
# @btime n2f_kernels = End2EndThermalImg.get_n2f_kernels(freqs, php, imghp);

# get PSF at fixed frequency, depth 
iF = 1
iZ = 1
@btime End2EndThermalImg.get_PSF_at_freq_and_z(freqs[iF], incidents[iF, iZ], surrogates[iF], geoms, n2f_kernels[iF], php, imghp);

# get PSFs at fixed depth
@btime End2EndThermalImg.get_PSFs_at_z(freqs, incidents[:, iZ], surrogates, geoms, n2f_kernels, php::PhysicsHyperParams, imghp::ImagingHyperParams);
