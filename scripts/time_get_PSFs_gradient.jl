using End2EndThermalImg
using BenchmarkTools
using Zygote

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

rechp = ReconstructionHyperParams(
    T_init_type = "uniform"
)

jhp = JobHyperParams(
    php = php,
    imghp = imghp,
    opthp = opthp,
    rechp = rechp
)

object = End2EndThermalImg.get_object(imghp)
geoms = End2EndThermalImg.initialize_geoms(php, opthp)
PSF_zcoords = End2EndThermalImg.get_PSF_zcoords(imghp)

B = End2EndThermalImg.get_black_body_spectrum(object.Tmap, php)
freqs = End2EndThermalImg.get_freq_chebpoints(php)
surrogates = End2EndThermalImg.load_surrogate_models(php)
incidents = End2EndThermalImg.get_incident_fields(freqs, PSF_zcoords, php);
n2f_kernels = End2EndThermalImg.get_n2f_kernels(freqs, php, imghp);

@info "Timing get_PSFs"
@time sum(sum(End2EndThermalImg.get_PSFs(freqs, incidents, surrogates, geoms, n2f_kernels, php, imghp)))
@time sum(sum(End2EndThermalImg.get_PSFs(freqs, incidents, surrogates, geoms, n2f_kernels, php, imghp)))

@info "Timing Zygote.gradient of get_PSFs with respect to geoms"
@time Zygote.gradient(geoms -> sum(sum(End2EndThermalImg.get_PSFs(freqs, incidents, surrogates, geoms, n2f_kernels, php, imghp))), geoms)
@time Zygote.gradient(geoms -> sum(sum(End2EndThermalImg.get_PSFs(freqs, incidents, surrogates, geoms, n2f_kernels, php, imghp))), geoms)
