using End2EndThermalImg
using BenchmarkTools

php = PhysicsHyperParams(
    λlb_μm = 8.0, 
    λub_μm = 12.0, 
    freq_order = 2, 
    focal_length_μm = 20000.0,
    num_unit_cells = 2048,
    unit_cell_length_μm = 4.0,
    pillar_width_lb_μm = 1.8,
    pillar_width_ub_μm = 2.7,
    pillar_width_order = 64,
    pillar_height_μm = 10.0,
    pillar_material = "Si_no_absorption",
    substrate_height_μm = 300.0,
    substrate_material = "Si_no_absorption",
    nG = 10
)

End2EndThermalImg.compute_and_save_surrogate_transmission_matrix(php)
End2EndThermalImg.plot_surrogate_models(php)