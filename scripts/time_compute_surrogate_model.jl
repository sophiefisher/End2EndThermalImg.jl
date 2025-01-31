using End2EndThermalImg
using BenchmarkTools

php = PhysicsHyperParams(λlb = 8.0, 
    λub = 12.0, 
    f_order = 1, 
    focal_length = 20000.0,
    num_unit_cells = 2048,
    unit_cell_length = 4.0,
    pillar_width_lb = 1.8,
    pillar_width_ub = 2.7,
    pillar_width_order = 2,
    pillar_height = 10.0,
    pillar_material = "Si_no_absorption",
    substrate_height = 300.0,
    substrate_material = "Si_no_absorption",
    nG = 10
)

@btime End2EndThermalImg.compute_and_save_surrogate_transmission_matrix(php)