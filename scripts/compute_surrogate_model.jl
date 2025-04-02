#=
using End2EndThermalImg
using BenchmarkTools

php = PhysicsHyperParams(
    λlb_μm = 8.0, 
    λub_μm = 12.0, 
    freq_order = 5, 
    focal_length_μm = 20000.0,
    num_unit_cells = 2048,
    unit_cell_length_μm = 4.0,
    pillar_width_lb_μm = 1.8,
    pillar_width_ub_μm = 2.7,
    pillar_width_order = 5,
    pillar_height_μm = 10.0,
    pillar_material = "Si_no_absorption",
    substrate_height_μm = 300.0,
    substrate_material = "Si_no_absorption",
    nG = 1000
)

End2EndThermalImg.compute_and_save_surrogate_transmission_matrix(php)
End2EndThermalImg.plot_surrogate_models(php)
=#

using Distributed
using Parameters
@everywhere using PythonCall
@everywhere using End2EndThermalImg

php = PhysicsHyperParams(
    λlb_μm = 8.0, 
    λub_μm = 12.0, 
    freq_order = 5, 
    focal_length_μm = 20000.0,
    num_unit_cells = 2048,
    unit_cell_length_μm = 4.0,
    pillar_width_lb_μm = 1.8,
    pillar_width_ub_μm = 2.7,
    pillar_width_order = 5,
    pillar_height_μm = 10.0,
    pillar_material = "Si_no_absorption",
    substrate_height_μm = 300.0,
    substrate_material = "Si_no_absorption",
    nG = 1000
)

@unpack pillar_height, unit_cell_length = php
get_pillar_ϵ = End2EndThermalImg.get_permittivity_function(php.pillar_material)
get_substrate_ϵ = End2EndThermalImg.get_permittivity_function(php.substrate_material)
freq_chebpoints = End2EndThermalImg.get_freq_chebpoints(php)
width_chebpoints = End2EndThermalImg.get_width_chebpoints(php)

freq = freq_chebpoints[1]
width = width_chebpoints[1]
λ_µm = End2EndThermalImg.convert_freq_unitless_to_λ_µm(freq, php)
pillar_ϵ = get_pillar_ϵ(λ_µm)
substrate_ϵ = get_substrate_ϵ(λ_µm)

@time pmap(1:4) do i
    println(myid())
    transmission = get_transmission(freq, width, pillar_height, pillar_ϵ, unit_cell_length, substrate_ϵ, php.nG)
end