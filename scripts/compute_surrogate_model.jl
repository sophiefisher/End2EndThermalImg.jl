using End2EndThermalImg
using BenchmarkTools

# freq_order is overridable via the FREQ_ORDER environment variable (set per-submission with
# `sbatch --export=ALL,FREQ_ORDER=<n> cluster_scripts/run_compute_surrogate_model.sh`), to compute
# surrogate models at several frequency orders (for the frequency-convergence study) without
# maintaining separate copies of this script. All other php fields match the production
# hyperparameters used in scripts/time_pipeline.jl and scripts/endtoend_opt.jl, so the resulting
# surrogate models are directly comparable/interchangeable with those runs.
freq_order = parse(Int, get(ENV, "FREQ_ORDER", "20"))

php = PhysicsHyperParams(
    λlb_μm = 8.0,
    λub_μm = 12.0,
    freq_order = freq_order,
    focal_length_μm = 2.0e4,
    num_unit_cells = 2048,
    unit_cell_length_μm = 4.0,
    pillar_width_lb_μm = 1.8,
    pillar_width_ub_μm = 2.7,
    pillar_width_order = 160,
    pillar_height_μm = 10.0,
    pillar_material = "Si_no_absorption",
    substrate_height_μm = 300.0,
    substrate_material = "Si_no_absorption",
    nG = 1000
)

@time End2EndThermalImg.compute_and_save_surrogate_transmission_matrix(php)
End2EndThermalImg.plot_surrogate_models(php)
