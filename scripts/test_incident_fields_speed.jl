using End2EndThermalImg
using Distributed

nw = End2EndThermalImg.setup_cluster_workers()
println("nworkers = ", nworkers(), "  workers = ", workers())
println("driver threads = ", Threads.nthreads())
println("threads per worker = ", [remotecall_fetch(Threads.nthreads, p) for p in workers()])
println("driver hostname = ", gethostname())
println("hostname per worker = ", [remotecall_fetch(gethostname, p) for p in workers()])

@everywhere begin
    php = PhysicsHyperParams(
        λlb_μm = 8.0,
        λub_μm = 12.0,
        freq_order = 20,
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

    imghp = ImagingHyperParams(
        objN = 32,
        imgN = 128,
        binN = 3,
        sampleN = 1,
        PSF_zlb_μm = -2.1e6,
        PSF_zub_μm = -0.9e6,
        PSF_zlen = 11,
        PSF_scale = 10.0^12,
        smoothness_order = 3.0,
        object_type = GaussianObject(Tlb = 136.15, Tub = 450.15, std_dev_T = 6.0, zlb_μm = -2.0e6, zub_μm = -1.0e6, std_dev_z = 6.0, php = php),
        noise_level = 0.02,
        php = php
    )
end

freqs = End2EndThermalImg.get_freq_chebpoints(php)
PSF_zcoords = End2EndThermalImg.get_PSF_zcoords(imghp)

@everywhere using Base.Threads

@everywhere function compute_incident_shard(shard_idx, nshards, freqs, PSF_zcoords, php)
    grid_inds = End2EndThermalImg.shard_grid_indices((length(freqs), length(PSF_zcoords)), nshards, shard_idx)
    shard_keys = Tuple.(grid_inds)
    shard_vals = Vector{Matrix{ComplexF64}}(undef, length(grid_inds))
    @threads for i in eachindex(grid_inds)
        idx = grid_inds[i]
        shard_vals[i] = End2EndThermalImg.get_incident_field(freqs[idx[1]], PSF_zcoords[idx[2]], php)
    end
    Dict(zip(shard_keys, shard_vals))
end

# distributed: shard the (freq, z) grid across every worker/node
nshards = nworkers()
println("\n--- distributed: pmap across $(nshards) worker(s) ---")
@time incidents_parts = pmap(i -> compute_incident_shard(i, nshards, freqs, PSF_zcoords, php), 1:nshards)

incidents_distributed = Matrix{Matrix{ComplexF64}}(undef, length(freqs), length(PSF_zcoords))
for part in incidents_parts, (key, val) in part
    incidents_distributed[key...] = val
end
incidents_parts = nothing
GC.gc()

# single-process: @threads on the driver only
println("\n--- single-process: @threads on driver ($(Threads.nthreads()) threads) ---")
@time incidents_single_process = End2EndThermalImg.get_incident_fields(freqs, PSF_zcoords, php)

println("\nResults match: ", incidents_distributed == incidents_single_process)
