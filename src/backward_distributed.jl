# Custom adjoints crossing the pmap/remotecall process boundary that Zygote cannot
# differentiate through automatically. Two distinct pieces:
#
# 1. make_image_from_3D_distributed's rrule: bridges the pmap boundary around
#    make_image_from_3D_shard (itself plain, ordinarily-differentiable code -- no custom rule
#    needed for it directly). Only provides the `object` gradient, needed every inner
#    reconstruction-loop call.
#
# 2. get_fftPSFs_from_precomputed_distributed's rrule (+ Δgeoms_shard/Δgeoms_distributed): the
#    distributed analogue of backward.jl's get_fftPSFs_from_precomputed rrule, which threads a
#    local Zygote.gradient per (freq, z) pair across @threads on one process. Here those same
#    (freq, z) pairs are additionally sharded across workers (reusing the resident-shard state
#    from forward_distributed.jl: assign_shard_indices!/init_resident_incidents! ->
#    _shard_grid_inds/_incidents_shard), each worker still threading its own shard locally, and
#    the driver reduces (sums) the per-worker partial Δgeoms via pmap. This provides the outer
#    `geoms` gradient, via Λ_dot_dFdo_distributed below.
function ChainRulesCore.rrule(::typeof(make_image_from_3D_distributed), object, weights, php, imghp)
    image = make_image_from_3D_distributed(object, weights, php, imghp)
    function make_image_from_3D_distributed_pullback(Δimage)
        Δ = ChainRulesCore.unthunk(Δimage)
        Δobject_parts = pmap(workers()) do p
            remotecall_fetch(p, object, Δ, weights, php, imghp) do o, Δlocal, w, ph, ih
                Zygote.gradient(oo -> dot(make_image_from_3D_shard(oo, w, ph, ih), Δlocal), o)[1]
            end
        end
        # object is always a (; Tmap, zmap) NamedTuple (see unflatten_object, prepare.jl), and
        # Zygote's gradient w.r.t. it is a plain NamedTuple with the same fields -- which has no
        # `+` method for `sum` to use directly, so combine field-by-field instead.
        Δobject = reduce((a, b) -> (Tmap = a.Tmap .+ b.Tmap, zmap = a.zmap .+ b.zmap), Δobject_parts)
        return NoTangent(), Δobject, NoTangent(), NoTangent(), NoTangent()
    end
    return image, make_image_from_3D_distributed_pullback
end

# runs on a worker: this worker's shard of the (freq, z) Zygote.gradient accumulation, mirroring
# get_fftPSFs_from_precomputed's rrule pullback in backward.jl (same per-pair local
# Zygote.gradient + per-thread accumulation), but restricted to _shard_grid_inds/_incidents_shard
# (set up by assign_shard_indices!/init_resident_incidents!, forward_distributed.jl) instead of
# the full (freq, z) grid. surrogates/n2f_kernels/plans_n2f/plans_PSF are reached via `Main.`,
# same reason as recompute_fftPSFs_shard! above (this file's own module-vs-Main note): they're
# resident globals set up by the calling script's @everywhere setup, in Main, not this module.
function Δgeoms_at_pair(iF, iZ, geoms, ΔfftPSFs, freqs, php, imghp, tid, PSF_scale)
    # surrogates/n2f_kernels/plans_n2f/plans_PSF are pulled into local variables *before*
    # entering the closure below, rather than referenced as `Main.x` directly inside it.
    # `Main.x` is a cross-module property access, which lowers to an actual
    # `getproperty(::Module, ::Symbol)` call; since it's inside the closure Zygote
    # differentiates, Zygote has to build a differentiation rule for that getproperty call
    # itself (even though the result doesn't depend on `geoms`) -- and that specific, rarely
    # exercised path was observed to corrupt Zygote's generated pullback code (`UndefVarError:
    # j not defined in Zygote`, an internal Zygote symbol, not anything of ours). A bare global
    # referenced from *within* Main (as in the original notebook prototype, where this function
    # was `@everywhere`'d into Main directly) lowers to a more primitive GlobalRef and never
    # hits this path -- porting into this module (which needs the `Main.` prefix to reach the
    # same globals) introduced it. Resolving these outside the closure sidesteps it: the
    # closure then only captures already-resolved local values.
    surrogate = Main.surrogates[iF]
    n2f_k = Main.n2f_kernels[iF]
    plan_n2f = Main.plans_n2f[tid]
    plan_PSF = Main.plans_PSF[tid]
    # incident/Δy hoisted out of the closure too, for the same reason as surrogate/n2f_k/plan_n2f/
    # plan_PSF above: leaving a Dict/array getindex inside the closure Zygote differentiates gives
    # it one more (unnecessary, geoms-independent) expression to build a pullback for, and this
    # closure has already shown that expanding what it captures directly can corrupt Zygote's
    # generated code (see the comment above) -- confirmed empirically to matter here: with these
    # two still indexed inside the closure, this function returned exactly 2x the correct
    # sensitivity (found via an isolated rrule-only finite-difference check).
    incident = _incidents_shard[(iF, iZ)]
    Δy = ΔfftPSFs[iF, iZ]
    Zygote.gradient(
        g -> real(dot(get_fftPSF(get_PSF_at_freq_and_z(
            freqs[iF], incident, surrogate, g,
            n2f_k, plan_n2f, php, imghp, PSF_scale
        ), plan_PSF), Δy)),
        geoms
    )[1]
end

function Δgeoms_shard(geoms, ΔfftPSFs, freqs, php, imghp, PSF_scale)
    Δgeoms_per_thread = [zeros(size(geoms)) for _ in 1:nthreads()]
    @threads for i in eachindex(_shard_grid_inds)
        idx = _shard_grid_inds[i]
        iF, iZ = Tuple(idx)
        tid = threadid()
        Δgeoms_per_thread[tid] .+= Δgeoms_at_pair(iF, iZ, geoms, ΔfftPSFs, freqs, php, imghp, tid, PSF_scale)
    end
    sum(Δgeoms_per_thread)
end

# driver-only: dispatches Δgeoms_shard to every worker and reduces (sums) the partial results.
# unthunk mirrors make_image_from_3D_distributed_pullback above -- cotangents can arrive as lazy
# Thunks, which either don't serialize across the pmap boundary or would force redundant work on
# every worker; materialize once here instead (a no-op if already concrete).
#
# Dispatch is a sequential loop of non-blocking `remotecall`s, not `pmap` -- pmap fires all 5
# sends concurrently (via @async-scheduled tasks), which hit the same ClusterSerializer message
# corruption that distributed.jl's broadcast_value works around for the same reason (see its
# comment), once ΔfftPSFs is a real, sizeable, nested payload rather than a small test one.
# `remotecall` (unlike `remotecall_fetch`/`remotecall_wait`) only blocks on dispatching the
# message, not on the remote computation finishing, so sending sequentially here doesn't
# serialize the ~20s-per-worker compute -- all 5 workers still run concurrently once dispatched;
# only the sends themselves are staggered.
function Δgeoms_distributed(geoms, ΔfftPSFs, freqs, php, imghp, PSF_scale)
    ΔfftPSFs_concrete = ChainRulesCore.unthunk(ΔfftPSFs)
    futures = [remotecall(Δgeoms_shard, p, geoms, ΔfftPSFs_concrete, freqs, php, imghp, PSF_scale) for p in workers()]
    parts = fetch.(futures)
    sum(parts)
end

# distributed analogue of get_fftPSFs_from_precomputed (backward.jl): same no-op forward pass
# (geoms only matters to the gradient), but the rrule crosses the pmap boundary via
# Δgeoms_distributed instead of only @threads.
get_fftPSFs_from_precomputed_distributed(geoms, fftPSFs_precomputed, freqs, php, imghp, PSF_scale) = fftPSFs_precomputed

function ChainRulesCore.rrule(::typeof(get_fftPSFs_from_precomputed_distributed), geoms, fftPSFs_precomputed, freqs, php, imghp, PSF_scale)
    function pullback(ΔfftPSFs)
        Δgeoms = Δgeoms_distributed(geoms, ΔfftPSFs, freqs, php, imghp, PSF_scale)
        return NoTangent(), Δgeoms, NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent()
    end
    return fftPSFs_precomputed, pullback
end

# distributed analogue of Λ_dot_dFdo (backward.jl): identical except it calls
# get_fftPSFs_from_precomputed_distributed instead of get_fftPSFs_from_precomputed --
# incidents/surrogates/n2f_kernels/plans_n2f are dropped from the signature entirely, since
# they're no longer needed by the driver at all (only inside Δgeoms_shard, which pulls them
# from each worker's own resident state).
function Λ_dot_dFdo_distributed(Λ, object, object_opt, geoms, fftPSFs_precomputed, freqs, weights, plans_PSF, noise, noise_type::AbstractNoiseType, php, imghp, PSF_scale)
    PSF_zcoords = ignore_derivatives(() -> get_PSF_zcoords(imghp))
    δ_Δz = ignore_derivatives(() -> get_discretized_δ_function(imghp.smoothness_order, imghp.PSF_Δz))
    C_interp_3D = ignore_derivatives(() -> get_C_interp_3D(object_opt.zmap, PSF_zcoords, δ_Δz, imghp))
    dBdT_opt = ignore_derivatives(() -> [dB_dT.(object_opt.Tmap, freq, php.wavcen) for freq in freqs])
    B_opt = ignore_derivatives(() -> get_black_body_spectrum(object_opt.Tmap, php))
    dC_dz_opt = ignore_derivatives(() -> [dC_dz.(object_opt.zmap, PSF_zcoord, imghp.smoothness_order, imghp.PSF_Δz) for PSF_zcoord in PSF_zcoords])

    fftPSFs = get_fftPSFs_from_precomputed_distributed(geoms, fftPSFs_precomputed, freqs, php, imghp, PSF_scale)
    plan_PSF = plans_PSF[1]

    v = make_noisy_image_from_3D(object, fftPSFs, weights, noise, noise_type, plan_PSF, php, imghp)
    image_opt = make_image_from_3D(object_opt, fftPSFs, weights, plan_PSF, php, imghp)
    image_diff = v .- image_opt

    Λ_T_2D = ignore_derivatives(() -> reshape(Λ[1:imghp.objN^2], imghp.objN, imghp.objN))

    # see Λ_dot_dFdo (backward.jl) for why this scatter needs ignore_derivatives specifically
    Λ_z_2D = ignore_derivatives() do
        Λ_z = Λ[imghp.objN^2 + 1:end]
        _, active_linear_indices = get_active_variables(object_opt.zmap, imghp)
        free_linear_indices = setdiff(1:imghp.objN^2, active_linear_indices)
        Λ_z_full = zeros(eltype(Λ_z), imghp.objN^2)
        Λ_z_full[free_linear_indices] .= Λ_z
        reshape(Λ_z_full, imghp.objN, imghp.objN)
    end

    term1 = sum(
        iF -> sum(
            iZ -> begin
                spectrum = Λ_T_2D .* C_interp_3D[:, :, iZ] .* dBdT_opt[iF]
                weights[iF] .* convolve_with_fftPSF_at_freq_and_z(fftPSFs[iF, iZ], spectrum, plan_PSF)
            end,
            eachindex(PSF_zcoords)),
        eachindex(freqs))
    term1 = 2 * dot(image_diff, term1)

    term2 = sum(
        iF -> sum(
            iZ -> begin
                spectrum = Λ_z_2D .* B_opt[iF] .* dC_dz_opt[iZ]
                weights[iF] .* convolve_with_fftPSF_at_freq_and_z(fftPSFs[iF, iZ], spectrum, plan_PSF)
            end,
            eachindex(PSF_zcoords)),
        eachindex(freqs))
    term2 = 2 * dot(image_diff, term2)

    # sign fix: see Λ_dot_dFdo (backward.jl) -- this returns -Λ⋅∇_o L(object_opt, geoms) as built
    # above, confirmed numerically (dot(Λ, ∇_o L) computed directly via Zygote came out to
    # -1.0002x the un-negated term1+term2).
    return -(term1 + term2)
end

# Full gradient of the outer objective w.r.t. [geoms; α; β], for the eventual outer optimization
# that will jointly optimize geoms alongside the regularization weights (see rechp's α_init/β_init,
# src/prepare.jl). The geoms block is just Λ_dot_dFdo_distributed, autodiffed through its geoms
# argument as before. The α/β blocks don't need autodiff: since the inner objective's only α/β-
# dependence is α*sum((Tmap - T_background)^2) + β*sum((zmap - z_middle)^2), its mixed partials
# w.r.t. (object_opt, α) and (object_opt, β) are just 2*(T_opt - T_background) and
# 2*(z_opt - z_middle) respectively (zero cross terms between the T/z blocks), evaluated at
# object_opt (the same point Λ itself is defined at). Contracting those against the matching block
# of Λ gives dL/dα, dL/dβ directly via the same adjoint identity used for dL/dgeoms -- no second
# Hessian solve needed. z_opt's active (bound-pinned) pixels are dropped the same way Λ_z/b_z
# already are (see get_active_variables / Λ_dot_dFdo_distributed above).
function outer_gradient_distributed(Λ, object, object_opt, geoms, fftPSFs_precomputed, freqs, weights, plans_PSF, noise, noise_type::AbstractNoiseType, rechp, php, imghp, PSF_scale)
    ∇geoms = Zygote.gradient(
        g -> Λ_dot_dFdo_distributed(Λ, object, object_opt, g, fftPSFs_precomputed, freqs, weights, plans_PSF, noise, noise_type, php, imghp, PSF_scale),
        geoms
    )[1]

    @unpack T_background, z_middle = rechp
    _, active_linear_indices = get_active_variables(object_opt.zmap, imghp)
    free_linear_indices = setdiff(1:imghp.objN^2, active_linear_indices)

    Λ_T = Λ[1:imghp.objN^2]
    Λ_z = Λ[imghp.objN^2 + 1:end]

    dL_dα = 2 * dot(Λ_T, vec(object_opt.Tmap .- T_background))
    dL_dβ = 2 * dot(Λ_z, vec(object_opt.zmap .- z_middle)[free_linear_indices])

    return [vec(∇geoms); dL_dα; dL_dβ]
end

# Distributed analogue of build_hessian_operator (backward.jl). Like Λ_dot_dFdo_distributed,
# only the per-call Hessian-vector product itself is distributed -- image_diff still comes from
# a concrete, gathered fftPSFs matrix on the driver, same reason as Λ_dot_dFdo_distributed's own
# forward pass. Unlike make_image_from_3D_shard (one map + one reduce), hessian_vector_product!
# has a sandwich structure: intermediate_TT/intermediate_Tz are *global* sums over every (iF, iZ)
# pair, and term_TT_1/term_zT_1/term_Tz_1/term_zz_1 each need a transpose-convolution of that
# same global intermediate at every (iF, iZ) pair again. So each call needs two remote
# round-trips instead of one:
#   1. hvp_forward_shard: each worker forward-convolves v_T_2D/v_z_2D against its own
#      _fftPSFs_shard; pmap+remotecall_fetch sums the partials into the full
#      intermediate_TT/intermediate_Tz on the driver.
#   2. hvp_backward_shard: the (tiny, objN x objN) global intermediates go back out to every
#      worker, which transpose-convolve them against its own shard; gather+sum again.
# sum_TT_2/sum_Tz_2/sum_zz_2 (the image_diff-based residual-curvature terms) are v-independent,
# so they're computed once via the same sharded map-reduce at setup time (init_hvp_shard!/
# init_hvp_resident_distributed!), not per call -- same idea as gathering fftPSFs once instead of
# per (freq, z) query. dBdT_opt/C_interp_3D/B_opt/dC_dz_opt/d2B_dT2_opt/d2C_dz2_opt only exist on
# the driver otherwise; since they're small (not shard-sized like fftPSFs), each worker just
# computes its own full copies locally from object_opt (passed as a plain argument, not
# broadcast) rather than sharding them too.

# runs on a worker: computes this worker's full copies of the v-independent quantities needed by
# hessian_vector_product_distributed! below, plus this worker's shard of the image_diff-based
# residual-curvature sums.
function init_hvp_shard!(object_opt, image_diff, freqs, PSF_zcoords, weights, php, imghp)
    δ_Δz = get_discretized_δ_function(imghp.smoothness_order, imghp.PSF_Δz)
    global _hvp_C_interp_3D = get_C_interp_3D(object_opt.zmap, PSF_zcoords, δ_Δz, imghp)
    global _hvp_dBdT_opt = [dB_dT.(object_opt.Tmap, freq, php.wavcen) for freq in freqs]
    global _hvp_d2B_dT2_opt = [d2B_dT2.(object_opt.Tmap, freq, php.wavcen) for freq in freqs]
    global _hvp_B_opt = get_black_body_spectrum(object_opt.Tmap, php)
    global _hvp_dC_dz_opt = [dC_dz.(object_opt.zmap, PSF_zcoord, imghp.smoothness_order, imghp.PSF_Δz) for PSF_zcoord in PSF_zcoords]
    global _hvp_d2C_dz2_opt = [d2C_dz2.(object_opt.zmap, PSF_zcoord, imghp.smoothness_order, imghp.PSF_Δz) for PSF_zcoord in PSF_zcoords]

    plan_PSF = Main.plans_PSF[1]
    sum_TT_2_partial = zeros(imghp.objN, imghp.objN)
    sum_Tz_2_partial = zeros(imghp.objN, imghp.objN)
    sum_zz_2_partial = zeros(imghp.objN, imghp.objN)
    for idx in _shard_grid_inds
        iF, iZ = Tuple(idx)
        conv = weights[iF] .* convolveT_with_fftPSF_at_freq_and_z(_fftPSFs_shard[(iF, iZ)], image_diff, plan_PSF)
        sum_TT_2_partial .+= conv .* _hvp_C_interp_3D[:, :, iZ] .* _hvp_d2B_dT2_opt[iF]
        sum_Tz_2_partial .+= conv .* _hvp_dC_dz_opt[iZ] .* _hvp_dBdT_opt[iF]
        sum_zz_2_partial .+= conv .* _hvp_B_opt[iF] .* _hvp_d2C_dz2_opt[iZ]
    end

    (sum_TT_2_partial, sum_Tz_2_partial, sum_zz_2_partial)
end

# driver-only: dispatches init_hvp_shard! to every worker and reduces (sums) the three partial
# v-independent residual-curvature sums, once, for the whole minres! solve.
function init_hvp_resident_distributed!(object_opt, image_diff, freqs, PSF_zcoords, weights, php, imghp)
    parts = pmap(p -> remotecall_fetch(init_hvp_shard!, p, object_opt, image_diff, freqs, PSF_zcoords, weights, php, imghp), workers())
    sum_TT_2 = sum(p[1] for p in parts)
    sum_Tz_2 = sum(p[2] for p in parts)
    sum_zz_2 = sum(p[3] for p in parts)
    (sum_TT_2, sum_Tz_2, sum_zz_2)
end

# round 1 (see the header comment above): intermediate_TT/intermediate_Tz don't share any
# computation with each other, so there's no redundancy to merge here -- two independent
# sum(...) calls, which correctly infer their output shape (imgN x imgN, NOT objN x objN --
# convolve_with_fftPSF_at_freq_and_z maps an objN x objN spectrum through the larger fftPSF
# kernel up to the image-domain size) from the actual values.
function hvp_forward_shard(v_T_2D, v_z_2D, weights)
    plan_PSF = Main.plans_PSF[1]
    intermediate_TT_partial = sum(_shard_grid_inds) do idx
        iF, iZ = Tuple(idx)
        spectrum_T = v_T_2D .* _hvp_dBdT_opt[iF] .* _hvp_C_interp_3D[:, :, iZ]
        weights[iF] .* convolve_with_fftPSF_at_freq_and_z(_fftPSFs_shard[(iF, iZ)], spectrum_T, plan_PSF)
    end
    intermediate_Tz_partial = sum(_shard_grid_inds) do idx
        iF, iZ = Tuple(idx)
        spectrum_z = v_z_2D .* _hvp_B_opt[iF] .* _hvp_dC_dz_opt[iZ]
        weights[iF] .* convolve_with_fftPSF_at_freq_and_z(_fftPSFs_shard[(iF, iZ)], spectrum_z, plan_PSF)
    end
    (intermediate_TT_partial, intermediate_Tz_partial)
end

# round 2 (see the header comment above): term_TT_1/term_zT_1 (and term_Tz_1/term_zz_1) DO share
# a convolution, so we merge those loops -- but convolveT_with_fftPSF_at_freq_and_z maps the
# imgN x imgN intermediate back down to objN x objN, so a zeros(size(intermediate_TT))-style
# pre-allocation would be the wrong shape. Peel the first shard index to seed each accumulator
# from an actual computed value instead of a guess.
function hvp_backward_shard(intermediate_TT, intermediate_Tz, weights)
    plan_PSF = Main.plans_PSF[1]
    grid_inds = collect(_shard_grid_inds)

    iF, iZ = Tuple(grid_inds[1])
    conv_TT = weights[iF] .* convolveT_with_fftPSF_at_freq_and_z(_fftPSFs_shard[(iF, iZ)], intermediate_TT, plan_PSF)
    conv_Tz = weights[iF] .* convolveT_with_fftPSF_at_freq_and_z(_fftPSFs_shard[(iF, iZ)], intermediate_Tz, plan_PSF)
    term_TT_1_partial = (_hvp_dBdT_opt[iF] .* _hvp_C_interp_3D[:, :, iZ]) .* conv_TT
    term_zT_1_partial = (_hvp_B_opt[iF] .* _hvp_dC_dz_opt[iZ]) .* conv_TT
    term_Tz_1_partial = (_hvp_dBdT_opt[iF] .* _hvp_C_interp_3D[:, :, iZ]) .* conv_Tz
    term_zz_1_partial = (_hvp_dC_dz_opt[iZ] .* _hvp_B_opt[iF]) .* conv_Tz

    for idx in grid_inds[2:end]
        iF, iZ = Tuple(idx)
        conv_TT = weights[iF] .* convolveT_with_fftPSF_at_freq_and_z(_fftPSFs_shard[(iF, iZ)], intermediate_TT, plan_PSF)
        conv_Tz = weights[iF] .* convolveT_with_fftPSF_at_freq_and_z(_fftPSFs_shard[(iF, iZ)], intermediate_Tz, plan_PSF)
        term_TT_1_partial .+= (_hvp_dBdT_opt[iF] .* _hvp_C_interp_3D[:, :, iZ]) .* conv_TT
        term_zT_1_partial .+= (_hvp_B_opt[iF] .* _hvp_dC_dz_opt[iZ]) .* conv_TT
        term_Tz_1_partial .+= (_hvp_dBdT_opt[iF] .* _hvp_C_interp_3D[:, :, iZ]) .* conv_Tz
        term_zz_1_partial .+= (_hvp_dC_dz_opt[iZ] .* _hvp_B_opt[iF]) .* conv_Tz
    end

    (term_TT_1_partial, term_zT_1_partial, term_Tz_1_partial, term_zz_1_partial)
end

# α/β: as in build_hessian_operator (backward.jl) -- required, and must equal the values given to
# the reconstruction, since the operator is H_data + diag(2α on T, 2β on the free z pixels).
function build_hessian_operator_distributed(object, object_opt, fftPSFs, freqs, weights, plans_PSF, noise, noise_type::AbstractNoiseType, α, β, php, imghp)
    PSF_zcoords = get_PSF_zcoords(imghp)
    _, active_linear_indices = get_active_variables(object_opt.zmap, imghp)
    free_linear_indices = setdiff(1:imghp.objN^2, active_linear_indices)
    N = 2 * imghp.objN^2 - length(active_linear_indices)

    plan_PSF = plans_PSF[1]
    noisy_image = make_noisy_image_from_3D(object, fftPSFs, weights, noise, noise_type, plan_PSF, php, imghp)
    image_opt = make_image_from_3D(object_opt, fftPSFs, weights, plan_PSF, php, imghp)
    image_diff = noisy_image .- image_opt

    sum_TT_2, sum_Tz_2, sum_zz_2 = init_hvp_resident_distributed!(object_opt, image_diff, freqs, PSF_zcoords, weights, php, imghp)

    function hessian_vector_product_distributed!(Hv::AbstractVector, v::AbstractVector)
        v_T = @view v[1:imghp.objN^2]
        v_z = @view v[imghp.objN^2+1:end]

        v_T_2D = reshape(v_T, imghp.objN, imghp.objN)

        v_z_full = zeros(eltype(v_z), imghp.objN^2)
        v_z_full[free_linear_indices] .= v_z
        v_z_2D = reshape(v_z_full, imghp.objN, imghp.objN)

        # round 1: forward-convolve v_T_2D/v_z_2D against each worker's shard, reduce to the
        # full (global) intermediate arrays
        parts1 = pmap(p -> remotecall_fetch(hvp_forward_shard, p, v_T_2D, v_z_2D, weights), workers())
        intermediate_TT = sum(p[1] for p in parts1)
        intermediate_Tz = sum(p[2] for p in parts1)

        # round 2: transpose-convolve the now-global intermediates against each worker's shard
        # again, reduce the four partial term arrays
        parts2 = pmap(p -> remotecall_fetch(hvp_backward_shard, p, intermediate_TT, intermediate_Tz, weights), workers())
        term_TT_1 = sum(p[1] for p in parts2)
        term_zT_1 = sum(p[2] for p in parts2)
        term_Tz_1 = sum(p[3] for p in parts2)
        term_zz_1 = sum(p[4] for p in parts2)

        term_TT_2 = sum_TT_2 .* v_T_2D
        term_Tz_2 = sum_Tz_2 .* v_z_2D
        term_zT_2 = sum_Tz_2 .* v_T_2D
        term_zz_2 = sum_zz_2 .* v_z_2D

        # 2α .* v_T / 2β .* v_z_full are the regularization Hessians of the T-T and z-z blocks
        # (v_z_full is zero at the active pixels, and Hv_zz is on the full grid before compaction)
        Hv_TT = vec(2 .* term_TT_1) .+ vec(-2 .* term_TT_2) .+ 2α .* v_T
        Hv_Tz = vec(2 .* term_Tz_1) .+ vec(-2 .* term_Tz_2)
        Hv[1:imghp.objN^2] .= Hv_TT .+ Hv_Tz

        Hv_zT = vec(2 .* term_zT_1) .+ vec(-2 .* term_zT_2)
        Hv_zz = vec(2 .* term_zz_1) .+ vec(-2 .* term_zz_2) .+ 2β .* v_z_full
        Hv_z_full = Hv_zT .+ Hv_zz
        Hv[imghp.objN^2+1:end] .= Hv_z_full[free_linear_indices]

        Hv
    end

    LinearMaps.LinearMap{Float64}(hessian_vector_product_distributed!, N; ismutating=true, issymmetric=true)
end
