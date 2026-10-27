function dB_dT(T, freq, wavcen)
    A = (h * c * 10^6)/(wavcen * kB)
    numerator = (2 * A * freq^4 * exp(A * freq / T) ) 
    denominator = T^2 * (exp(A * freq / T) - 1)^2
    numerator/denominator
end

function d2B_dT2(T, freq, wavcen)
    A = (h * c * 10^6)/(wavcen * kB)
    numerator = A * freq^4 * (-2*T + A * freq * coth(A*freq/(2*T) ) )
    denominator = T^4 * (-1 + cosh(A*freq/T) )
    numerator/denominator
end  

function df_δ_dz(z, smoothness_order)
    if smoothness_order == Inf
        return z > 0.0 ? exp(-1 / z) / z^2 : 0.0
    else
        return z > 0.0 ? (smoothness_order + 1)*z^smoothness_order : 0.0
    end
end

function d2f_δ_dz2(z, smoothness_order)
    if smoothness_order == Inf
        return z > 0.0 ? exp(-1 / z) * (1 - 2z) / z^4 : 0.0
    else
        return z > 0.0 ? (smoothness_order + 1)*smoothness_order*z^(smoothness_order - 1) : 0.0
    end
end

function dg_δ_dz(z, smoothness_order)
    g_denom = (f_δ(z, smoothness_order) + f_δ(1-z, smoothness_order))
    numerator_term1 = g_denom*df_δ_dz(z, smoothness_order)
    numerator_term2 = f_δ(z, smoothness_order)*(df_δ_dz(z, smoothness_order) - df_δ_dz(1-z, smoothness_order))
    denominator = g_denom ^2
    (numerator_term1 - numerator_term2)/denominator
end

function d2g_δ_dz2(z, smoothness_order)
    D = f_δ(z, smoothness_order) + f_δ(1-z, smoothness_order)
    dD_dz = df_δ_dz(z, smoothness_order) - df_δ_dz(1-z, smoothness_order)
    d2D_dz2 = d2f_δ_dz2(z, smoothness_order) + d2f_δ_dz2(1-z, smoothness_order)
    N = df_δ_dz(z, smoothness_order)*D - f_δ(z, smoothness_order)*dD_dz
    dN_dz = d2f_δ_dz2(z, smoothness_order)*D - f_δ(z, smoothness_order)*d2D_dz2
    (dN_dz*D - 2*N*dD_dz)/D^3
end

function dδ_dz(z, smoothness_order, Δz)
    (1/Δz)*((-1/Δz)*dg_δ_dz(z/Δz, smoothness_order) + (1/Δz)*dg_δ_dz(-z/Δz, smoothness_order))
end

function d2δ_dz2(z, smoothness_order, Δz)
    (1/Δz)*((-1/Δz^2)*d2g_δ_dz2(z/Δz, smoothness_order) + (-1/Δz^2)*d2g_δ_dz2(-z/Δz, smoothness_order))
end

function dC_dz(z, PSF_zcoord, smoothness_order, Δz)
    diff = PSF_zcoord - z
    -dδ_dz(diff, smoothness_order, Δz)
end

function d2C_dz2(z, PSF_zcoord, smoothness_order, Δz)
    diff = PSF_zcoord - z
    d2δ_dz2(diff, smoothness_order, Δz)
end

function get_fftPSFs_from_precomputed(geoms, fftPSFs_precomputed, freqs, incidents, surrogates, n2f_kernels, plans_n2f, plans_PSF, php, imghp, PSF_scale)
    fftPSFs_precomputed
end

function ChainRulesCore.rrule(::typeof(get_fftPSFs_from_precomputed), geoms, fftPSFs_precomputed, freqs, incidents, surrogates, n2f_kernels, plans_n2f, plans_PSF, php, imghp, PSF_scale)
    PSF_zlen = imghp.PSF_zlen
    nF = length(freqs)
    inds = CartesianIndices((1:nF, 1:PSF_zlen))

    function pullback(ΔfftPSFs)
        Δgeoms_per_thread = [zeros(size(geoms)) for _ in 1:nthreads()]
        @threads for idx in eachindex(inds)
            tid = threadid()
            I = inds[idx]
            iF, iZ = I[1], I[2]
            # everything the closure needs hoisted into locals first (matching Δgeoms_at_pair,
            # backward_distributed.jl): leaving array getindex/cross-module calls inside the
            # closure Zygote differentiates gives it more (geoms-independent) expressions to
            # build a pullback for, and this exact pattern has been confirmed (in the distributed
            # twin) to silently corrupt the computed gradient by exactly a factor of 2.
            incident = incidents[iF, iZ]
            surrogate = surrogates[iF]
            n2f_k = n2f_kernels[iF]
            plan_n2f = plans_n2f[tid]
            plan_PSF = plans_PSF[tid]
            Δy = ΔfftPSFs[iF, iZ]
            Δgeoms_local = Zygote.gradient(
                g -> real(dot(get_fftPSF(get_PSF_at_freq_and_z(
                    freqs[iF], incident, surrogate, g,
                    n2f_k, plan_n2f, php, imghp, PSF_scale
                ), plan_PSF), Δy)),
                geoms
            )[1]
            Δgeoms_per_thread[tid] .+= Δgeoms_local
        end
        Δgeoms = sum(Δgeoms_per_thread)
        return NoTangent(), Δgeoms, NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent()
    end

    return fftPSFs_precomputed, pullback
end

# Computes Λ ⋅ dF/do, the adjoint-weighted sensitivity of the forward image-formation map F
# (parameterized by the metasurface geometry o = geoms) needed for the outer geoms-optimization
# gradient. Λ is the concatenation of the reconstruction adjoint state for T and z:
# Λ = [Λ_T; Λ_z], where Λ_T has length objN^2 - (# active T pixels) and Λ_z has length
# objN^2 - (# active z pixels): pixels pinned at their bounds (T_LOWER_BOUND for T,
# PSF_zlb/PSF_zub for z) are dropped from the reconstruction's optimization variables entirely, so
# they're simply absent from Λ_T/Λ_z rather than zero.
function Λ_dot_dFdo(Λ, object, object_opt, geoms, fftPSFs_precomputed, freqs, incidents, surrogates, n2f_kernels, plans_n2f, plans_PSF, weights, noise, noise_type::AbstractNoiseType, php, imghp, PSF_scale)
    PSF_zcoords = ignore_derivatives(() -> get_PSF_zcoords(imghp))
    δ_Δz = ignore_derivatives(() -> get_discretized_δ_function(imghp.smoothness_order, imghp.PSF_Δz))
    C_interp_3D = ignore_derivatives(() -> get_C_interp_3D(object_opt.zmap, PSF_zcoords, δ_Δz, imghp))
    dBdT_opt = ignore_derivatives(() -> [dB_dT.(object_opt.Tmap, freq, php.wavcen) for freq in freqs])
    B_opt = ignore_derivatives(() -> get_black_body_spectrum(object_opt.Tmap, php))
    dC_dz_opt = ignore_derivatives(() -> [dC_dz.(object_opt.zmap, PSF_zcoord, imghp.smoothness_order, imghp.PSF_Δz) for PSF_zcoord in PSF_zcoords])

    # fftPSFs is generated from geoms via get_fftPSFs_from_precomputed, whose custom rrule
    # above returns fftPSFs_precomputed as-is on the forward pass (cheap) but computes the
    # correct d(fftPSFs)/d(geoms) on the backward pass by threading a local Zygote.gradient
    # per (freq, z) pair -- this is what makes Λ_dot_dFdo differentiable w.r.t. geoms via
    # ordinary Zygote.gradient, without hand-deriving that piece.
    fftPSFs = get_fftPSFs_from_precomputed(geoms, fftPSFs_precomputed, freqs, incidents, surrogates, n2f_kernels, plans_n2f, plans_PSF, php, imghp, PSF_scale)
    plan_PSF = plans_PSF[1]

    v = make_noisy_image_from_3D(object, fftPSFs, weights, noise, noise_type, plan_PSF, php, imghp)
    image_opt = make_image_from_3D(object_opt, fftPSFs, weights, plan_PSF, php, imghp)
    image_diff = v .- image_opt

    # scatter Λ_T/Λ_z back into their free (non-active) pixels, in ascending linear-index order,
    # leaving zeros at the active positions. ignore_derivatives isn't just for consistency with the
    # rest of this function -- the mutating `.=` scatter can't be differentiated through by Zygote
    # at all.
    free_T, free_z = ignore_derivatives() do
        idx = get_free_and_active_indices(object_opt, imghp)
        idx.free_T, idx.free_z
    end
    Λ_T_2D = ignore_derivatives() do
        Λ_T = Λ[1:length(free_T)]
        Λ_T_full = zeros(eltype(Λ_T), imghp.objN^2)
        Λ_T_full[free_T] .= Λ_T
        reshape(Λ_T_full, imghp.objN, imghp.objN)
    end
    Λ_z_2D = ignore_derivatives() do
        Λ_z = Λ[length(free_T) + 1:end]
        Λ_z_full = zeros(eltype(Λ_z), imghp.objN^2)
        Λ_z_full[free_z] .= Λ_z
        reshape(Λ_z_full, imghp.objN, imghp.objN)
    end

    # first term: 2 * dot(image_diff, sum over (freq, z) of weights[iF] .* convolve(fftPSFs[iF, iZ], spectrum)),
    # where spectrum = Λ_T_2D .* C_interp_3D[:,:,iZ] .* dBdT_opt[iF].
    term1 = sum(
        iF -> sum(
            iZ -> begin
                spectrum = Λ_T_2D .* C_interp_3D[:, :, iZ] .* dBdT_opt[iF]
                weights[iF] .* convolve_with_fftPSF_at_freq_and_z(fftPSFs[iF, iZ], spectrum, plan_PSF)
            end,
            eachindex(PSF_zcoords)),
        eachindex(freqs))
    term1 = 2 * dot(image_diff, term1)

    # second term: same structure as term1, but spectrum = Λ_z_2D .* B_opt[iF] .* dC_dz_opt[iZ].
    term2 = sum(
        iF -> sum(
            iZ -> begin
                spectrum = Λ_z_2D .* B_opt[iF] .* dC_dz_opt[iZ]
                weights[iF] .* convolve_with_fftPSF_at_freq_and_z(fftPSFs[iF, iZ], spectrum, plan_PSF)
            end,
            eachindex(PSF_zcoords)),
        eachindex(freqs))
    term2 = 2 * dot(image_diff, term2)

    # term1/term2 as built above equal -Λ⋅∇_o L(object_opt, geoms), not +Λ⋅∇_o L: L's data term is
    # (image - noisy_image)^2, so d(image)/do enters with an extra minus relative to image_diff =
    # noisy_image - image_opt used here, which term1/term2's convolution-adjoint construction
    # doesn't include. Confirmed numerically (not just by this derivation): dot(Λ, ∇_o L(object_opt))
    # computed directly via Zygote.gradient of reconstruction_objective_distributed came out to
    # -1.0002x this function's un-negated return value.
    return -(term1 + term2)
end

# lb/ub passed explicitly (rather than looked up from imghp) so this works for both zmap
# (imghp.PSF_zlb/PSF_zub) and Tmap (T_LOWER_BOUND, Inf -- T has no upper bound, and a finite Tmap
# pixel is never == Inf, so this is "no upper bound" with no special-casing). The == exact-equality
# check relies on the scale/unscale (to_scaled/to_physical) round-trip in
# forward[_distributed].jl staying exact for the bound values -- already relied on for z, and for
# T's own current pinning check (src/optimize.jl).
function get_active_variables(map, lb, ub)
    cartesian_indices = findall(x -> x == lb || x == ub, map)
    linear_indices = LinearIndices(map)[cartesian_indices]
    return cartesian_indices, linear_indices
end

# Shared free/active linear-index bookkeeping for both Tmap (pinned at T_LOWER_BOUND) and zmap
# (pinned at PSF_zlb/PSF_zub), needed by every adjoint/Hessian call site that touches Λ, b, or the
# outer gradient -- Λ/b/Hv are each the concatenation of a (now-compacted) T block and z block.
function get_free_and_active_indices(object_opt, imghp)
    _, active_T = get_active_variables(object_opt.Tmap, T_LOWER_BOUND, Inf)
    _, active_z = get_active_variables(object_opt.zmap, imghp.PSF_zlb, imghp.PSF_zub)
    free_T = setdiff(1:imghp.objN^2, active_T)
    free_z = setdiff(1:imghp.objN^2, active_z)
    (; free_T, free_z, active_T, active_z)
end

# Builds a matrix-free LinearMaps.LinearMap for the Hessian of the inner reconstruction
# objective w.r.t. the active object variables (Tmap/zmap pixels), evaluated at object_opt -- the
# same active-variable space that Λ lives in in Λ_dot_dFdo/Λ_dot_dFdo_distributed (Λ_T has length
# objN^2 - (# active T pixels), Λ_z has length objN^2 - (# active z pixels), see the comment above
# Λ_dot_dFdo). The inner objective is the data misfit plus α||T - T_bg||^2 + β||z - z_middle||^2,
# so the operator is H_data + diag(2α on the free T pixels, 2β on the free z pixels): the 2α/2β
# terms are the T-T and z-z blocks' regularization parts. α and β must therefore be the same
# values given to reconstruct_Tmap_and_zmap[_distributed]; when they are large relative to the
# (tiny) data-misfit curvature the regularization dominates the true Hessian. Intended for
# IterativeSolvers.minres!(x0, A, b) to solve the adjoint equation A * Λ = b for Λ without ever
# forming A explicitly. Not differentiated through (unlike Λ_dot_dFdo/get_fftPSFs_from_precomputed
# above), so fftPSFs is just a concrete matrix here -- no rrule/geoms machinery needed.
function build_hessian_operator(object, object_opt, fftPSFs, freqs, weights, plans_PSF, noise, noise_type::AbstractNoiseType, α, β, php, imghp)
    PSF_zcoords = get_PSF_zcoords(imghp)
    idx = get_free_and_active_indices(object_opt, imghp)
    free_T, free_z = idx.free_T, idx.free_z
    N = length(free_T) + length(free_z)

    # quantities that only depend on object_opt/freqs/php/imghp (not on v), so compute them once
    # here rather than re-deriving them on every hessian_vector_product! call / minres! iteration
    δ_Δz = get_discretized_δ_function(imghp.smoothness_order, imghp.PSF_Δz)
    C_interp_3D = get_C_interp_3D(object_opt.zmap, PSF_zcoords, δ_Δz, imghp)
    dBdT_opt = [dB_dT.(object_opt.Tmap, freq, php.wavcen) for freq in freqs]
    d2B_dT2_opt = [d2B_dT2.(object_opt.Tmap, freq, php.wavcen) for freq in freqs]
    B_opt = get_black_body_spectrum(object_opt.Tmap, php)
    dC_dz_opt = [dC_dz.(object_opt.zmap, PSF_zcoord, imghp.smoothness_order, imghp.PSF_Δz) for PSF_zcoord in PSF_zcoords]
    d2C_dz2_opt = [d2C_dz2.(object_opt.zmap, PSF_zcoord, imghp.smoothness_order, imghp.PSF_Δz) for PSF_zcoord in PSF_zcoords]
    plan_PSF = plans_PSF[1]

    noisy_image = make_noisy_image_from_3D(object, fftPSFs, weights, noise, noise_type, plan_PSF, php, imghp)
    image_opt = make_image_from_3D(object_opt, fftPSFs, weights, plan_PSF, php, imghp)
    image_diff = noisy_image .- image_opt

    # conv_image_diff[iF, iZ] doesn't depend on v either -- it's the transpose-convolution of the
    # fixed image_diff residual against each (freq, z) fftPSF. The four residual-curvature terms
    # below (H_TT/H_Tz/H_zT/H_zz's "_2" parts) all start from this identical convolution, so
    # precompute it once here instead of 4 separate times per (iF, iZ) on every single call.
    conv_image_diff = [weights[iF] .* convolveT_with_fftPSF_at_freq_and_z(fftPSFs[iF, iZ], image_diff, plan_PSF)
                        for iF in eachindex(freqs), iZ in eachindex(PSF_zcoords)]

    # the residual-curvature terms are fully v-independent apart from the final .* v_T_2D /
    # .* v_z_2D, so precompute the double sums themselves too. sum_Tz_2 covers both H_Tz's and
    # H_zT's residual-curvature term, since dC_dz_opt[iZ] .* dBdT_opt[iF] is the same factor in both.
    sum_TT_2 = sum(iF -> sum(iZ -> conv_image_diff[iF, iZ] .* C_interp_3D[:, :, iZ] .* d2B_dT2_opt[iF], eachindex(PSF_zcoords)), eachindex(freqs))
    sum_Tz_2 = sum(iF -> sum(iZ -> conv_image_diff[iF, iZ] .* dC_dz_opt[iZ] .* dBdT_opt[iF], eachindex(PSF_zcoords)), eachindex(freqs))
    sum_zz_2 = sum(iF -> sum(iZ -> conv_image_diff[iF, iZ] .* B_opt[iF] .* d2C_dz2_opt[iZ], eachindex(PSF_zcoords)), eachindex(freqs))

    function hessian_vector_product!(Hv::AbstractVector, v::AbstractVector)
        v_T = @view v[1:length(free_T)]
        v_z = @view v[length(free_T)+1:end]

        # scatter v_T/v_z back into their free (non-active) pixels, leaving zeros at the active
        # (bound-pinned) positions -- mirrors Λ_T_2D/Λ_z_2D in Λ_dot_dFdo/Λ_dot_dFdo_distributed
        v_T_full = zeros(eltype(v_T), imghp.objN^2)
        v_T_full[free_T] .= v_T
        v_T_2D = reshape(v_T_full, imghp.objN, imghp.objN)

        v_z_full = zeros(eltype(v_z), imghp.objN^2)
        v_z_full[free_z] .= v_z
        v_z_2D = reshape(v_z_full, imghp.objN, imghp.objN)

        intermediate_TT = sum(
            iF -> sum(
                iZ -> begin
                    spectrum = v_T_2D .* dBdT_opt[iF] .* C_interp_3D[:, :, iZ]
                    weights[iF] .* convolve_with_fftPSF_at_freq_and_z(fftPSFs[iF, iZ], spectrum, plan_PSF)
                end,
                eachindex(PSF_zcoords)),
            eachindex(freqs))

        intermediate_Tz = sum(
            iF -> sum(
                iZ -> begin
                    spectrum = v_z_2D .* B_opt[iF] .* dC_dz_opt[iZ]
                    weights[iF] .* convolve_with_fftPSF_at_freq_and_z(fftPSFs[iF, iZ], spectrum, plan_PSF)
                end,
                eachindex(PSF_zcoords)),
            eachindex(freqs))

        # term_TT_1 (H_TT's Gauss-Newton part) and term_zT_1 both apply
        # convolveT_with_fftPSF_at_freq_and_z to the SAME intermediate_TT at every (iF, iZ) --
        # compute that convolution once per (iF, iZ) here and accumulate into both, instead of
        # two separate double loops each redoing it
        term_TT_1 = zeros(imghp.objN, imghp.objN)
        term_zT_1 = zeros(imghp.objN, imghp.objN)
        for iF in eachindex(freqs), iZ in eachindex(PSF_zcoords)
            conv = weights[iF] .* convolveT_with_fftPSF_at_freq_and_z(fftPSFs[iF, iZ], intermediate_TT, plan_PSF)
            term_TT_1 .+= (dBdT_opt[iF] .* C_interp_3D[:, :, iZ]) .* conv
            term_zT_1 .+= (B_opt[iF] .* dC_dz_opt[iZ]) .* conv
        end

        # same deal for term_Tz_1/term_zz_1: both apply convolveT_with_fftPSF_at_freq_and_z to
        # intermediate_Tz, so share that convolution across one loop instead of two
        term_Tz_1 = zeros(imghp.objN, imghp.objN)
        term_zz_1 = zeros(imghp.objN, imghp.objN)
        for iF in eachindex(freqs), iZ in eachindex(PSF_zcoords)
            conv = weights[iF] .* convolveT_with_fftPSF_at_freq_and_z(fftPSFs[iF, iZ], intermediate_Tz, plan_PSF)
            term_Tz_1 .+= (dBdT_opt[iF] .* C_interp_3D[:, :, iZ]) .* conv
            term_zz_1 .+= (dC_dz_opt[iZ] .* B_opt[iF]) .* conv
        end

        # residual-curvature terms are fully v-independent apart from the final .* v_T_2D /
        # .* v_z_2D, so reuse the sums precomputed above instead of recomputing any convolutions
        term_TT_2 = sum_TT_2 .* v_T_2D
        term_Tz_2 = sum_Tz_2 .* v_z_2D
        term_zT_2 = sum_Tz_2 .* v_T_2D
        term_zz_2 = sum_zz_2 .* v_z_2D

        # 2α .* v_T_full / 2β .* v_z_full are the regularization Hessians of the T-T and z-z
        # blocks (both zero at their respective active pixels, and both blocks are on the full
        # grid before compaction)
        Hv_TT = vec(2 .* term_TT_1) .+ vec(-2 .* term_TT_2) .+ 2α .* v_T_full
        Hv_Tz = vec(2 .* term_Tz_1) .+ vec(-2 .* term_Tz_2)

        Hv_zT = vec(2 .* term_zT_1) .+ vec(-2 .* term_zT_2)
        Hv_zz = vec(2 .* term_zz_1) .+ vec(-2 .* term_zz_2) .+ 2β .* v_z_full

        # combine each block (full objN x objN grid), then compact down to the free (non-active)
        # pixels -- mirrors the free_T/free_z scatter used to build v_T_2D/v_z_2D, in reverse
        Hv_T_full = Hv_TT .+ Hv_Tz
        Hv_z_full = Hv_zT .+ Hv_zz
        Hv[1:length(free_T)] .= Hv_T_full[free_T]
        Hv[length(free_T)+1:end] .= Hv_z_full[free_z]

        Hv
    end

    LinearMaps.LinearMap{Float64}(hessian_vector_product!, N; ismutating=true, issymmetric=true)
end