function convolveT(out, kernel, plan)
    outL = size(out, 1)
    kerL = size(kernel, 1)
    inpL = kerL - outL

    out_pad = [zeros(inpL, inpL) zeros(inpL, outL); zeros(outL, inpL) out]
    arr_pad = plan * ((plan \ out_pad) .* kernel )

    arr = arr_pad[1:inpL, 1:inpL]
    arr
end

convolveT_with_fftPSF_at_freq_and_z(fftPSF_at_freq_and_z, image, plan_PSF) = real.(convolveT(image, fftPSF_at_freq_and_z, plan_PSF))

function dB_dT(T, freq, wavcen)
    A = (ħ * c * 10^6)/(wavcen * kB)
    numerator = (2 * A * freq^4 * exp(A * freq / T) ) 
    denominator = T^2 * (exp(A * freq / T) - 1)^2
    numerator/denominator
end

function df_δ_dz(z, smoothness_order)
    if smoothness_order == Inf
        return z > 0.0 ? exp(-1 / z) / z^2 : 0.0
    else
        return z > 0.0 ? (smoothness_order + 1)*z^smoothness_order : 0.0
    end
end

function dg_δ_dz(z, smoothness_order)
    g_denom = (f_δ(z, smoothness_order) + f_δ(1-z, smoothness_order))
    numerator_term1 = g_denom*df_δ_dz(z, smoothness_order)
    numerator_term2 = f_δ(z, smoothness_order)*(df_δ_dz(z, smoothness_order) - df_δ_dz(1-z, smoothness_order))
    denominator = g_denom ^2
    (numerator_term1 - numerator_term2)/denominator
end

function dδ_dz(z, smoothness_order, Δz)
    (1/Δz)*((-1/Δz)*dg_δ_dz(z/Δz, smoothness_order) + (1/Δz)*dg_δ_dz(-z/Δz, smoothness_order))
end

function dC_dz(z, PSF_zcoord, smoothness_order, Δz)
    diff = PSF_zcoord - z
    -dδ_dz(diff, smoothness_order, Δz)
end

function get_fftPSFs_from_precomputed(geoms, fftPSFs_precomputed, freqs, incidents, surrogates, n2f_kernels, plans_n2f, plans_PSF, php, imghp)
    fftPSFs_precomputed
end

function ChainRulesCore.rrule(::typeof(get_fftPSFs_from_precomputed), geoms, fftPSFs_precomputed, freqs, incidents, surrogates, n2f_kernels, plans_n2f, plans_PSF, php, imghp)
    PSF_zlen = imghp.PSF_zlen
    nF = length(freqs)
    inds = CartesianIndices((1:nF, 1:PSF_zlen))

    function pullback(ΔfftPSFs)
        Δgeoms_per_thread = [zeros(size(geoms)) for _ in 1:nthreads()]
        @threads for idx in eachindex(inds)
            tid = threadid()
            I = inds[idx]
            iF, iZ = I[1], I[2]
            Δgeoms_local = Zygote.gradient(
                g -> real(dot(End2EndThermalImg.get_fftPSF(End2EndThermalImg.get_PSF_at_freq_and_z(
                    freqs[iF], incidents[iF, iZ], surrogates[iF], g,
                    n2f_kernels[iF], plans_n2f[tid], php, imghp
                ), plans_PSF[tid]), ΔfftPSFs[iF, iZ])),
                geoms
            )[1]
            Δgeoms_per_thread[tid] .+= Δgeoms_local
        end
        Δgeoms = sum(Δgeoms_per_thread)
        return NoTangent(), Δgeoms, NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent()
    end

    return fftPSFs_precomputed, pullback
end

# TODO: should convolutions be computed only once?
# function Λ_dot_dFdo(Λ, object, object_opt, fftPSFs_precomputed, freqs, incidents, surrogates, geoms, α, β, noise, noise_scale, n2f_kernels, plans_n2f, plans_PSF, weights, php, imghp, rechp)
#     PSF_zcoords = ignore_derivatives( ()-> get_PSF_zcoords(imghp))
#     δ_Δz = ignore_derivatives( ()-> get_discretized_δ_function(imghp.smoothness_order, imghp.PSF_Δz))
#     C_opt = ignore_derivatives( ()-> get_C_interp_3D(object_opt.zmap, PSF_zcoords, δ_Δz, imghp))
#     dB_dT_opt = ignore_derivatives( ()-> [dB_dT.(object_opt.Tmap, freq, php.wavcen) for freq in freqs])
#     B_opt = ignore_derivatives( ()-> get_black_body_spectrum(object_opt.Tmap, php))
#     dC_dz_opt = ignore_derivatives( ()-> [dC_dz.(object_opt.zmap, PSF_zcoord, imghp.smoothness_order, imghp.PSF_Δz) for PSF_zcoord in PSF_zcoords])

#     fftPSFs = get_fftPSFs_from_precomputed(geoms, fftPSFs_precomputed, freqs, incidents, surrogates, n2f_kernels, plans_n2f, plans_PSF, php, imghp)
#     plan_PSF = plans_PSF[1]
#     v = make_noisy_image_from_3D(object, fftPSFs, weights, noise, noise_scale, plan_PSF, php, imghp)
#     image_opt = make_image_from_3D(object_opt, fftPSFs, weights, plan_PSF, php, imghp)
#     image_diff = v .- image_opt

#     # first compute Λ_T dot df/dT
#     Λ_T_2D = ignore_derivatives(() -> reshape(Λ[1:imghp.objN^2], imghp.objN, imghp.objN))

#     term1 = sum(
#     iF -> begin
#         ΛdB = ignore_derivatives(() -> Λ_T_2D .* dB_dT_opt[iF])
#         sum(
#             iZ -> -2 * weights[iF] * dot(ΛdB .* C_opt[:,:,iZ], convolveT_with_fftPSF_at_freq_and_z(fftPSFs[iF, iZ], image_diff, plan_PSF)),
#             eachindex(PSF_zcoords))
#     end,
#     eachindex(freqs))
    
#     # then compute Λ_z dot df/dz
#     Λ_z_2D = ignore_derivatives(() -> reshape(Λ[imghp.objN^2 + 1:end], imghp.objN, imghp.objN))
#     term2 = sum(
#     iF -> begin
#         ΛB = ignore_derivatives(() -> Λ_z_2D .* B_opt[iF])
#         sum(
#             iZ -> -2 * weights[iF] * dot(ΛB .* dC_dz_opt[iZ], convolveT_with_fftPSF_at_freq_and_z(fftPSFs[iF, iZ], image_diff, plan_PSF)),
#             eachindex(PSF_zcoords))
#     end,
#     eachindex(freqs))

#     term1 + term2
# end

function Λ_dot_dFdo(Λ, object, object_opt, fftPSFs_precomputed, freqs, incidents, surrogates, geoms, α, β, noise, noise_scale, n2f_kernels, plans_n2f, plans_PSF, weights, php, imghp, rechp)
    PSF_zcoords = ignore_derivatives(() -> get_PSF_zcoords(imghp))
    δ_Δz = ignore_derivatives(() -> get_discretized_δ_function(imghp.smoothness_order, imghp.PSF_Δz))
    C_opt = ignore_derivatives(() -> get_C_interp_3D(object_opt.zmap, PSF_zcoords, δ_Δz, imghp))
    dB_dT_opt = ignore_derivatives(() -> [dB_dT.(object_opt.Tmap, freq, php.wavcen) for freq in freqs])
    B_opt = ignore_derivatives(() -> get_black_body_spectrum(object_opt.Tmap, php))
    dC_dz_opt = ignore_derivatives(() -> [dC_dz.(object_opt.zmap, PSF_zcoord, imghp.smoothness_order, imghp.PSF_Δz) for PSF_zcoord in PSF_zcoords])
    fftPSFs = get_fftPSFs_from_precomputed(geoms, fftPSFs_precomputed, freqs, incidents, surrogates, n2f_kernels, plans_n2f, plans_PSF, php, imghp)
    plan_PSF = plans_PSF[1]
    v = make_noisy_image_from_3D(object, fftPSFs, weights, noise, noise_scale, plan_PSF, php, imghp)
    image_opt = make_image_from_3D(object_opt, fftPSFs, weights, plan_PSF, php, imghp)
    image_diff = v .- image_opt
    Λ_T_2D = ignore_derivatives(() -> reshape(Λ[1:imghp.objN^2], imghp.objN, imghp.objN))
    Λ_z_2D = ignore_derivatives(() -> reshape(Λ[imghp.objN^2 + 1:end], imghp.objN, imghp.objN))
    sum(
        iF -> begin
            ΛdB = ignore_derivatives(() -> Λ_T_2D .* dB_dT_opt[iF])
            ΛB  = ignore_derivatives(() -> Λ_z_2D .* B_opt[iF])
            sum(
                iZ -> begin
                    Λweight = ignore_derivatives(() -> ΛdB .* C_opt[:,:,iZ] .+ ΛB .* dC_dz_opt[iZ])
                    conv = convolveT_with_fftPSF_at_freq_and_z(fftPSFs[iF, iZ], image_diff, plan_PSF)
                    -2 * weights[iF] * dot(Λweight, conv)
                end,
                eachindex(PSF_zcoords))
        end,
        eachindex(freqs))
end

function get_active_variables(zmap, imghp)
    cartesian_indices = findall(z -> z == imghp.PSF_zlb || z == imghp.PSF_zub, zmap)
    linear_indices = LinearIndices(zmap)[cartesian_indices]
    return cartesian_indices, linear_indices
end