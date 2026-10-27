const DEFAULT_FIGSIZE_X = 3
const DEFAULT_FIGSIZE_Y = 2
const DEFAULT_FONTSIZE = 12

function plot_geoms(geoms_flat, php; figsize_x = DEFAULT_FIGSIZE_X, figsize_y = DEFAULT_FIGSIZE_Y, fontsize = DEFAULT_FONTSIZE, cmap = ColorMap(colorschemes[:grays].colors), savefile = "")
    geoms = unflatten_square_matrix(geoms_flat)
    matplotlib.rcParams["font.size"] = fontsize
    fig, ax = subplots(figsize=(figsize_x, figsize_y))
    im = ax.imshow(geoms, cmap=cmap, vmin=php.pillar_width_lb, vmax=php.pillar_width_ub)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Pillar widths (unitless)")
    ax.axis("off")
    fig.tight_layout()
    
    if !isempty(savefile)
        fig.savefig(savefile)
    end
    plotclose()
    return fig, ax
end

function plot_objective_function_history(objective_history; figsize_x = DEFAULT_FIGSIZE_X, figsize_y = DEFAULT_FIGSIZE_Y, fontsize = DEFAULT_FONTSIZE, savefile = "")
    matplotlib.rcParams["font.size"] = fontsize
    fig, ax = subplots(figsize=(figsize_x, figsize_y))
    ax.plot(1:length(objective_history), objective_history)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Objective function")
    fig.tight_layout()
    if !isempty(savefile)
        fig.savefig(savefile)
    end
    plotclose()
    return fig, ax
end

function plot_alpha_beta_history(αs, βs; figsize_x = DEFAULT_FIGSIZE_X, figsize_y = DEFAULT_FIGSIZE_Y, fontsize = DEFAULT_FONTSIZE, log_scale = false, savefile = "")
    matplotlib.rcParams["font.size"] = fontsize
    fig, axes = subplots(1, 2; figsize=(2 * figsize_x, figsize_y))

    # log_scale: α/β can span many orders of magnitude within a single history (e.g. a
    # resolve_initial_α_β search sweeping 1e-7 to 1e-5 and beyond) -- linear axes would be close
    # to unreadable there. Defaults to false so train_outer_optimization_distributed's own, more
    # gradually-changing Adam-driven history keeps its existing linear appearance.
    if log_scale
        axes[0].semilogy(1:length(αs), αs, marker="o", color="tab:blue")
        axes[1].semilogy(1:length(βs), βs, marker="o", color="tab:orange")
    else
        axes[0].plot(1:length(αs), αs, marker="o", color="tab:blue")
        axes[1].plot(1:length(βs), βs, marker="o", color="tab:orange")
    end
    axes[0].set_xlabel("Iteration")
    axes[0].set_ylabel("α")
    axes[0].grid(true, which = log_scale ? "both" : "major", alpha=0.3)

    axes[1].set_xlabel("Iteration")
    axes[1].set_ylabel("β")
    axes[1].grid(true, which = log_scale ? "both" : "major", alpha=0.3)

    fig.tight_layout()
    if !isempty(savefile)
        fig.savefig(savefile)
    end
    plotclose()
    return fig, axes
end

function compute_PSF_grid_dimensions(num_PSFs)
    factors = []
    for i in 1:num_PSFs
        if num_PSFs % i == 0
            push!(factors, i)
        end
    end
    # Find the factor pair closest to square
    best_diff = Inf
    numy, numx = 1, num_PSFs
    for factor in factors
        other_factor = num_PSFs ÷ factor
        diff = abs(factor - other_factor)
        if diff < best_diff
            best_diff = diff
            numy, numx = min(factor, other_factor), max(factor, other_factor)
        end
    end
    numy, numx
end

unitless_z_to_meters(z, php) = z * php.wavcen / 10^6

# TODO: convert z coordinates to meters (?) actually: think i'm going to leave unitless
function plot_PSFs_at_z(PSFs_at_z, zcoord, freqs; figsize_x = DEFAULT_FIGSIZE_X, figsize_y = DEFAULT_FIGSIZE_Y, fontsize = DEFAULT_FONTSIZE, cmap = "viridis", savefile = "", vmin = nothing, vmax = nothing, extra_title = "", log_scale = true)
    matplotlib.rcParams["font.size"] = fontsize
    num_PSFs = length(PSFs_at_z)
    numy, numx = compute_PSF_grid_dimensions(num_PSFs)
    if isnothing(vmin) && isnothing(vmax)
        vmin = minimum([minimum(PSF) for PSF in PSFs_at_z])
        vmax = maximum([maximum(PSF) for PSF in PSFs_at_z])
    end

    fig, axes = subplots(numy, numx, figsize=(numx*figsize_x, numy*figsize_y), constrained_layout=true)
    
    im = nothing
    for (i, PSF) in enumerate(PSFs_at_z)
        if num_PSFs == 1
            ax = axes
        elseif numy == 1 || numx == 1
            ax = axes[i-1]
        else
            row = (i-1) ÷ numx
            col = (i-1) % numx
            ax = axes[row, col]
        end
        norm = log_scale ? matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax) : nothing
        im = ax.imshow(PSF, norm=norm, cmap=cmap)
        ax.axis("off")
        ax.set_title(L"\nu = %$(round(freqs[i],digits=2))", fontsize = fontsize)
    end
    
    cbar = fig.colorbar(im, ax=axes, location="right", shrink=0.9)
    fig.text(0.5, 1.05, "z = $zcoord $(extra_title)", ha="center", va="top", transform=fig.transFigure, fontsize = fontsize)
    if !isempty(savefile)
        fig.savefig(savefile, bbox_inches="tight")
    end
    plotclose()
    return fig, axes
end

function plot_object(object; figsize_x = DEFAULT_FIGSIZE_X,
                                     figsize_y = DEFAULT_FIGSIZE_Y,
                                     fontsize = 10)
    matplotlib.rcParams["font.size"] = fontsize
    fig, ax = subplots(2, 1, figsize=(figsize_x, 2*figsize_y))
    im0 = ax[0].imshow(object.Tmap, cmap = "magma")
    cbar = fig.colorbar(im0, ax=ax[0], label="T (Kelvin)")
    ax[0].set_title(L"$T(x,y)$")
    ax[0].axis("off")

    im1 = ax[1].imshow(object.zmap, cmap = ColorMap(colorschemes[:devon].colors))
    cbar = fig.colorbar(im1, ax=ax[1], label="z (unitless)")
    ax[1].set_title(L"$z(x,y)$")
    ax[1].axis("off")

    fig.tight_layout()
    plotclose()
    return fig, ax
end

function image_axis!(ax)
    ax.set_xticks(Float64[])
    ax.set_yticks(Float64[])
end

function plot_reconstruction_results(object, object_opt;
                                     figsize_x = DEFAULT_FIGSIZE_X,
                                     figsize_y = DEFAULT_FIGSIZE_Y,
                                     fontsize = 10,
                                     extra_title = "",
                                     savefile = "")

    matplotlib.rcParams["font.size"] = fontsize
    fig, ax = subplots(2, 3, figsize=(3*figsize_x, 2*figsize_y))

    Tmin = minimum(object.Tmap)
    Tmax = maximum(object.Tmap)
    zmin = minimum(object.zmap)
    zmax = maximum(object.zmap)

    im00 = ax[0,0].imshow(object.Tmap, cmap = "magma", vmin = Tmin, vmax = Tmax)
    ax[0,0].set_title(L"$T(x,y)$ [ground truth]")
    fig.colorbar(im00, ax=ax[0,0], label="T (Kelvin)")
    image_axis!(ax[0,0])

    im10 = ax[1,0].imshow(object.zmap, cmap = ColorMap(colorschemes[:devon].colors), vmin = zmin, vmax = zmax)
    ax[1,0].set_title(L"$z(x,y)$ [ground truth]")
    fig.colorbar(im10, ax=ax[1,0], label="z (unitless)")
    image_axis!(ax[1,0])

    im01 = ax[0,1].imshow(object_opt.Tmap, cmap = "magma", vmin = Tmin, vmax = Tmax)
    ax[0,1].set_title(L"$T_{est}(x,y)$ [reconstrcted]")
    fig.colorbar(im01, ax=ax[0,1], label="T (Kelvin)")
    image_axis!(ax[0,1])

    im11 = ax[1,1].imshow(object_opt.zmap, cmap = ColorMap(colorschemes[:devon].colors), vmin = zmin, vmax = zmax)
    ax[1,1].set_title(L"$z_{est}(x,y)$ [reconstructed]")
    fig.colorbar(im11, ax=ax[1,1], label="z (unitless)")
    image_axis!(ax[1,1])

    # per-pixel signed relative error, as a percentage of the true value at that pixel -- safe
    # here since neither Tmap nor zmap ever approaches zero (both are generated with lb/ub bounds
    # bracketed away from 0; see e.g. UniformlyRandomObject/GaussianObject, src/prepare.jl).
    error_T = 100 .* (object_opt.Tmap .- object.Tmap) ./ object.Tmap
    error_z = 100 .* (object_opt.zmap .- object.zmap) ./ object.zmap

    # relative L2 error -- ||opt - true|| / ||true||, one shared scalar per map (not a pointwise
    # ratio at each pixel like error_T/error_z above). relerr_T^2 + relerr_z^2 is exactly
    # L_reconstruction_error (src/optimize.jl) -- the actual scalar the reconstruction and the
    # outer bilevel objective are built around.
    relerr_T = norm(object_opt.Tmap .- object.Tmap) / norm(object.Tmap)
    relerr_z = norm(object_opt.zmap .- object.zmap) / norm(object.zmap)

    elim_T = maximum(abs, error_T)
    im02 = ax[0,2].imshow(error_T, cmap = "RdBu", vmin = -elim_T, vmax = elim_T)
    ax[0,2].set_title("per-pixel relative error (%)\nrelative L2 error: $(round(100*relerr_T, digits=2))%")
    fig.colorbar(im02, ax=ax[0,2], label="%")
    image_axis!(ax[0,2])

    elim_z = maximum(abs, error_z)
    im12 = ax[1,2].imshow(error_z, cmap = "RdBu", vmin = -elim_z, vmax = elim_z)
    ax[1,2].set_title("per-pixel relative error (%)\nrelative L2 error: $(round(100*relerr_z, digits=2))%")
    fig.colorbar(im12, ax=ax[1,2], label="%")
    image_axis!(ax[1,2])

    fig.text(0.5, 1.05, "$(extra_title)", ha="center", va="top", transform=fig.transFigure, fontsize = fontsize + 2)
    fig.tight_layout()
    if !isempty(savefile)
        fig.savefig(savefile, bbox_inches="tight")
    end
    plotclose()
    return fig, ax
end

function plot_image(image, imghp; figsize_x = DEFAULT_FIGSIZE_X, figsize_y = DEFAULT_FIGSIZE_Y, fontsize = DEFAULT_FONTSIZE, scale_image_to_object = false, extra_title = "")
    matplotlib.rcParams["font.size"] = fontsize
    image_scale_factor = scale_image_to_object ? imghp.imgN / imghp.objN : 1
    fig, ax = subplots(1, 1, figsize=(image_scale_factor*figsize_x, image_scale_factor*figsize_y))

    im = ax.imshow(image)
    ax.set_title("Image $(extra_title)")
    ax.axis("off")
    cbar = fig.colorbar(im, ax=ax, cmap = "viridis")

    fig.tight_layout()
    plotclose()
    return fig, ax
end