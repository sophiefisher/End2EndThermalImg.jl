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
function plot_PSFs_at_z(PSFs_at_z, zcoord, freqs; figsize_x = DEFAULT_FIGSIZE_X, figsize_y = DEFAULT_FIGSIZE_Y, fontsize = DEFAULT_FONTSIZE, cmap = "viridis", savefile = "", vmin = nothing, vmax = nothing, extra_title = "")
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
        im = ax.imshow(PSF, norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), cmap=cmap)
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

function plot_reconstruction_results(object, object_opt;
                                     figsize_x = DEFAULT_FIGSIZE_X,
                                     figsize_y = DEFAULT_FIGSIZE_Y,
                                     fontsize = 10,
                                     extra_title = "")

    matplotlib.rcParams["font.size"] = fontsize
    fig, ax = subplots(2, 3, figsize=(3*figsize_x, 2*figsize_y))
    
    Tmin = minimum(object.Tmap)
    Tmax = maximum(object.Tmap)
    zmin = minimum(object.zmap)
    zmax = maximum(object.zmap)

    im00 = ax[0,0].imshow(object.Tmap, cmap = "magma", vmin = Tmin, vmax = Tmax)
    ax[0,0].set_title(L"$T(x,y)$ [ground truth]")
    fig.colorbar(im00, ax=ax[0,0])
    ax[0,0].axis("off")

    im10 = ax[1,0].imshow(object.zmap, cmap = ColorMap(colorschemes[:devon].colors), vmin = zmin, vmax = zmax)
    ax[1,0].set_title(L"$z(x,y)$ [ground truth]")
    fig.colorbar(im10, ax=ax[1,0])
    ax[1,0].axis("off")

    im01 = ax[0,1].imshow(object_opt.Tmap, cmap = "magma", vmin = Tmin, vmax = Tmax)
    ax[0,1].set_title(L"$T_{est}(x,y)$ [reconstrcted]")
    fig.colorbar(im01, ax=ax[0,1])
    ax[0,1].axis("off")

    im11 = ax[1,1].imshow(object_opt.zmap, cmap = ColorMap(colorschemes[:devon].colors), vmin = zmin, vmax = zmax)
    ax[1,1].set_title(L"$z_{est}(x,y)$ [reconstructed]")
    fig.colorbar(im11, ax=ax[1,1])
    ax[1,1].axis("off")

    im02 = ax[0,2].imshow((object.Tmap .- object_opt.Tmap).^2 ./ object.Tmap.^2, cmap = ColorMap(colorschemes[:grays].colors), vmin = 0)
    ax[0,2].set_title(L"Relative square error $\frac{(T_i - T_{est_i})^2}{T_i^2}$")
    fig.colorbar(im02, ax=ax[0,2])
    ax[0,2].axis("off")

    im12 = ax[1,2].imshow((object.zmap .- object_opt.zmap).^2 ./ object.zmap.^2, cmap = ColorMap(colorschemes[:grays].colors), vmin = 0)
    ax[1,2].set_title(L"Relative square error $\frac{(z_i - z_{est_i})^2}{z_i^2}$")
    fig.colorbar(im12, ax=ax[1,2])
    ax[1,2].axis("off")

    MSE_T = sum((object.Tmap .- object_opt.Tmap).^2) / sum(object.Tmap.^2)
    MSE_z = sum((object.zmap .- object_opt.zmap).^2) / sum(object.zmap.^2)

    ax[0,2].text(1.6, 0.5, "Relative MSE \n = $(round(MSE_T, digits=4)) \n \n Relative RMSE \n = $(round(sqrt(MSE_T)*100, digits=4))%",
                 transform = ax[0,2].transAxes,
                 va = "center", ha = "left", fontsize = fontsize)

    ax[1,2].text(1.6, 0.5, "Relative MSE \n = $(round(MSE_z, digits=4)) \n \n Relative RMSE \n = $(round(sqrt(MSE_z)*100, digits=4))%",

                 transform = ax[1,2].transAxes,
                 va = "center", ha = "left", fontsize = fontsize)

    fig.text(0.5, 1.05, "$(extra_title)", ha="center", va="top", transform=fig.transFigure, fontsize = fontsize + 2)
    fig.tight_layout()
    plotclose()
    return fig, ax
end

function plot_noisy_image(noisy_image, imghp; figsize_x = DEFAULT_FIGSIZE_X, figsize_y = DEFAULT_FIGSIZE_Y, fontsize = DEFAULT_FONTSIZE)
    matplotlib.rcParams["font.size"] = fontsize
    image_scale_factor = imghp.imgN / imghp.objN
    fig, ax = subplots(1, 1, figsize=(image_scale_factor*figsize_x, image_scale_factor*figsize_y))

    im = ax.imshow(noisy_image)
    ax.set_title("Noisy image")
    ax.axis("off")
    cbar = fig.colorbar(im, ax=ax, cmap = "viridis")

    fig.tight_layout()
    plotclose()
    return fig, ax
end