using DelimitedFiles

# Regenerates objdata/MIT<N>.csv (Tmap source) in place, and generates a new objdata/MIT<N>_zmap.csv
# (depth-map source) alongside it, for every size already present in objdata/. Both stay in the
# same raw [0,1] convention LoadedObject (src/prepare.jl) expects, so they get scaled into physical
# units later via whatever Tlb/Tub/zlb_μm/zub_μm a LoadedObject instance is built with -- nothing
# about the scale is baked in here.
#
# Tmap: the existing MIT<N>.csv files have background = 0, letters = 1, with some anti-aliased
# edge pixels in between. Remapped via new = 0.5 + 0.5*old, so background moves to 0.5 (the
# midpoint) and letters stay at 1.0, preserving the original anti-aliasing gradient at the edges.
#
# zmap: the existing files carry no information about which pixels belong to which letter -- just a
# single combined mask. This is the geometric MIT logo mark, not plain letterforms: rendering
# MIT64.csv (see session notes) shows the "M" is actually 3 disjoint vertical bars, "I" is a short
# bar with a fainter (anti-aliased) lower tail, and "T" is a wide top bar over a narrower stem --
# so naive connected-component labeling of the whole mask doesn't give 3 blobs (M's 3 bars count
# separately). Confirmed instead, by projecting column-wise maxima, that every size from 16 up has
# exactly 5 disjoint column groups, consistently ordered M-bar1, M-bar2, M-bar3, I, T left to right
# -- so grouping columns (not pixels) and merging the first 3 groups into "M" recovers the letters
# directly from the existing data, no new source image needed. Each letter gets a single fixed
# depth level; background gets 0.5, matching the Tmap convention above. Unlike Tmap, boundary/anti-
# aliased pixels are not smoothly blended here -- zmap represents 4 genuinely discrete depths
# (background + M + I + T), so every nonzero raw pixel in a letter's column range is assigned that
# letter's single level.
#
# Size 8 is too coarse for this: the gaps between M's 3 bars vanish entirely (confirmed: only 3
# column groups survive, not 5), so there's no reliable way to recover which pixels are which
# letter at that resolution. zmap generation is skipped for size 8 (Tmap's background-shift still
# applies to it, since that doesn't depend on letter separation) -- not a practical loss, since
# production uses objN = 32.
const BACKGROUND_LEVEL = 0.5
const LETTER_LEVELS = (M = 0.0, I = 0.75, T = 1.0)
const EXPECTED_GROUPS = 5 # 3 M-bars + I + T

# Column j "belongs" to the glyph if any pixel in it is non-background; a pixel is itself
# foreground if its raw value is non-zero (background is always exactly 0.0 in the source data,
# confirmed across every size -- even the faintest anti-aliased edge pixel seen is 0.48, so there's
# no ambiguous middle ground to threshold against).
function column_groups(raw)
    colmax = [maximum(@view raw[:, j]) for j in axes(raw, 2)]
    support = colmax .> 0.01
    groups = UnitRange{Int}[]
    start = nothing
    for j in eachindex(support)
        if support[j] && start === nothing
            start = j
        elseif !support[j] && start !== nothing
            push!(groups, start:(j - 1))
            start = nothing
        end
    end
    start !== nothing && push!(groups, start:length(support))
    groups
end

function generate_zmap(raw, size_N)
    groups = column_groups(raw)
    if length(groups) != EXPECTED_GROUPS
        @warn "generate_MIT_logo($size_N): expected $EXPECTED_GROUPS column groups (3 M-bars + I + T), " *
              "found $(length(groups)) -- skipping zmap generation at this size"
        return nothing
    end

    letter_cols = (M = reduce(vcat, collect.(groups[1:3])), I = collect(groups[4]), T = collect(groups[5]))
    zmap = fill(BACKGROUND_LEVEL, size(raw))
    for letter in (:M, :I, :T)
        cols = getfield(letter_cols, letter)
        level = getfield(LETTER_LEVELS, letter)
        for j in cols, i in axes(raw, 1)
            raw[i, j] > 0.0 && (zmap[i, j] = level)
        end
    end
    zmap
end

function generate_MIT_logo(size_N)
    Tmap_path = "objdata/MIT$(size_N).csv"
    isfile(Tmap_path) || error("generate_MIT_logo: $Tmap_path doesn't exist")
    raw = Matrix{Float64}(readdlm(Tmap_path, ',', Float64))
    minimum(raw) == 0.0 ||
        error("generate_MIT_logo($size_N): minimum value in $Tmap_path is $(minimum(raw)), expected exactly " *
              "0.0 -- this file looks like it was already regenerated (background already moved off 0); " *
              "restore the original 0/1-background file before re-running this script")

    zmap = generate_zmap(raw, size_N)

    Tmap_new = BACKGROUND_LEVEL .+ (1 - BACKGROUND_LEVEL) .* raw
    writedlm(Tmap_path, Tmap_new, ',')
    println("MIT$(size_N): regenerated $Tmap_path (background -> $BACKGROUND_LEVEL)")

    if zmap !== nothing
        zmap_path = "objdata/MIT$(size_N)_zmap.csv"
        writedlm(zmap_path, zmap, ',')
        println("MIT$(size_N): wrote $zmap_path (levels M=$(LETTER_LEVELS.M) I=$(LETTER_LEVELS.I) " *
                "T=$(LETTER_LEVELS.T), background=$BACKGROUND_LEVEL)")
    end
end

sizes = sort([parse(Int, m.captures[1]) for f in readdir("objdata")
              for m in [match(r"^MIT(\d+)\.csv$", f)] if m !== nothing])
println("Regenerating MIT logos for sizes: ", sizes)
for size_N in sizes
    generate_MIT_logo(size_N)
end
