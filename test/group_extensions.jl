include("shared/mtktestset.jl")

@mtktestset("HomotopyContinuation Extension Test", "extensions/homotopy_continuation.jl")
@mtktestset("BifurcationKit Extension Test", "extensions/bifurcationkit.jl")
@mtktestset("Initialization maps AD", "extensions/initialization_maps_ad.jl")
# `ad.jl` is not included here: it runs via ModelingToolkitBase's Extensions group.
@mtktestset("AD Initials cotangent", "extensions/ad_initials.jl")
