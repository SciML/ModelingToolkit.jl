include("shared/mtktestset.jl")

@mtktestset("HomotopyContinuation Extension Test", "extensions/homotopy_continuation.jl")
@mtktestset("BifurcationKit Extension Test", "extensions/bifurcationkit.jl")
@mtktestset("Initialization maps AD", "extensions/initialization_maps_ad.jl")
# @mtktestset("Auto Differentiation Test", "extensions/ad.jl")
