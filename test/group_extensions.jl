include("shared/mtktestset.jl")

@mtktestset("HomotopyContinuation Extension Test", "extensions/homotopy_continuation.jl")
@mtktestset("BifurcationKit Extension Test", "extensions/bifurcationkit.jl")
@mtktestset("Initialization maps AD", "extensions/initialization_maps_ad.jl")
# Parent `@mtktestset("Auto Differentiation Test", "extensions/ad.jl")` stays
# disabled: that file is re-enabled in ModelingToolkitBase's Extensions group,
# and `@mtktestset` would re-run the same include (duplicate). The Initials
# testset that needed full MTK `mtkcompile` was not moved here — it still fails
# under ModelingToolkit (see job REPORT); leaving it out of this change.
