include("shared/mtktestset.jl")

@safetestset "FMI Extension Test" include("fmi/fmi.jl")
@safetestset "FMI Regressions" include("fmi/regressions.jl")
