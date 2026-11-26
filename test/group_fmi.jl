include("shared/mtktestset.jl")

@safetestset "FMI Regressions" include("fmi/regressions.jl")
@safetestset "FMI Extension Test" include("fmi/fmi.jl")
@safetestset "FMI Reference FMUs" include("fmi/reference_fmus.jl")
