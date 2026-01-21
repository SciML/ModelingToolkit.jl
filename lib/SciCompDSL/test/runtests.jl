using SafeTestsets, Pkg, Test

const MTKBasePath = joinpath(dirname(dirname(@__DIR__)), "ModelingToolkitBase")
const MTKBasePkgSpec = PackageSpec(; path = MTKBasePath)

const MTKPath = dirname(dirname(dirname(@__DIR__)))
const MTKPkgSpec = PackageSpec(; path = MTKPath)

Pkg.develop([MTKBasePkgSpec, MTKPkgSpec])

@safetestset "Model parsing - MTKBase" include("model_parsing.jl")
@safetestset "Model parsing - MTK" begin
    using ModelingToolkit
    import ModelingToolkitBase
    include("model_parsing.jl")
end
