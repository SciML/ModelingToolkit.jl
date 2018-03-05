using SciCompDSL, Base.Test

@testset "Parsing Test" begin include("variable_parsing.jl") end
@testset "System Construction Test" begin include("system_construction.jl") end
@testset "Default Values" begin include("default_values.jl") end
