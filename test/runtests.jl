using ModelingToolkit, Test

@testset "Parsing Test" begin include("variable_parsing.jl") end
@testset "Differentiation Test" begin include("derivatives.jl") end
@testset "Simplify Test" begin include("simplify.jl") end
@testset "System Construction Test" begin include("system_construction.jl") end
