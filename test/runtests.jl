using ModelingToolkit, Test

@testset "Parsing Test" begin include("variable_parsing.jl") end
@testset "Basic Variables and Operations" begin include("basic_variables_and_operations.jl") end
@testset "Differentiation Test" begin include("derivatives.jl") end
@testset "Domain Test" begin include("domains.jl") end
@testset "Simplify Test" begin include("simplify.jl") end
@testset "Components Test" begin include("components.jl") end
@testset "System Construction Test" begin include("system_construction.jl") end
