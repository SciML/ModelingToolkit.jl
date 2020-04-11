using ModelingToolkit, Test

@testset "Parsing Test" begin include("variable_parsing.jl") end
@testset "Differentiation Test" begin include("derivatives.jl") end
@testset "Simplify Test" begin include("simplify.jl") end
@testset "Operation Overloads Test" begin include("operation_overloads.jl") end
@testset "Direct Usage Test" begin include("direct.jl") end
@testset "ODESystem Test" begin include("odesystem.jl") end
@testset "Mass Matrix Test" begin include("mass_matrix.jl") end
@testset "SDESystem Test" begin include("sdesystem.jl") end
@testset "NonlinearSystem Test" begin include("nonlinearsystem.jl") end
@testset "OptimizationSystem Test" begin include("optimizationsystem.jl") end
@testset "Build Targets Test" begin include("build_targets.jl") end
@testset "Domain Test" begin include("domains.jl") end
@testset "Constraints Test" begin include("constraints.jl") end
@testset "PDE Construction Test" begin include("pde.jl") end
@testset "Distributed Test" begin include("distributed.jl") end
#@testset "Latexify recipes Test" begin include("latexify.jl") end
