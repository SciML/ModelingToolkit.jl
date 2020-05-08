using SafeTestsets, Test

@safetestset "Parsing Test" begin include("variable_parsing.jl") end
@safetestset "Differentiation Test" begin include("derivatives.jl") end
@safetestset "Simplify Test" begin include("simplify.jl") end
@safetestset "Operation Overloads Test" begin include("operation_overloads.jl") end
@safetestset "Direct Usage Test" begin include("direct.jl") end
@safetestset "Build Function Test" begin include("build_function.jl") end
@safetestset "ODESystem Test" begin include("odesystem.jl") end
@safetestset "Mass Matrix Test" begin include("mass_matrix.jl") end
@safetestset "SDESystem Test" begin include("sdesystem.jl") end
@safetestset "NonlinearSystem Test" begin include("nonlinearsystem.jl") end
@safetestset "OptimizationSystem Test" begin include("optimizationsystem.jl") end
@safetestset "ReactionSystem Test" begin include("reactionsystem.jl") end
@safetestset "Build Targets Test" begin include("build_targets.jl") end
@safetestset "Domain Test" begin include("domains.jl") end
@safetestset "Constraints Test" begin include("constraints.jl") end
@safetestset "PDE Construction Test" begin include("pde.jl") end
@safetestset "Lowering Integration Test" begin include("lowering_solving.jl") end
@safetestset "Test Big System Usage" begin include("bigsystem.jl") end
#@testset "Latexify recipes Test" begin include("latexify.jl") end
@testset "Distributed Test" begin include("distributed.jl") end
