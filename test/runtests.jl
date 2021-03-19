using SafeTestsets, Test

@safetestset "Symbolic parameters test" begin include("symbolic_parameters.jl") end
@safetestset "Parsing Test" begin include("variable_parsing.jl") end
@safetestset "Differentiation Test" begin include("derivatives.jl") end
@safetestset "Simplify Test" begin include("simplify.jl") end
@safetestset "Operation Overloads Test" begin include("operation_overloads.jl") end
@safetestset "Direct Usage Test" begin include("direct.jl") end
@safetestset "System Linearity Test" begin include("linearity.jl") end
@safetestset "ODESystem Test" begin include("odesystem.jl") end
@safetestset "LabelledArrays Test" begin include("labelledarrays.jl") end
@safetestset "Mass Matrix Test" begin include("mass_matrix.jl") end
@safetestset "SteadyStateSystem Test" begin include("steadystatesystems.jl") end
@safetestset "SDESystem Test" begin include("sdesystem.jl") end
@safetestset "NonlinearSystem Test" begin include("nonlinearsystem.jl") end
@safetestset "OptimizationSystem Test" begin include("optimizationsystem.jl") end
@safetestset "ReactionSystem Test" begin include("reactionsystem.jl") end
@safetestset "JumpSystem Test" begin include("jumpsystem.jl") end
@safetestset "ControlSystem Test" begin include("controlsystem.jl") end
@safetestset "Domain Test" begin include("domains.jl") end
@safetestset "Modelingtoolkitize Test" begin include("modelingtoolkitize.jl") end
@safetestset "Constraints Test" begin include("constraints.jl") end
@safetestset "Reduction Test" begin include("reduction.jl") end
@safetestset "Components Test" begin include("components.jl") end
@safetestset "PDE Construction Test" begin include("pde.jl") end
@safetestset "Lowering Integration Test" begin include("lowering_solving.jl") end
@safetestset "Test Big System Usage" begin include("bigsystem.jl") end
@safetestset "Depdendency Graph Test" begin include("dep_graphs.jl") end
@safetestset "Function Registration Test" begin include("function_registration.jl") end
@safetestset "Precompiled Modules Test" begin include("precompile_test.jl") end
@testset "Distributed Test" begin include("distributed.jl") end
@safetestset "Variable Utils Test" begin include("variable_utils.jl") end
println("Last test requires gcc available in the path!")
@safetestset "C Compilation Test" begin include("ccompile.jl") end
@safetestset "Latexify recipes Test" begin include("latexify.jl") end
@safetestset "StructuralTransformations" begin include("structural_transformation/runtests.jl") end
@testset "Serialization" begin include("serialization.jl") end
@testset "print_tree" begin include("print_tree.jl") end
