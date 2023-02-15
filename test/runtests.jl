using SafeTestsets, Test

@safetestset "AliasGraph Test" begin include("alias.jl") end
@safetestset "Pantelides Test" begin include("pantelides.jl") end
@safetestset "Linear Algebra Test" begin include("linalg.jl") end
@safetestset "AbstractSystem Test" begin include("abstractsystem.jl") end
@safetestset "Variable Scope Tests" begin include("variable_scope.jl") end
@safetestset "Symbolic Parameters Test" begin include("symbolic_parameters.jl") end
@safetestset "Parsing Test" begin include("variable_parsing.jl") end
@safetestset "Simplify Test" begin include("simplify.jl") end
@safetestset "Direct Usage Test" begin include("direct.jl") end
@safetestset "System Linearity Test" begin include("linearity.jl") end
@safetestset "Linearization Tests" begin include("linearize.jl") end
@safetestset "Input Output Test" begin include("input_output_handling.jl") end
@safetestset "Clock Test" begin include("clock.jl") end
@safetestset "DiscreteSystem Test" begin include("discretesystem.jl") end
@safetestset "ODESystem Test" begin include("odesystem.jl") end
@safetestset "Unitful Quantities Test" begin include("units.jl") end
@safetestset "LabelledArrays Test" begin include("labelledarrays.jl") end
@safetestset "Mass Matrix Test" begin include("mass_matrix.jl") end
@safetestset "SteadyStateSystem Test" begin include("steadystatesystems.jl") end
@safetestset "SDESystem Test" begin include("sdesystem.jl") end
@safetestset "NonlinearSystem Test" begin include("nonlinearsystem.jl") end
@safetestset "PDE Construction Test" begin include("pde.jl") end
@safetestset "JumpSystem Test" begin include("jumpsystem.jl") end
@safetestset "Constraints Test" begin include("constraints.jl") end
@safetestset "Reduction Test" begin include("reduction.jl") end
@safetestset "ODAEProblem Test" begin include("odaeproblem.jl") end
@safetestset "Components Test" begin include("components.jl") end
@safetestset "print_tree" begin include("print_tree.jl") end
@safetestset "Error Handling" begin include("error_handling.jl") end
@safetestset "StructuralTransformations" begin include("structural_transformation/runtests.jl") end
@safetestset "State Selection Test" begin include("state_selection.jl") end
@safetestset "Symbolic Event Test" begin include("symbolic_events.jl") end
@safetestset "Stream Connnect Test" begin include("stream_connectors.jl") end
@safetestset "Lowering Integration Test" begin include("lowering_solving.jl") end
@safetestset "Test Big System Usage" begin include("bigsystem.jl") end
@safetestset "Depdendency Graph Test" begin include("dep_graphs.jl") end
@safetestset "Function Registration Test" begin include("function_registration.jl") end
@safetestset "Precompiled Modules Test" begin include("precompile_test.jl") end
@testset "Distributed Test" begin include("distributed.jl") end
@safetestset "Variable Utils Test" begin include("variable_utils.jl") end
@safetestset "Variable Metadata Test" begin include("test_variable_metadata.jl") end
@safetestset "DAE Jacobians Test" begin include("dae_jacobian.jl") end
@safetestset "Jacobian Sparsity" begin include("jacobiansparsity.jl") end
println("Last test requires gcc available in the path!")
@safetestset "C Compilation Test" begin include("ccompile.jl") end
@testset "Serialization" begin include("serialization.jl") end
@safetestset "Modelingtoolkitize Test" begin include("modelingtoolkitize.jl") end
@safetestset "OptimizationSystem Test" begin include("optimizationsystem.jl") end
@safetestset "FuncAffect Test" begin include("funcaffect.jl") end
@safetestset "Constants Test" begin include("constants.jl") end
# Reference tests go Last
@safetestset "Latexify recipes Test" begin include("latexify.jl") end
