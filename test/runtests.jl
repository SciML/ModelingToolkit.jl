using SafeTestsets, Pkg, Test

const GROUP = get(ENV, "GROUP", "All")

function activate_extensions_env()
    Pkg.activate("extensions")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
end

@time begin
    if GROUP == "All" || GROUP == "InterfaceI"
        @safetestset "Linear Algebra Test" include("linalg.jl")
        @safetestset "AbstractSystem Test" include("abstractsystem.jl")
        @safetestset "Variable Scope Tests" include("variable_scope.jl")
        @safetestset "Symbolic Parameters Test" include("symbolic_parameters.jl")
        @safetestset "Parsing Test" include("variable_parsing.jl")
        @safetestset "Simplify Test" include("simplify.jl")
        @safetestset "Direct Usage Test" include("direct.jl")
        @safetestset "System Linearity Test" include("linearity.jl")
        @safetestset "Linearization Tests" include("linearize.jl")
        @safetestset "Input Output Test" include("input_output_handling.jl")
        @safetestset "Clock Test" include("clock.jl")
        @safetestset "DiscreteSystem Test" include("discretesystem.jl")
        @safetestset "ODESystem Test" include("odesystem.jl")
        @safetestset "Unitful Quantities Test" include("units.jl")
        @safetestset "LabelledArrays Test" include("labelledarrays.jl")
        @safetestset "Mass Matrix Test" include("mass_matrix.jl")
        @safetestset "SteadyStateSystem Test" include("steadystatesystems.jl")
        @safetestset "SDESystem Test" include("sdesystem.jl")
        @safetestset "NonlinearSystem Test" include("nonlinearsystem.jl")
        @safetestset "PDE Construction Test" include("pde.jl")
        @safetestset "JumpSystem Test" include("jumpsystem.jl")
        @safetestset "Constraints Test" include("constraints.jl")
        @safetestset "Reduction Test" include("reduction.jl")
        @safetestset "Split Parameters Test" include("split_parameters.jl")
        @safetestset "ODAEProblem Test" include("odaeproblem.jl")
        @safetestset "StaticArrays Test" include("static_arrays.jl")
        @safetestset "Components Test" include("components.jl")
        @safetestset "Model Parsing Test" include("model_parsing.jl")
        @safetestset "print_tree" include("print_tree.jl")
        @safetestset "Error Handling" include("error_handling.jl")
        @safetestset "StructuralTransformations" include("structural_transformation/runtests.jl")
        @safetestset "State Selection Test" include("state_selection.jl")
        @safetestset "Symbolic Event Test" include("symbolic_events.jl")
        @safetestset "Stream Connect Test" include("stream_connectors.jl")
        @safetestset "Domain Connect Test" include("domain_connectors.jl")
        @safetestset "Lowering Integration Test" include("lowering_solving.jl")
        @safetestset "Test Big System Usage" include("bigsystem.jl")
        @safetestset "Dependency Graph Test" include("dep_graphs.jl")
        @safetestset "Function Registration Test" include("function_registration.jl")
        @safetestset "Precompiled Modules Test" include("precompile_test.jl")
        @safetestset "Variable Utils Test" include("variable_utils.jl")
        @safetestset "Variable Metadata Test" include("test_variable_metadata.jl")
        @safetestset "DAE Jacobians Test" include("dae_jacobian.jl")
        @safetestset "Jacobian Sparsity" include("jacobiansparsity.jl")
        @safetestset "Modelingtoolkitize Test" include("modelingtoolkitize.jl")
        @safetestset "OptimizationSystem Test" include("optimizationsystem.jl")
        @safetestset "FuncAffect Test" include("funcaffect.jl")
        @safetestset "Constants Test" include("constants.jl")
        @safetestset "Inverse Models Test" include("inversemodel.jl")
    end

    if GROUP == "All" || GROUP == "InterfaceII"
        println("C compilation test requires gcc available in the path!")
        @safetestset "C Compilation Test" include("ccompile.jl")
        @testset "Distributed Test" include("distributed.jl")
        @testset "Serialization" include("serialization.jl")
    end

    if GROUP == "All" || GROUP == "RegressionI"
        @safetestset "Latexify recipes Test" include("latexify.jl")
    end

    if GROUP == "All" || GROUP == "Extensions"
        activate_extensions_env()
        @safetestset "BifurcationKit Extension Test" include("extensions/bifurcationkit.jl")
    end
end
