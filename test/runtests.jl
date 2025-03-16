using SafeTestsets, Pkg, Test
# https://github.com/JuliaLang/julia/issues/54664
import REPL

const GROUP = get(ENV, "GROUP", "All")

function activate_fmi_env()
    Pkg.activate("fmi")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
end

function activate_extensions_env()
    Pkg.activate("extensions")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
end

function activate_downstream_env()
    Pkg.activate("downstream")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
end

@time begin
    if GROUP == "All" || GROUP == "InterfaceI"
        @testset "InterfaceI" begin
            @safetestset "Linear Algebra Test" include("linalg.jl")
            @safetestset "AbstractSystem Test" include("abstractsystem.jl")
            @safetestset "Variable Scope Tests" include("variable_scope.jl")
            @safetestset "Symbolic Parameters Test" include("symbolic_parameters.jl")
            @safetestset "Parsing Test" include("variable_parsing.jl")
            @safetestset "Simplify Test" include("simplify.jl")
            @safetestset "Direct Usage Test" include("direct.jl")
            @safetestset "System Linearity Test" include("linearity.jl")
            @safetestset "Input Output Test" include("input_output_handling.jl")
            @safetestset "Clock Test" include("clock.jl")
            @safetestset "ODESystem Test" include("odesystem.jl")
            @safetestset "Dynamic Quantities Test" include("dq_units.jl")
            @safetestset "Unitful Quantities Test" include("units.jl")
            @safetestset "Mass Matrix Test" include("mass_matrix.jl")
            @safetestset "Reduction Test" include("reduction.jl")
            @safetestset "Split Parameters Test" include("split_parameters.jl")
            @safetestset "StaticArrays Test" include("static_arrays.jl")
            @safetestset "Components Test" include("components.jl")
            @safetestset "Model Parsing Test" include("model_parsing.jl")
            @safetestset "Error Handling" include("error_handling.jl")
            @safetestset "StructuralTransformations" include("structural_transformation/runtests.jl")
            @safetestset "State Selection Test" include("state_selection.jl")
            @safetestset "Symbolic Event Test" include("symbolic_events.jl")
            @safetestset "Stream Connect Test" include("stream_connectors.jl")
            @safetestset "Domain Connect Test" include("domain_connectors.jl")
            @safetestset "Lowering Integration Test" include("lowering_solving.jl")
            @safetestset "Dependency Graph Test" include("dep_graphs.jl")
            @safetestset "Function Registration Test" include("function_registration.jl")
            @safetestset "Precompiled Modules Test" include("precompile_test.jl")
            @safetestset "DAE Jacobians Test" include("dae_jacobian.jl")
            @safetestset "Jacobian Sparsity" include("jacobiansparsity.jl")
            @safetestset "Modelingtoolkitize Test" include("modelingtoolkitize.jl")
            @safetestset "FuncAffect Test" include("funcaffect.jl")
            @safetestset "Constants Test" include("constants.jl")
            @safetestset "Parameter Dependency Test" include("parameter_dependencies.jl")
            @safetestset "Equation Type Accessors Test" include("equation_type_accessors.jl")
            @safetestset "System Accessor Functions Test" include("accessor_functions.jl")
            @safetestset "Equations with complex values" include("complex.jl")
        end
    end

    if GROUP == "All" || GROUP == "Initialization"
        @safetestset "Guess Propagation" include("guess_propagation.jl")
        @safetestset "Hierarchical Initialization Equations" include("hierarchical_initialization_eqs.jl")
        @safetestset "InitializationSystem Test" include("initializationsystem.jl")
        @safetestset "Initial Values Test" include("initial_values.jl")
    end

    if GROUP == "All" || GROUP == "InterfaceII"
        @testset "InterfaceII" begin
            @safetestset "Code Generation Test" include("code_generation.jl")
            @safetestset "IndexCache Test" include("index_cache.jl")
            @safetestset "Variable Utils Test" include("variable_utils.jl")
            @safetestset "Variable Metadata Test" include("test_variable_metadata.jl")
            @safetestset "OptimizationSystem Test" include("optimizationsystem.jl")
            @safetestset "Discrete System" include("discrete_system.jl")
            @safetestset "Implicit Discrete System" include("implicit_discrete_system.jl")
            @safetestset "SteadyStateSystem Test" include("steadystatesystems.jl")
            @safetestset "SDESystem Test" include("sdesystem.jl")
            @safetestset "DDESystem Test" include("dde.jl")
            @safetestset "NonlinearSystem Test" include("nonlinearsystem.jl")
            @safetestset "SCCNonlinearProblem Test" include("scc_nonlinear_problem.jl")
            @safetestset "PDE Construction Test" include("pdesystem.jl")
            @safetestset "JumpSystem Test" include("jumpsystem.jl")
            @safetestset "BVProblem Test" include("bvproblem.jl")
            @safetestset "print_tree" include("print_tree.jl")
            @safetestset "Constraints Test" include("constraints.jl")
            @safetestset "IfLifting Test" include("if_lifting.jl")
            @safetestset "Analysis Points Test" include("analysis_points.jl")
            @safetestset "Causal Variables Connection Test" include("causal_variables_connection.jl")
            @safetestset "Debugging Test" include("debugging.jl")
        end
    end

    if GROUP == "All" || GROUP == "SymbolicIndexingInterface"
        @safetestset "SymbolicIndexingInterface test" include("symbolic_indexing_interface.jl")
        @safetestset "SciML Problem Input Test" include("sciml_problem_inputs.jl")
        @safetestset "MTKParameters Test" include("mtkparameters.jl")
    end

    if GROUP == "All" || GROUP == "Extended"
        @safetestset "Test Big System Usage" include("bigsystem.jl")
        println("C compilation test requires gcc available in the path!")
        @safetestset "C Compilation Test" include("ccompile.jl")
        @testset "Distributed Test" include("distributed.jl")
        @testset "Serialization" include("serialization.jl")
    end

    if GROUP == "All" || GROUP == "RegressionI"
        @safetestset "Latexify recipes Test" include("latexify.jl")
    end

    if GROUP == "All" || GROUP == "Downstream"
        activate_downstream_env()
        @safetestset "Linearization Tests" include("downstream/linearize.jl")
        @safetestset "Linearization Dummy Derivative Tests" include("downstream/linearization_dd.jl")
        @safetestset "Inverse Models Test" include("downstream/inversemodel.jl")
        @safetestset "Analysis Points Test" include("downstream/analysis_points.jl")
        @safetestset "Analysis Points Test" include("downstream/test_disturbance_model.jl")
    end

    if GROUP == "All" || GROUP == "FMI"
        activate_fmi_env()
        @safetestset "FMI Extension Test" include("fmi/fmi.jl")
    end

    if GROUP == "All" || GROUP == "Extensions"
        activate_extensions_env()
        @safetestset "HomotopyContinuation Extension Test" include("extensions/homotopy_continuation.jl")
        @safetestset "Auto Differentiation Test" include("extensions/ad.jl")
        @safetestset "LabelledArrays Test" include("labelledarrays.jl")
        @safetestset "BifurcationKit Extension Test" include("extensions/bifurcationkit.jl")
        @safetestset "InfiniteOpt Extension Test" include("extensions/test_infiniteopt.jl")
    end
end
