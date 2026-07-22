using SafeTestsets, Pkg, Test
# https://github.com/JuliaLang/julia/issues/54664
import REPL

# The centralized SciML sublibrary CI (SublibraryCI.yml -> sublibrary-tests.yml@v1)
# emits GROUP="ModelingToolkitBase" for the [Core] section and
# GROUP="ModelingToolkitBase_<Section>" for every other section in test_groups.toml.
# Strip that "<pkg>_" prefix so the bare section names below ("InterfaceI", "QA", …)
# drive the dispatch. A bare "<pkg>" maps to "Core"; anything else (e.g. "All" for
# local runs) is passed through unchanged.
const _G = get(ENV, "GROUP", "All")
const _SUB = "ModelingToolkitBase"
const GROUP = _G == _SUB ? "Core" :
    (startswith(_G, _SUB * "_") ? _G[(length(_SUB) + 2):end] : _G)

function activate_extensions_env()
    Pkg.activate("extensions")
    Pkg.develop([PackageSpec(path = dirname(@__DIR__))])
    return Pkg.instantiate()
end

function activate_optimization_env()
    Pkg.activate("optimization")
    Pkg.develop([PackageSpec(path = dirname(@__DIR__))])
    return Pkg.instantiate()
end

function activate_qa_env()
    Pkg.activate("qa")
    Pkg.develop([PackageSpec(path = dirname(@__DIR__))])
    return Pkg.instantiate()
end

function activate_downstream_env()
    Pkg.activate("downstream")
    Pkg.develop([PackageSpec(path = dirname(@__DIR__))])
    return Pkg.instantiate()
end

@time begin
    if GROUP == "All" || GROUP == "InterfaceI"
        @testset "InterfaceI" begin
            @safetestset "AbstractSystem Test" include("abstractsystem.jl")
            @safetestset "Variable Scope Tests" include("variable_scope.jl")
            @safetestset "Parsing Test" include("variable_parsing.jl")
            @safetestset "System Linearity Test" include("linearity.jl")
            @safetestset "Variable binding semantics" include("binding_semantics.jl")
            @safetestset "Input Output Test" include("input_output_handling.jl")
            @safetestset "System Building Tests" include("system_building.jl")
            @safetestset "Simple `mtkcompile`" include("simple_mtkcompile.jl")
            @safetestset "`@mtkcomplete`" include("mtkcomplete.jl")
            @safetestset "ODESystem Test" include("odesystem.jl")
            @safetestset "Dynamic Quantities Test" include("dq_units.jl")
            @safetestset "Mass Matrix Test" include("mass_matrix.jl")
            @safetestset "Split Parameters Test" include("split_parameters.jl")
            @safetestset "StaticArrays Test" include("static_arrays.jl")
            @safetestset "Components Test" include("components.jl")
            @safetestset "Error Handling" include("error_handling.jl")
            @safetestset "Basic transformations" include("basic_transformations.jl")
            @safetestset "Change of variables" include("changeofvariables.jl")
            @safetestset "Symbolic Event Test" include("symbolic_events.jl")
            @safetestset "Stream Connect Test" include("stream_connectors.jl")
            @safetestset "Domain Connect Test" include("domain_connectors.jl")
            @safetestset "Dependency Graph Test" include("dep_graphs.jl")
            @safetestset "Function Registration Test" include("function_registration.jl")
            @safetestset "Precompiled Modules Test" include("precompile_test.jl")
            @safetestset "DAE Jacobians Test" include("dae_jacobian.jl")
            @safetestset "Jacobian Sparsity" include("jacobiansparsity.jl")
            @safetestset "Modelingtoolkitize Test" include("modelingtoolkitize.jl")
            @safetestset "Constants Test" include("constants.jl")
            @safetestset "Parameter Bindings Test" include("parameter_bindings.jl")
            @safetestset "Equation Type Accessors Test" include("equation_type_accessors.jl")
            @safetestset "System Accessor Functions Test" include("accessor_functions.jl")
            @safetestset "Equations with complex values" include("complex.jl")
        end
    end

    if GROUP == "All" || GROUP == "Initialization"
        @safetestset "Guess Propagation" include("guess_propagation.jl")
        @safetestset "InitializationSystem Test" include("initializationsystem.jl")
        @safetestset "Initial Values Test" include("initial_values.jl")
    end

    if GROUP == "All" || GROUP == "InterfaceII"
        @testset "InterfaceII" begin
            @safetestset "Code Generation Test" include("code_generation.jl")
            @safetestset "IndexCache Test" include("index_cache.jl")
            @safetestset "Variable Utils Test" include("variable_utils.jl")
            @safetestset "Variable Metadata Test" include("test_variable_metadata.jl")
            @safetestset "Discrete System" include("discrete_system.jl")
            @safetestset "Implicit Discrete System" include("implicit_discrete_system.jl")
            @safetestset "SteadyStateSystem Test" include("steadystatesystems.jl")
            @safetestset "SDESystem Test" include("sdesystem.jl")
            @safetestset "DDESystem Test" include("dde.jl")
            @safetestset "NonlinearSystem Test" include("nonlinearsystem.jl")
            @safetestset "Homotopy lowering" include("homotopy_lowering.jl")
            @safetestset "Homotopy problem construction & sweep" include("homotopy_problem.jl")
            @safetestset "Homotopy OMC parity" include("homotopy_omc_parity.jl")
            @safetestset "Homotopy initialization routing" include("homotopy_initialization.jl")
            @safetestset "PDE Construction Test" include("pdesystem.jl")
            @safetestset "JumpSystem Test" include("jumpsystem.jl")
            @safetestset "Poissonians Test" include("poissonians.jl")
            @safetestset "Extend SDE/Jump Test" include("extend_sde_jump.jl")
            @safetestset "print_tree" include("print_tree.jl")
            @safetestset "Analysis Points Test" include("analysis_points.jl")
            @safetestset "Causal Variables Connection Test" include("causal_variables_connection.jl")
            @safetestset "Debugging Test" include("debugging.jl")
            @safetestset "Namespacing test" include("namespacing.jl")
            @safetestset "LinearProblem Tests" include("linearproblem.jl")
            @safetestset "Optimal Control + Constraints Tests" include("bvproblem.jl")
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
        @safetestset "Multithreading test" include("multithreading.jl")
        @testset "Distributed Test" include("distributed.jl")
        @testset "Serialization" include("serialization.jl")
    end

    if GROUP == "All" || GROUP == "RegressionI"
        @safetestset "Latexify recipes Test" include("latexify.jl")
    end

    if GROUP == "All" || GROUP == "Extensions"
        activate_extensions_env()
        @safetestset "HomotopyContinuation Extension Test" include("extensions/homotopy_continuation.jl")
        @safetestset "LabelledArrays Test" include("extensions/labelledarrays.jl")
        @safetestset "BifurcationKit Extension Test" include("extensions/bifurcationkit.jl")
        # @safetestset "Auto Differentiation Test" include("extensions/ad.jl")
    end

    if GROUP == "All" || GROUP == "Optimization"
        activate_optimization_env()
        @safetestset "OptimizationSystem Test" include("optimization/optimizationsystem.jl")
        @safetestset "InfiniteOpt Extension Test" include("optimization/test_infiniteopt.jl")
        @safetestset "Dynamic Optimization Collocation Solvers" include("optimization/dynamic_optimization.jl")
    end

    if GROUP == "All" || GROUP == "QA"
        activate_qa_env()
        @safetestset "JET Tests" include("qa/jet.jl")
        @safetestset "Aqua Tests" include("qa/aqua.jl")
    end
end
