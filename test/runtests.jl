using SafeTestsets, Pkg, Test
# https://github.com/JuliaLang/julia/issues/54664
import REPL

const MTKBasePath = joinpath(dirname(@__DIR__), "lib", "ModelingToolkitBase")
const MTKBasePkgSpec = PackageSpec(; path = MTKBasePath)
Pkg.develop([MTKBasePkgSpec])

const GROUP = get(ENV, "GROUP", "All")

function activate_fmi_env()
    Pkg.activate("fmi")
    Pkg.develop([MTKBasePkgSpec, PackageSpec(path = dirname(@__DIR__))])
    return Pkg.instantiate()
end

function activate_extensions_env()
    Pkg.activate(joinpath(MTKBasePath, "test", "extensions"))
    Pkg.develop([MTKBasePkgSpec, PackageSpec(path = dirname(@__DIR__))])
    return Pkg.instantiate()
end

function activate_downstream_env()
    Pkg.activate("downstream")
    Pkg.develop([MTKBasePkgSpec, PackageSpec(path = dirname(@__DIR__))])
    return Pkg.instantiate()
end

macro mtktestset(name, file)
    return quote
        @safetestset $name begin
            using ModelingToolkit
            import ModelingToolkitBase
            include(joinpath(pkgdir(ModelingToolkitBase), "test", $file))
        end
    end
end

@time begin
    if GROUP == "All" || GROUP == "InterfaceI"
        @testset "InterfaceI" begin
            @mtktestset("Input Output Test", "input_output_handling.jl")
            @safetestset "Clock Test" include("clock.jl")
            @mtktestset("Variable binding semantics", "binding_semantics.jl")
            @mtktestset("ODESystem Test", "odesystem.jl")
            @mtktestset("Dynamic Quantities Test", "dq_units.jl")
            @safetestset "Reduction Test" include("reduction.jl")
            @mtktestset("Split Parameters Test", "split_parameters.jl")
            @mtktestset("Components Test", "components.jl")
            @safetestset "StructuralTransformations" include("structural_transformation/runtests.jl")
            @mtktestset("Basic transformations", "basic_transformations.jl")
            @mtktestset("Change of variables", "changeofvariables.jl")
            @safetestset "State Selection Test" include("state_selection.jl")
            @mtktestset("Symbolic Event Test", "symbolic_events.jl")
            @mtktestset("Stream Connect Test", "stream_connectors.jl")
            @mtktestset("Jacobian Sparsity", "jacobiansparsity.jl")
            @mtktestset("Modelingtoolkitize Test", "modelingtoolkitize.jl")
            @mtktestset("Constants Test", "constants.jl")
            @mtktestset("System Accessor Functions Test", "accessor_functions.jl")
        end
    end

    if GROUP == "All" || GROUP == "Initialization"
        @mtktestset("Guess Propagation", "guess_propagation.jl")
        @safetestset "Hierarchical Initialization Equations" include("hierarchical_initialization_eqs.jl")
        @mtktestset("InitializationSystem Test", "initializationsystem.jl")
        @mtktestset("Initial Values Test", "initial_values.jl")
    end

    if GROUP == "All" || GROUP == "InterfaceII"
        @testset "InterfaceII" begin
            @mtktestset("Code Generation Test", "code_generation.jl")
            @mtktestset("Discrete System", "discrete_system.jl")
            @mtktestset("Implicit Discrete System", "implicit_discrete_system.jl")
            @mtktestset("SDESystem Test", "sdesystem.jl")
            @mtktestset("DDESystem Test", "dde.jl")
            @mtktestset("NonlinearSystem Test", "nonlinearsystem.jl")
            @safetestset "SCCNonlinearProblem Test" include("scc_nonlinear_problem.jl")
            @safetestset "IfLifting Test" include("if_lifting.jl")
            @mtktestset("Analysis Points Test", "analysis_points.jl")
            @mtktestset("Causal Variables Connection Test", "causal_variables_connection.jl")
            @safetestset "Subsystem replacement" include("substitute_component.jl")
            @safetestset "Linearization Tests" include("linearize.jl")
            @safetestset "Fractional Differential Equations Tests" include("fractional_to_ordinary.jl")
            @safetestset "SemilinearODEProblem tests" include("semilinearodeproblem.jl")
            @mtktestset("OptimizationSystem Test", "optimizationsystem.jl")
        end
    end

    if GROUP == "All" || GROUP == "SymbolicIndexingInterface"
        @mtktestset("SciML Problem Input Test", "sciml_problem_inputs.jl")
        @mtktestset("MTKParameters Test", "mtkparameters.jl")
    end

    if GROUP == "All" || GROUP == "Downstream"
        activate_downstream_env()
        @safetestset "Linearization Dummy Derivative Tests" include("downstream/linearization_dd.jl")
        @safetestset "Inverse Models Test" include("downstream/inversemodel.jl")
        @safetestset "Disturbance model Test" include("downstream/test_disturbance_model.jl")
    end

    if GROUP == "All" || GROUP == "FMI"
        activate_fmi_env()
        @safetestset "FMI Extension Test" include("fmi/fmi.jl")
    end

    if GROUP == "All" || GROUP == "Extensions"
        activate_extensions_env()
        @mtktestset("HomotopyContinuation Extension Test", "extensions/homotopy_continuation.jl")
        @mtktestset("BifurcationKit Extension Test", "extensions/bifurcationkit.jl")
        @mtktestset("InfiniteOpt Extension Test", "extensions/test_infiniteopt.jl")
        # @mtktestset("Auto Differentiation Test", "extensions/ad.jl")
        @mtktestset("Dynamic Optimization Collocation Solvers", "extensions/dynamic_optimization.jl")
    end
end
