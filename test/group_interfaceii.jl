include("shared/mtktestset.jl")

@testset "InterfaceII" begin
    @mtktestset("Code Generation Test", "code_generation.jl")
    @mtktestset("Discrete System", "discrete_system.jl")
    @mtktestset("Implicit Discrete System", "implicit_discrete_system.jl")
    @mtktestset("SDESystem Test", "sdesystem.jl")
    @mtktestset("DDESystem Test", "dde.jl")
    @mtktestset("NonlinearSystem Test", "nonlinearsystem.jl")
    @safetestset "SCCNonlinearProblem Test" include("scc_nonlinear_problem.jl")
    @safetestset "IfLifting Test" include("if_lifting.jl")
    @safetestset "Simplification determinism" include("nondeterminism_simplification.jl")
    @mtktestset("Analysis Points Test", "analysis_points.jl")
    @mtktestset("Causal Variables Connection Test", "causal_variables_connection.jl")
    @safetestset "Subsystem replacement" include("substitute_component.jl")
    @safetestset "Linearization Tests" include("linearize.jl")
    @safetestset "Fractional Differential Equations Tests" include("fractional_to_ordinary.jl")
    @safetestset "SemilinearODEProblem tests" include("semilinearodeproblem.jl")
end
