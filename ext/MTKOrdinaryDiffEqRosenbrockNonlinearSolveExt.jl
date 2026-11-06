module MTKOrdinaryDiffEqRosenbrockNonlinearSolveExt

using ModelingToolkit
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqNonlinearSolve: OrdinaryDiffEqNonlinearSolve
using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    daeprob = ModelingToolkit.precompile_dae_problem()
    sccprob = ModelingToolkit.precompile_scc_dae_problem()
    @compile_workload begin
        solve(daeprob, Rodas5P())
        solve(sccprob, Rodas5P())
    end
end

end
