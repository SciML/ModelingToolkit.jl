module MTKOrdinaryDiffEqRosenbrockNonlinearSolveExt

using ModelingToolkit
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqNonlinearSolve: OrdinaryDiffEqNonlinearSolve
using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    daeprob = ModelingToolkit.precompile_dae_problem()
    @compile_workload begin
        solve(daeprob, Rodas5P())
    end
end

end
