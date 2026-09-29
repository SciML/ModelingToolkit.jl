module MTKOrdinaryDiffEqRosenbrockExt

using ModelingToolkit
using OrdinaryDiffEqRosenbrock: Rodas5P
using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    odeprob = ModelingToolkit.precompile_ode_problem()
    @compile_workload begin
        solve(odeprob, Rodas5P())
    end
end

end
