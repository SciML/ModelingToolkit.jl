module MTKOrdinaryDiffEqTsit5Ext

using ModelingToolkit
using OrdinaryDiffEqTsit5: Tsit5
using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    odeprob = ModelingToolkit.precompile_ode_problem()
    @compile_workload begin
        solve(odeprob, Tsit5())
    end
end

end
