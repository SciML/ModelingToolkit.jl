module MTKOrdinaryDiffEqDefaultExt

using ModelingToolkit
using OrdinaryDiffEqDefault: OrdinaryDiffEqDefault
using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    odeprob = ModelingToolkit.precompile_ode_problem()
    daeprob = ModelingToolkit.precompile_dae_problem()
    @compile_workload begin
        solve(odeprob)
        solve(daeprob)
    end
end

end
