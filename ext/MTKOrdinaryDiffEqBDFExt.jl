module MTKOrdinaryDiffEqBDFExt

using ModelingToolkit
using OrdinaryDiffEqBDF: FBDF
using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    odeprob = ModelingToolkit.precompile_ode_problem()
    daeprob = ModelingToolkit.precompile_dae_problem()
    @compile_workload begin
        solve(odeprob, FBDF())
        solve(daeprob, FBDF())
    end
end

end
