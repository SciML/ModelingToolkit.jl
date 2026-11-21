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
        solve(odeprob, FBDF(); abstol = 1.0e-6, reltol = 1.0e-6)
        solve(daeprob, FBDF(); abstol = 1.0e-6, reltol = 1.0e-6)
    end
end

end
