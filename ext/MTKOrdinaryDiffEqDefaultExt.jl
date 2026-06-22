module MTKOrdinaryDiffEqDefaultExt

using ModelingToolkit
using OrdinaryDiffEqDefault
using PrecompileTools
using ModelingToolkit: t_nounits, D_nounits

@setup_workload begin
    @parameters a = 1.0 b = 1.0
    @variables x(t_nounits) y(t_nounits)
    prob = ODEProblem(
        mtkcompile(
            System(
                [D_nounits(x) ~ a * y, D_nounits(y) ~ -b * x],
                t_nounits; name = :precompile_default
            )
        ),
        [x => 1.0, y => 0.0], (0.0, 1.0)
    )
    @compile_workload begin
        solve(prob)
    end
end

end
