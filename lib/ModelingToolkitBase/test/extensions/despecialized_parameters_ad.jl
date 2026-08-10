using ForwardDiff
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using SciMLSensitivity
using SciMLStructures
using Test
using Zygote

@testset "despecialized parameter sensitivities" begin
    @parameters a = 2.0
    @variables x(t) = 1.0
    sys = mtkcompile(System([D(x) ~ -a * x], t; name = :despecialized_parameter_ad))
    prob = ODEProblem(sys, [], (0.0, 1.0))

    function terminal_value(a_value)
        p = SciMLStructures.replace(SciMLStructures.Tunable(), prob.p, [a_value])
        remade = remake(prob; p)
        sol = solve(
            remade, Tsit5(); saveat = [1.0], abstol = 1.0e-10, reltol = 1.0e-10
        )
        return sol[x][end]
    end

    expected = -exp(-2)
    @test ForwardDiff.derivative(terminal_value, 2.0) ≈ expected rtol = 1.0e-6
    @test only(Zygote.gradient(terminal_value, 2.0)) ≈ expected rtol = 1.0e-6
end
