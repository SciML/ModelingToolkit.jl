using ForwardDiff
using ModelingToolkitBase
using ModelingToolkitBase: SymbolicDiscreteCallback, t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using SciMLBase
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

@testset "discrete parameter reconstruction preserves dual values" begin
    @parameters c = 2.0
    @discretes d(t) = 1.0
    @variables x(t)
    ev = SymbolicDiscreteCallback(
        0.5 => [d ~ Pre(d) + 1.0]; discrete_parameters = [d]
    )
    sys = mtkcompile(
        System(
            [D(x) ~ -c * x * d], t, [x], [c, d];
            initial_conditions = [x => 1.0], discrete_events = [ev]
        )
    )

    for specialize in (
            SciMLBase.AutoDespecialize, SciMLBase.AutoSpecialize, SciMLBase.FullSpecialize
        )
        prob = ODEProblem{true, specialize}(sys, [], (0.0, 1.0))

        terminal_value_p = c_value -> solve(
            remake(prob; p = [c => c_value]), Tsit5();
            abstol = 1.0e-12, reltol = 1.0e-12
        ).u[end][1]
        terminal_value_u0 = x0 -> solve(
            remake(prob; u0 = [x => x0]), Tsit5();
            abstol = 1.0e-12, reltol = 1.0e-12
        ).u[end][1]

        @test ForwardDiff.derivative(terminal_value_p, 2.2) ≈
            -1.5 * exp(-1.5 * 2.2) rtol = 1.0e-6
        @test ForwardDiff.derivative(terminal_value_u0, 1.0) ≈
            exp(-1.5 * 2.0) rtol = 1.0e-6
    end
end
