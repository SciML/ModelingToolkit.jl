using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Zygote
using SymbolicIndexingInterface
using SciMLStructures
using OrdinaryDiffEq
using SciMLSensitivity
using ForwardDiff

@variables x(t)[1:3] y(t)
@parameters p[1:3, 1:3] q
eqs = [D(x) ~ p * x
       D(y) ~ sum(p) + q * y]
u0 = [x => zeros(3),
    y => 1.0]
ps = [p => zeros(3, 3),
    q => 1.0]
tspan = (0.0, 10.0)
@mtkbuild sys = ODESystem(eqs, t)
prob = ODEProblem(sys, u0, tspan, ps)
sol = solve(prob, Tsit5())

mtkparams = parameter_values(prob)
new_p = rand(10)
gs = gradient(new_p) do new_p
    new_params = SciMLStructures.replace(SciMLStructures.Tunable(), mtkparams, new_p)
    new_prob = remake(prob, p = new_params)
    new_sol = solve(new_prob, Tsit5())
    sum(new_sol)
end

@testset "Issue#2997" begin
    pars = @parameters y0 mh Tγ0 Th0 h ργ0
    vars = @variables x(t)
    @named sys = ODESystem([D(x) ~ y0],
        t,
        vars,
        pars;
        defaults = [
            y0 => mh * 3.1 / (2.3 * Th0),
            mh => 123.4,
            Th0 => (4 / 11)^(1 / 3) * Tγ0,
            Tγ0 => (15 / π^2 * ργ0 * (2 * h)^2 / 7)^(1 / 4) / 5
        ])
    sys = structural_simplify(sys)

    function x_at_0(θ)
        prob = ODEProblem(sys, [sys.x => 1.0], (0.0, 1.0), [sys.ργ0 => θ[1], sys.h => θ[2]])
        return prob.u0[1]
    end

    @test ForwardDiff.gradient(x_at_0, [0.3, 0.7]) == zeros(2)
end
