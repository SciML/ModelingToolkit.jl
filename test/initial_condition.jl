using ModelingToolkit
using OrdinaryDiffEq

# This example adds a default initial condition for x, and makes sure x is chosen as the state variables
# The user then supplies an initial condition for y which is in conflict with the default initial condition for x
@variables t x(t)=1 y(t)
D = Differential(t)
eqs = [D(x) ~ -x
    y ~ x]

@named sys = ODESystem(eqs, t)
ssys = structural_simplify(sys)

if VERSION >= v"1.8"
    @test_throws "The user-provided initial condition for y(t) = 2.0 is conflicting with another initial condition (1.0) that takes precedent." begin
        prob = ODEProblem(ssys, [y => 2.0], (0.0, 1.0))
        sol = solve(prob, Tsit5())
        @test sol(0.0, idxs = y) == 2.0
    end
end
