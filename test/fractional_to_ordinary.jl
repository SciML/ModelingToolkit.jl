using ModelingToolkit, OrdinaryDiffEq, ODEInterfaceDiffEq, SpecialFunctions
using Test

# Testing for α < 1
# Uses example 1 from Section 7 of https://arxiv.org/pdf/2506.04188
@independent_variables t
@variables x(t)
D = Differential(t)
tspan = (0., 1.)

function expect(t, α)
    return (3/2*t^(α/2) - t^4)^2
end

α = 0.5
eqs = (9*gamma(1 + α)/4) - (3*t^(4 - α/2)*gamma(5 + α/2)/gamma(5 - α/2))
eqs += (gamma(9)*t^(8 - α)/gamma(9 - α)) + (3/2*t^(α/2)-t^4)^3 - x^(3/2)
sys = fractional_to_ordinary(eqs, x, α, 10^-7, 1)

prob = ODEProblem(sys, [], tspan)
sol = solve(prob, radau5(), abstol = 1e-10, reltol = 1e-10)

time = 0
while(time <= 1)
    @test isapprox(expect(time, α), sol(time, idxs=x), atol=1e-3)
    time += 0.1
end

α = 0.3
eqs = (9*gamma(1 + α)/4) - (3*t^(4 - α/2)*gamma(5 + α/2)/gamma(5 - α/2))
eqs += (gamma(9)*t^(8 - α)/gamma(9 - α)) + (3/2*t^(α/2)-t^4)^3 - x^(3/2)
sys = fractional_to_ordinary(eqs, x, α, 10^-7, 1)

prob = ODEProblem(sys, [], tspan)
sol = solve(prob, radau5(), abstol = 1e-10, reltol = 1e-10)

time = 0
while(time <= 1)
    @test isapprox(expect(time, α), sol(time, idxs=x), atol=1e-3)
    time += 0.1
end

α = 0.9
eqs = (9*gamma(1 + α)/4) - (3*t^(4 - α/2)*gamma(5 + α/2)/gamma(5 - α/2))
eqs += (gamma(9)*t^(8 - α)/gamma(9 - α)) + (3/2*t^(α/2)-t^4)^3 - x^(3/2)
sys = fractional_to_ordinary(eqs, x, α, 10^-7, 1)

prob = ODEProblem(sys, [], tspan)
sol = solve(prob, radau5(), abstol = 1e-10, reltol = 1e-10)

time = 0
while(time <= 1)
    @test isapprox(expect(time, α), sol(time, idxs=x), atol=1e-3)
    time += 0.1
end

# Testing for example 2 of Section 7
@independent_variables t
@variables x(t) y(t)
D = Differential(t)
tspan = (0., 220.)

sys = fractional_to_ordinary([1 - 4*x + x^2 * y, 3*x - x^2 * y], [x, y], [1.3, 0.8], 10^-6, 220; initials=[[1.2, 1], 2.8])
prob = ODEProblem(sys, [], tspan)
sol = solve(prob, radau5(), abstol = 1e-8, reltol = 1e-8)

@test isapprox(1.0097684171, sol(220, idxs=x), atol=1e-3)
@test isapprox(2.1581264031, sol(220, idxs=y), atol=1e-3)