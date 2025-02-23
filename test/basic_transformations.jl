using ModelingToolkit, OrdinaryDiffEq, Test

@independent_variables t
@parameters α β γ δ
@variables x(t) y(t)
D = Differential(t)

eqs = [D(x) ~ α * x - β * x * y,
    D(y) ~ -δ * y + γ * x * y]

sys = ODESystem(eqs)

u0 = [x => 1.0,
    y => 1.0]

p = [α => 1.5,
    β => 1.0,
    δ => 3.0,
    γ => 1.0]

tspan = (0.0, 10.0)
prob = ODEProblem(sys, u0, tspan, p)
sol = solve(prob, Tsit5())

sys2 = liouville_transform(sys)
@variables trJ

u0 = [x => 1.0,
    y => 1.0,
    trJ => 1.0]

prob = ODEProblem(sys2, u0, tspan, p, jac = true)
sol = solve(prob, Tsit5())
@test sol[end, end] ≈ 1.0742818931017244

@testset "Change independent variable" begin
    # Autonomous 1st order (t doesn't appear in the equations)
    @independent_variables t
    @variables x(t) y(t) z(t)
    eqs = [
        D(x) ~ y
        D(y) ~ 2*D(x)
        z ~ x + D(y)
    ]
    @named sys1 = ODESystem(eqs, t)

    @independent_variables s
    @named sys2 = ModelingToolkit.change_independent_variable(sys1, s, s^2)

    sys1 = structural_simplify(sys1)
    sys2 = structural_simplify(sys2)
    prob1 = ODEProblem(sys1, unknowns(sys1) .=> 1.0, (0.0, 1.0))
    prob2 = ODEProblem(sys2, unknowns(sys2) .=> 1.0, (0.0, 1.0))
    sol1 = solve(prob1, Tsit5(); reltol = 1e-10, abstol = 1e-10)
    sol2 = solve(prob2, Tsit5(); reltol = 1e-10, abstol = 1e-10)
    ts = range(0.0, 1.0, length = 50)
    ss = .√(ts)
    @test isapprox(sol1(ts), sol2(ss); atol = 1e-8)
end
