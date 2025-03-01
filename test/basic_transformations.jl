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
    @independent_variables t
    @variables x(t) y(t) z(t) s(t)
    eqs = [
        D(x) ~ y
        D(D(y)) ~ 2 * x * D(y)
        z ~ x + D(y)
        D(s) ~ 1 / (2*s)
    ]
    @named sys1 = ODESystem(eqs, t)
    sys1 = complete(sys1)

    @independent_variables s
    sys2 = ModelingToolkit.change_independent_variable(sys1, s)

    sys1 = structural_simplify(sys1; allow_symbolic = true)
    sys2 = structural_simplify(sys2; allow_symbolic = true)
    prob1 = ODEProblem(sys1, [sys1.x => 1.0, sys1.y => 1.0, Differential(t)(sys1.y) => 0.0, sys1.s => 1.0], (1.0, 4.0))
    prob2 = ODEProblem(sys2, [sys2.x => 1.0, sys2.y => 1.0, Differential(s)(sys2.y) => 0.0], (1.0, 2.0))
    sol1 = solve(prob1, Tsit5(); reltol = 1e-10, abstol = 1e-10)
    sol2 = solve(prob2, Tsit5(); reltol = 1e-10, abstol = 1e-10)
    ts = range(0.0, 1.0, length = 50)
    ss = .√(ts)
    @test all(isapprox.(sol1(ts, idxs=sys1.x), sol2(ss, idxs=sys2.x); atol = 1e-7)) &&
          all(isapprox.(sol1(ts, idxs=sys1.y), sol2(ss, idxs=sys2.y); atol = 1e-7))
end

@testset "Change independent variable (Friedmann equation)" begin
    @independent_variables t
    D = Differential(t)
    @variables a(t) ρr(t) ρm(t) ρΛ(t) ρ(t) P(t) ϕ(t)
    @parameters Ωr0 Ωm0 ΩΛ0
    eqs = [
        ρr ~ 3/(8*Num(π)) * Ωr0 / a^4
        ρm ~ 3/(8*Num(π)) * Ωm0 / a^3
        ρΛ ~ 3/(8*Num(π)) * ΩΛ0
        ρ ~ ρr + ρm + ρΛ
        D(a) ~ √(8*Num(π)/3*ρ*a^4)
        D(D(ϕ)) ~ -3*D(a)/a*D(ϕ)
    ]
    @named M1 = ODESystem(eqs, t)
    M1 = complete(M1)

    @independent_variables a
    M2 = ModelingToolkit.change_independent_variable(M1, a)
    M2 = structural_simplify(M2; allow_symbolic = true)
    @test length(unknowns(M2)) == 2
end
