using ModelingToolkit, OrdinaryDiffEq, Test

@independent_variables t
D = Differential(t)

@testset "Liouville transform" begin
    @parameters α β γ δ
    @variables x(t) y(t)
    eqs = [D(x) ~ α * x - β * x * y, D(y) ~ -δ * y + γ * x * y]
    @named sys = ODESystem(eqs, t)
    sys = complete(sys)

    u0 = [x => 1.0, y => 1.0]
    p = [α => 1.5, β => 1.0, δ => 3.0, γ => 1.0]
    tspan = (0.0, 10.0)
    prob = ODEProblem(sys, u0, tspan, p)
    sol = solve(prob, Tsit5())

    sys2 = liouville_transform(sys)
    sys2 = complete(sys2)

    u0 = [x => 1.0, y => 1.0, sys2.trJ => 1.0]
    prob = ODEProblem(sys2, u0, tspan, p, jac = true)
    sol = solve(prob, Tsit5())
    @test sol[end, end] ≈ 1.0742818931017244
end

@testset "Change independent variable" begin
    @variables x(t) y(t) z(t) s(t)
    eqs = [
        D(x) ~ y
        D(D(y)) ~ 2 * x * D(y)
        z ~ x + D(y)
        D(s) ~ 1 / (2*s)
    ]
    M1 = ODESystem(eqs, t; name = :M) |> complete
    M2 = ModelingToolkit.change_independent_variable(M1, M1.s)

    M1 = structural_simplify(M1; allow_symbolic = true)
    M2 = structural_simplify(M2; allow_symbolic = true)
    prob1 = ODEProblem(M1, [M1.x => 1.0, M1.y => 1.0, Differential(M1.t)(M1.y) => 0.0, M1.s => 1.0], (1.0, 4.0))
    prob2 = ODEProblem(M2, [M2.x => 1.0, M2.y => 1.0, Differential(M2.s)(M2.y) => 0.0], (1.0, 2.0))
    sol1 = solve(prob1, Tsit5(); reltol = 1e-10, abstol = 1e-10)
    sol2 = solve(prob2, Tsit5(); reltol = 1e-10, abstol = 1e-10)
    ts = range(0.0, 1.0, length = 50)
    ss = .√(ts)
    @test all(isapprox.(sol1(ts, idxs=M1.x), sol2(ss, idxs=M2.x); atol = 1e-7)) &&
          all(isapprox.(sol1(ts, idxs=M1.y), sol2(ss, idxs=M2.y); atol = 1e-7))
end

@testset "Change independent variable (Friedmann equation)" begin
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

    M2 = ModelingToolkit.change_independent_variable(M1, M1.a)
    M2 = structural_simplify(M2; allow_symbolic = true)
    @test length(unknowns(M2)) == 2
end
