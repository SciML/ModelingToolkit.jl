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

@testset "Change independent variable (trivial)" begin
    @variables x(t) y(t)
    eqs1 = [D(D(x)) ~ D(x) + x, D(y) ~ 1]
    M1 = ODESystem(eqs1, t; name = :M) |> complete
    M2 = change_independent_variable(M1, M1.y)
    eqs2 = substitute(equations(M2), M2.y => M1.t) # system should be equivalent when parametrized with y (since D(y) ~ 1), so substitute back ...
    @test eqs1[1] == only(eqs2) # ... and check that the equations are unmodified
end

@testset "Change independent variable" begin
    @variables x(t) y(t) z(t) s(t)
    eqs = [
        D(x) ~ y
        D(D(y)) ~ 2 * x * D(y)
        z ~ x + D(y)
        D(s) ~ 1 / (2*s)
    ]
    initialization_eqs = [x ~ 1.0, y ~ 1.0, D(y) ~ 0.0]
    M1 = ODESystem(eqs, t; initialization_eqs, name = :M) |> complete
    M2 = change_independent_variable(M1, M1.s; dummies = true)

    M1 = structural_simplify(M1; allow_symbolic = true)
    M2 = structural_simplify(M2; allow_symbolic = true)
    prob1 = ODEProblem(M1, [M1.s => 1.0], (1.0, 4.0), [])
    prob2 = ODEProblem(M2, [], (1.0, 2.0), [])
    sol1 = solve(prob1, Tsit5(); reltol = 1e-10, abstol = 1e-10)
    sol2 = solve(prob2, Tsit5(); reltol = 1e-10, abstol = 1e-10)
    ts = range(0.0, 1.0, length = 50)
    ss = .√(ts)
    @test all(isapprox.(sol1(ts, idxs=M1.x), sol2(ss, idxs=M2.x); atol = 1e-7)) &&
          all(isapprox.(sol1(ts, idxs=M1.y), sol2(ss, idxs=M2.y); atol = 1e-7))
end

@testset "Change independent variable (Friedmann equation)" begin
    @independent_variables t
    D = Differential(t)
    @variables a(t) ȧ(t) Ω(t) ϕ(t)
    a, ȧ = GlobalScope.([a, ȧ])
    species(w; kw...) = ODESystem([D(Ω) ~ -3(1 + w) * D(a)/a * Ω], t, [Ω], []; kw...)
    @named r = species(1//3)
    @named m = species(0)
    @named Λ = species(-1)
    eqs = [
        Ω ~ r.Ω + m.Ω + Λ.Ω
        D(a) ~ ȧ
        ȧ ~ √(Ω) * a^2
        D(D(ϕ)) ~ -3*D(a)/a*D(ϕ)
    ]
    M1 = ODESystem(eqs, t, [Ω, a, ȧ, ϕ], []; name = :M)
    M1 = compose(M1, r, m, Λ)
    M1 = complete(M1; flatten = false)

    # Apply in two steps, where derivatives are defined at each step: first t -> a, then a -> b
    M2 = change_independent_variable(M1, M1.a; dummies = true)
    a, ȧ, Ω, Ωr, Ωm, ΩΛ, ϕ, a_t, a_tt = M2.a, M2.ȧ, M2.Ω, M2.r.Ω, M2.m.Ω, M2.Λ.Ω, M2.ϕ, M2.a_t, M2.a_tt
    Da = Differential(a)
    @test Set(equations(M2)) == Set([
        a_t ~ ȧ # 1st order dummy equation
        a_tt ~ Da(ȧ) * a_t # 2nd order dummy equation
        Ω ~ Ωr + Ωm + ΩΛ
        ȧ ~ √(Ω) * a^2
        a_tt*Da(ϕ) + a_t^2*(Da^2)(ϕ) ~ -3*a_t^2/a*Da(ϕ)
        a_t*Da(Ωr) ~ -4*Ωr*a_t/a
        a_t*Da(Ωm) ~ -3*Ωm*a_t/a
        a_t*Da(ΩΛ) ~ 0
    ])

    @variables b(M2.a)
    M3 = change_independent_variable(M2, b, [Differential(M2.a)(b) ~ exp(-b), M2.a ~ exp(b)])

    M1 = structural_simplify(M1)
    M2 = structural_simplify(M2; allow_symbolic = true)
    M3 = structural_simplify(M3; allow_symbolic = true)
    @test length(unknowns(M3)) == length(unknowns(M2)) == length(unknowns(M1)) - 1
end

@testset "Change independent variable (simple)" begin
    @variables x(t)
    Mt = ODESystem([D(x) ~ 2*x], t; name = :M) |> complete
    Mx = change_independent_variable(Mt, Mt.x; dummies = true)
    @test (@variables x x_t(x) x_tt(x); Set(equations(Mx)) == Set([x_t ~ 2*x, x_tt ~ 2*x_t]))
end

@testset "Change independent variable (free fall)" begin
    @variables x(t) y(t)
    @parameters g v # gravitational acceleration and constant horizontal velocity
    Mt = ODESystem([D(D(y)) ~ -g, D(x) ~ v], t; name = :M) |> complete # gives (x, y) as function of t, ...
    Mx = change_independent_variable(Mt, Mt.x; dummies = false) # ... but we want y as a function of x
    Mx = structural_simplify(Mx; allow_symbolic = true)
    Dx = Differential(Mx.x)
    prob = ODEProblem(Mx, [Mx.y => 0.0, Dx(Mx.y) => 1.0], (0.0, 20.0), [g => 9.81, v => 10.0]) # 1 = dy/dx = (dy/dt)/(dx/dt) means equal initial horizontal and vertical velocities
    sol = solve(prob, Tsit5(); reltol = 1e-5)
    @test all(isapprox.(sol[Mx.y], sol[Mx.x - g*(Mx.x/v)^2/2]; atol = 1e-10)) # compare to analytical solution (x(t) = v*t, y(t) = v*t - g*t^2/2)
end

@testset "Change independent variable (crazy analytical example)" begin
    @independent_variables t
    D = Differential(t)
    @variables x(t) y(t)
    M1 = ODESystem([ # crazy non-autonomous non-linear 2nd order ODE
        D(D(y)) ~ D(x)^2 + D(y^3) |> expand_derivatives # expand D(y^3) # TODO: make this test 3rd order
        D(x) ~ x^4 + y^5 + t^6
    ], t; name = :M) |> complete
    M2 = change_independent_variable(M1, M1.x; dummies = true)

    # Compare to pen-and-paper result
    @independent_variables x
    Dx = Differential(x)
    @variables x_t(x) x_tt(x) y(x) t(x)
    @test Set(equations(M2)) == Set([
        x_t^2*(Dx^2)(y) + x_tt*Dx(y) ~ x_t^2 + 3*y^2*Dx(y)*x_t # from D(D(y))
        x_t ~ x^4 + y^5 + t^6 # 1st order dummy equation
        x_tt ~ 4*x^3*x_t + 5*y^4*Dx(y)*x_t + 6*t^5 # 2nd order dummy equation
    ])
end

@testset "Change independent variable (errors)" begin
    @variables x(t) y z(y) w(t) v(t)
    M = ODESystem([D(x) ~ 0, v ~ x], t; name = :M)
    @test_throws "incomplete" change_independent_variable(M, M.x)
    M = complete(M)
    @test_throws "singular" change_independent_variable(M, M.x)
    @test_throws "structurally simplified" change_independent_variable(structural_simplify(M), y)
    @test_throws "not specified" change_independent_variable(M, w)
    @test_throws "not specified" change_independent_variable(M, v)
    @test_throws "not a function of the independent variable" change_independent_variable(M, y)
    @test_throws "not a function of the independent variable" change_independent_variable(M, z)
end
