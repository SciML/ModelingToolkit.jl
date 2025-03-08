using ModelingToolkit, OrdinaryDiffEq, DataInterpolations, Test

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
    M1 = ODESystem(eqs1, t; name = :M)
    M2 = change_independent_variable(M1, y)
    @variables y x(y) yˍt(y)
    Dy = Differential(y)
    @test Set(equations(M2)) == Set([yˍt^2*(Dy^2)(x) + yˍt*Dy(yˍt)*Dy(x) ~ x + Dy(x)*yˍt, yˍt ~ 1])
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
    M1 = ODESystem(eqs, t; initialization_eqs, name = :M)
    M2 = change_independent_variable(M1, s)

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

    # Apply in two steps, where derivatives are defined at each step: first t -> a, then a -> b
    M2 = change_independent_variable(M1, M1.a)
    M2c = complete(M2) # just for the following equation comparison (without namespacing)
    a, ȧ, Ω, Ωr, Ωm, ΩΛ, ϕ, aˍt = M2c.a, M2c.ȧ, M2c.Ω, M2c.r.Ω, M2c.m.Ω, M2c.Λ.Ω, M2c.ϕ, M2c.aˍt
    Da = Differential(a)
    @test Set(equations(M2)) == Set([
        aˍt ~ ȧ # dummy equation
        Ω ~ Ωr + Ωm + ΩΛ
        ȧ ~ √(Ω) * a^2
        Da(aˍt)*Da(ϕ)*aˍt + aˍt^2*(Da^2)(ϕ) ~ -3*aˍt^2/a*Da(ϕ)
        aˍt*Da(Ωr) ~ -4*Ωr*aˍt/a
        aˍt*Da(Ωm) ~ -3*Ωm*aˍt/a
        aˍt*Da(ΩΛ) ~ 0
    ])

    @variables b(M2.a)
    M3 = change_independent_variable(M2, b, [Differential(M2.a)(b) ~ exp(-b), M2.a ~ exp(b)])

    M1 = structural_simplify(M1)
    M2 = structural_simplify(M2; allow_symbolic = true)
    M3 = structural_simplify(M3; allow_symbolic = true)
    @test length(unknowns(M3)) == length(unknowns(M2)) == length(unknowns(M1)) - 1
end

@testset "Change independent variable (simple)" begin
    @variables x(t) y1(t) # y(t)[1:1] # TODO: use array variables y(t)[1:2] when fixed: https://github.com/JuliaSymbolics/Symbolics.jl/issues/1383
    Mt = ODESystem([D(x) ~ 2*x, D(y1) ~ y1], t; name = :M)
    Mx = change_independent_variable(Mt, x)
    @variables x xˍt(x) xˍtt(x) y1(x) # y(x)[1:1] # TODO: array variables
    Dx = Differential(x)
    @test (Set(equations(Mx)) == Set([xˍt ~ 2*x, xˍt*Dx(y1) ~ y1]))
end

@testset "Change independent variable (free fall with 1st order horizontal equation)" begin
    @variables x(t) y(t)
    @parameters g v # gravitational acceleration and constant horizontal velocity
    Mt = ODESystem([D(D(y)) ~ -g, D(x) ~ v], t; name = :M) # gives (x, y) as function of t, ...
    Mx = change_independent_variable(Mt, x; add_old_diff = true) # ... but we want y as a function of x
    Mx = structural_simplify(Mx; allow_symbolic = true)
    Dx = Differential(Mx.x)
    prob = ODEProblem(Mx, [Mx.y => 0.0, Dx(Mx.y) => 1.0, Mx.t => 0.0], (0.0, 20.0), [g => 9.81, v => 10.0]) # 1 = dy/dx = (dy/dt)/(dx/dt) means equal initial horizontal and vertical velocities
    sol = solve(prob, Tsit5(); reltol = 1e-5)
    @test all(isapprox.(sol[Mx.y], sol[Mx.x - g*(Mx.t)^2/2]; atol = 1e-10)) # compare to analytical solution (x(t) = v*t, y(t) = v*t - g*t^2/2)
end

@testset "Change independent variable (free fall with 2nd order horizontal equation)" begin
    @variables x(t) y(t)
    @parameters g # gravitational acceleration
    Mt = ODESystem([D(D(y)) ~ -g, D(D(x)) ~ 0], t; name = :M) # gives (x, y) as function of t, ...
    Mx = change_independent_variable(Mt, x; add_old_diff = true) # ... but we want y as a function of x
    Mx = structural_simplify(Mx; allow_symbolic = true)
    Dx = Differential(Mx.x)
    prob = ODEProblem(Mx, [Mx.y => 0.0, Dx(Mx.y) => 1.0, Mx.t => 0.0, Mx.xˍt => 10.0], (0.0, 20.0), [g => 9.81]) # 1 = dy/dx = (dy/dt)/(dx/dt) means equal initial horizontal and vertical velocities
    sol = solve(prob, Tsit5(); reltol = 1e-5)
    @test all(isapprox.(sol[Mx.y], sol[Mx.x - g*(Mx.t)^2/2]; atol = 1e-10)) # compare to analytical solution (x(t) = v*t, y(t) = v*t - g*t^2/2)
end

@testset "Change independent variable (crazy analytical example)" begin
    @independent_variables t
    D = Differential(t)
    @variables x(t) y(t)
    M1 = ODESystem([ # crazy non-autonomous non-linear 2nd order ODE
        D(D(y)) ~ D(x)^2 + D(y^3) |> expand_derivatives # expand D(y^3) # TODO: make this test 3rd order
        D(x) ~ x^4 + y^5 + t^6
    ], t; name = :M)
    M2 = change_independent_variable(M1, x)

    # Compare to pen-and-paper result
    @independent_variables x
    Dx = Differential(x)
    @variables xˍt(x) xˍtt(x) y(x) t(x)
    @test Set(equations(M2)) == Set([
        xˍt^2*(Dx^2)(y) + xˍt*Dx(xˍt)*Dx(y) ~ xˍt^2 + 3*y^2*Dx(y)*xˍt # from D(D(y))
        xˍt ~ x^4 + y^5 + t^6 # dummy equation
    ])
end

@testset "Change independent variable (registered function / callable parameter)" begin
    @independent_variables t
    D = Differential(t)
    @variables x(t) y(t)
    @parameters f::LinearInterpolation (fc::LinearInterpolation)(..) # non-callable and callable
    callme(interp::LinearInterpolation, input) = interp(input)
    @register_symbolic callme(interp::LinearInterpolation, input)
    M1 = ODESystem([
        D(x) ~ 2*t
        D(y) ~ 1*fc(t) + 2*fc(x) + 3*fc(y) + 1*callme(f, t) + 2*callme(f, x) + 3*callme(f, y)
    ], t; name = :M)

    # Ensure that interpolations are called with the same variables
    M2 = change_independent_variable(M1, x, [t ~ √(x)])
    @variables x xˍt(x) y(x) t(x)
    Dx = Differential(x)
    @test Set(equations(M2)) == Set([
        t ~ √(x)
        xˍt ~ 2*t
        xˍt*Dx(y) ~ 1*fc(t) + 2*fc(x) + 3*fc(y) + 1*callme(f, t) + 2*callme(f, x) + 3*callme(f, y)
    ])

    _f = LinearInterpolation([1.0, 1.0], [-100.0, +100.0]) # constant value 1
    M2s = structural_simplify(M2; allow_symbolic = true)
    prob = ODEProblem(M2s, [M2s.y => 0.0], (1.0, 4.0), [fc => _f, f => _f])
    sol = solve(prob, Tsit5(); abstol = 1e-5)
    @test isapprox(sol(4.0, idxs=M2.y), 12.0; atol = 1e-5) # Anal solution is D(y) ~ 12 => y(t) ~ 12*t + C => y(x) ~ 12*√(x) + C. With y(x=1)=0 => 12*(√(x)-1), so y(x=4) ~ 12
end

@testset "Change independent variable (errors)" begin
    @variables x(t) y z(y) w(t) v(t)
    M = ODESystem([D(x) ~ 1, v ~ x], t; name = :M)
    @test_throws "structurally simplified" change_independent_variable(structural_simplify(M), y)
    @test_throws "not a function of the independent variable" change_independent_variable(M, y)
    @test_throws "not a function of the independent variable" change_independent_variable(M, z)
    @variables x(..) # require explicit argument
    M = ODESystem([D(x(t)) ~ x(t-1)], t; name = :M)
    @test_throws "DDE" change_independent_variable(M, x(t))
end
