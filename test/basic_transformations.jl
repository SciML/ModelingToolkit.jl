using ModelingToolkit, OrdinaryDiffEq, DataInterpolations, DynamicQuantities, Test
using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput
using SymbolicUtils: symtype

@independent_variables t
D = Differential(t)

@testset "Liouville transform" begin
    @parameters α β γ δ
    @variables x(t) y(t)
    eqs = [D(x) ~ α * x - β * x * y, D(y) ~ -δ * y + γ * x * y]
    @named sys = System(eqs, t)
    sys = complete(sys)

    u0 = [x => 1.0, y => 1.0]
    p = [α => 1.5, β => 1.0, δ => 3.0, γ => 1.0]
    tspan = (0.0, 10.0)
    prob = ODEProblem(sys, [u0; p], tspan)
    sol = solve(prob, Tsit5())

    sys2 = liouville_transform(sys)
    sys2 = complete(sys2)

    u0 = [x => 1.0, y => 1.0, sys2.trJ => 1.0]
    prob = ODEProblem(sys2, [u0; p], tspan, jac = true)
    sol = solve(prob, Tsit5())
    @test sol[end, end] ≈ 1.0742818931017244
end

@testset "Change independent variable (trivial)" begin
    @variables x(t) y(t)
    eqs1 = [D(D(x)) ~ D(x) + x, D(y) ~ 1]
    M1 = System(eqs1, t; name = :M)
    M2 = change_independent_variable(M1, y)
    @variables y x(y) yˍt(y)
    Dy = Differential(y)
    @test Set(equations(M2)) == Set([
        yˍt^2 * (Dy^2)(x) + yˍt * Dy(yˍt) * Dy(x) ~ x + Dy(x) * yˍt,
        yˍt ~ 1
    ])
end

@testset "Change independent variable" begin
    @variables x(t) y(t) z(t) s(t)
    eqs = [
        D(x) ~ y,
        D(D(y)) ~ 2 * x * D(y),
        z ~ x + D(y),
        D(s) ~ 1 / (2 * s)
    ]
    initialization_eqs = [x ~ 1.0, y ~ 1.0, D(y) ~ 0.0]
    M1 = System(eqs, t; initialization_eqs, name = :M)
    M2 = change_independent_variable(M1, s)

    M1 = mtkcompile(M1; allow_symbolic = true)
    M2 = mtkcompile(M2; allow_symbolic = true)
    prob1 = ODEProblem(M1, [M1.s => 1.0], (1.0, 4.0))
    prob2 = ODEProblem(M2, [], (1.0, 2.0))
    sol1 = solve(prob1, Tsit5(); reltol = 1e-10, abstol = 1e-10)
    sol2 = solve(prob2, Tsit5(); reltol = 1e-10, abstol = 1e-10)
    ts = range(0.0, 1.0, length = 50)
    ss = .√(ts)
    @test all(isapprox.(sol1(ts, idxs = M1.x), sol2(ss, idxs = M2.x); atol = 1e-7)) &&
          all(isapprox.(sol1(ts, idxs = M1.y), sol2(ss, idxs = M2.y); atol = 1e-7))
end

@testset "Change independent variable (Friedmann equation)" begin
    @independent_variables t
    D = Differential(t)
    @variables a(t) ȧ(t) Ω(t) ϕ(t)
    a, ȧ = GlobalScope.([a, ȧ])
    species(w; kw...) = System([D(Ω) ~ -3(1 + w) * D(a) / a * Ω], t, [Ω], []; kw...)
    @named r = species(1 // 3)
    @named m = species(0)
    @named Λ = species(-1)
    eqs = [
        Ω ~ r.Ω + m.Ω + Λ.Ω,
        D(a) ~ ȧ,
        ȧ ~ √(Ω) * a^2,
        D(D(ϕ)) ~ -3 * D(a) / a * D(ϕ)
    ]
    M1 = System(eqs, t, [Ω, a, ȧ, ϕ], []; name = :M)
    M1 = compose(M1, r, m, Λ)

    # Apply in two steps, where derivatives are defined at each step: first t -> a, then a -> b
    M2 = change_independent_variable(M1, M1.a)
    M2c = complete(M2) # just for the following equation comparison (without namespacing)
    a, ȧ, Ω, ϕ, aˍt = M2c.a, M2c.ȧ, M2c.Ω, M2c.ϕ, M2c.aˍt
    Ωr, Ωm, ΩΛ = M2c.r.Ω, M2c.m.Ω, M2c.Λ.Ω
    Da = Differential(a)
    @test Set(equations(M2)) == Set([
        aˍt ~ ȧ, # dummy equation
        Ω ~ Ωr + Ωm + ΩΛ,
        ȧ ~ √(Ω) * a^2,
        Da(aˍt) * Da(ϕ) * aˍt + aˍt^2 * (Da^2)(ϕ) ~ -3 * aˍt^2 / a * Da(ϕ),
        aˍt * Da(Ωr) ~ -4 * Ωr * aˍt / a,
        aˍt * Da(Ωm) ~ -3 * Ωm * aˍt / a,
        aˍt * Da(ΩΛ) ~ 0
    ])

    @variables b(M2.a)
    extraeqs = [Differential(M2.a)(b) ~ exp(-b), M2.a ~ exp(b)]
    M3 = change_independent_variable(M2, b, extraeqs)

    M1 = mtkcompile(M1)
    M2 = mtkcompile(M2; allow_symbolic = true)
    M3 = mtkcompile(M3; allow_symbolic = true)
    @test length(unknowns(M3)) == length(unknowns(M2)) == length(unknowns(M1)) - 1
end

@testset "Change independent variable (simple)" begin
    @variables x(t) y1(t) # y(t)[1:1] # TODO: use array variables y(t)[1:2] when fixed: https://github.com/JuliaSymbolics/Symbolics.jl/issues/1383
    Mt = System([D(x) ~ 2 * x, D(y1) ~ y1], t; name = :M)
    Mx = change_independent_variable(Mt, x)
    @variables x xˍt(x) xˍtt(x) y1(x) # y(x)[1:1] # TODO: array variables
    Dx = Differential(x)
    @test Set(equations(Mx)) == Set([xˍt ~ 2 * x, xˍt * Dx(y1) ~ y1])
end

@testset "Change independent variable (free fall with 1st order horizontal equation)" begin
    @variables x(t) y(t)
    @parameters g=9.81 v # gravitational acceleration and constant horizontal velocity
    Mt = System([D(D(y)) ~ -g, D(x) ~ v], t; name = :M) # gives (x, y) as function of t, ...
    Mx = change_independent_variable(Mt, x; add_old_diff = true) # ... but we want y as a function of x
    Mx = mtkcompile(Mx; allow_symbolic = true)
    Dx = Differential(Mx.x)
    u0 = [Mx.y => 0.0, Dx(Mx.y) => 1.0, Mx.t => 0.0]
    p = [v => 10.0]
    prob = ODEProblem(Mx, [u0; p], (0.0, 20.0)) # 1 = dy/dx = (dy/dt)/(dx/dt) means equal initial horizontal and vertical velocities
    sol = solve(prob, Tsit5(); reltol = 1e-5)
    @test all(isapprox.(sol[Mx.y], sol[Mx.x - g * (Mx.t) ^ 2 / 2]; atol = 1e-10)) # compare to analytical solution (x(t) = v*t, y(t) = v*t - g*t^2/2)
end

@testset "Change independent variable (free fall with 2nd order horizontal equation)" begin
    @variables x(t) y(t)
    @parameters g = 9.81 # gravitational acceleration
    Mt = System([D(D(y)) ~ -g, D(D(x)) ~ 0], t; name = :M) # gives (x, y) as function of t, ...
    Mx = change_independent_variable(Mt, x; add_old_diff = true) # ... but we want y as a function of x
    Mx = mtkcompile(Mx; allow_symbolic = true)
    Dx = Differential(Mx.x)
    u0 = [Mx.y => 0.0, Dx(Mx.y) => 1.0, Mx.t => 0.0, Mx.xˍt => 10.0]
    prob = ODEProblem(Mx, u0, (0.0, 20.0)) # 1 = dy/dx = (dy/dt)/(dx/dt) means equal initial horizontal and vertical velocities
    sol = solve(prob, Tsit5(); reltol = 1e-5)
    @test all(isapprox.(sol[Mx.y], sol[Mx.x - g * (Mx.t) ^ 2 / 2]; atol = 1e-10)) # compare to analytical solution (x(t) = v*t, y(t) = v*t - g*t^2/2)
end

@testset "Change independent variable (crazy 3rd order nonlinear system)" begin
    @independent_variables t
    D = Differential(t)
    @variables x(t) y(t)
    eqs = [
        (D^3)(y) ~ D(x)^2 + (D^2)(y^2) |> expand_derivatives,
        D(x)^2 + D(y)^2 ~ x^4 + y^5 + t^6
    ]
    M1 = System(eqs, t; name = :M)
    M2 = change_independent_variable(M1, x; add_old_diff = true)
    @test_nowarn mtkcompile(M2)

    # Compare to pen-and-paper result
    @variables x xˍt(x) xˍt(x) y(x) t(x)
    Dx = Differential(x)
    areequivalent(eq1,
        eq2) = isequal(expand(eq1.lhs - eq2.lhs), 0) &&
               isequal(expand(eq1.rhs - eq2.rhs), 0)
    eq1lhs = xˍt^3 * (Dx^3)(y) + xˍt^2 * Dx(y) * (Dx^2)(xˍt) +
             xˍt * Dx(y) * (Dx(xˍt))^2 +
             3 * xˍt^2 * (Dx^2)(y) * Dx(xˍt)
    eq1rhs = xˍt^2 + 2 * xˍt^2 * Dx(y)^2 +
             2 * xˍt^2 * y * (Dx^2)(y) +
             2 * y * Dx(y) * Dx(xˍt) * xˍt
    eq1 = eq1lhs ~ eq1rhs
    eq2 = xˍt^2 + xˍt^2 * Dx(y)^2 ~ x^4 + y^5 + t^6
    eq3 = Dx(t) ~ 1 / xˍt
    @test areequivalent(equations(M2)[1], eq1)
    @test areequivalent(equations(M2)[2], eq2)
    @test areequivalent(equations(M2)[3], eq3)
end

@testset "Change independent variable (registered function / callable parameter)" begin
    @independent_variables t
    D = Differential(t)
    @variables x(t) y(t)
    @parameters f::LinearInterpolation (fc::LinearInterpolation)(..) # non-callable and callable
    callme(interp::LinearInterpolation, input) = interp(input)
    @register_symbolic callme(interp::LinearInterpolation, input)
    eqs = [
        D(x) ~ 2t,
        D(y) ~ 1fc(t) + 2fc(x) + 3fc(y) + 1callme(f, t) + 2callme(f, x) + 3callme(f, y)
    ]
    M1 = System(eqs, t; name = :M)

    # Ensure that interpolations are called with the same variables
    M2 = change_independent_variable(M1, x, [t ~ √(x)])
    @variables x xˍt(x) y(x) t(x)
    Dx = Differential(x)
    @test Set(equations(M2)) == Set([
        t ~ √(x),
        xˍt ~ 2t,
        xˍt * Dx(y) ~
        1fc(t) + 2fc(x) + 3fc(y) +
        1callme(f, t) + 2callme(f, x) + 3callme(f, y)
    ])

    _f = LinearInterpolation([1.0, 1.0], [-100.0, +100.0]) # constant value 1
    M2s = mtkcompile(M2; allow_symbolic = true)
    prob = ODEProblem(M2s, [M2s.y => 0.0, fc => _f, f => _f], (1.0, 4.0))
    sol = solve(prob, Tsit5(); abstol = 1e-5)
    @test isapprox(sol(4.0, idxs = M2.y), 12.0; atol = 1e-5) # Anal solution is D(y) ~ 12 => y(t) ~ 12*t + C => y(x) ~ 12*√(x) + C. With y(x=1)=0 => 12*(√(x)-1), so y(x=4) ~ 12
end

@testset "Change independent variable (errors)" begin
    @variables x(t) y z(y) w(t) v(t)
    M = System([D(x) ~ 1, v ~ x], t; name = :M)
    Ms = mtkcompile(M)
    @test_throws "structurally simplified" change_independent_variable(Ms, y)
    @test_throws "not a function of" change_independent_variable(M, y)
    @test_throws "not a function of" change_independent_variable(M, z)
    @variables x(..) # require explicit argument
    M = System([D(x(t)) ~ x(t - 1)], t; name = :M)
    @test_throws "DDE" change_independent_variable(M, x(t))
end

@testset "Change independent variable w/ units (free fall with 2nd order horizontal equation)" begin
    @independent_variables t_units [unit = u"s"]
    D_units = Differential(t_units)
    @variables x(t_units) [unit = u"m"] y(t_units) [unit = u"m"]
    @parameters g=9.81 [unit = u"m * s^-2"] # gravitational acceleration
    Mt = System([D_units(D_units(y)) ~ -g, D_units(D_units(x)) ~ 0], t_units; name = :M) # gives (x, y) as function of t, ...
    Mx = change_independent_variable(Mt, x; add_old_diff = true) # ... but we want y as a function of x
    Mx = mtkcompile(Mx; allow_symbolic = true)
    Dx = Differential(Mx.x)
    u0 = [Mx.y => 0.0, Dx(Mx.y) => 1.0, Mx.t_units => 0.0, Mx.xˍt_units => 10.0]
    prob = ODEProblem(Mx, u0, (0.0, 20.0)) # 1 = dy/dx = (dy/dt)/(dx/dt) means equal initial horizontal and vertical velocities
    sol = solve(prob, Tsit5(); reltol = 1e-5)
    # compare to analytical solution (x(t) = v*t, y(t) = v*t - g*t^2/2)
    @test all(isapprox.(sol[Mx.y], sol[Mx.x - g * (Mx.t_units) ^ 2 / 2]; atol = 1e-10))
end

@testset "Change independent variable, no equations" begin
    # make this "look" like the standard library RealInput
    @mtkmodel Input begin
        @variables begin
            u(t)
        end
    end
    @named input_sys = Input()
    input_sys = complete(input_sys)
    # test no failures
    @test change_independent_variable(input_sys, input_sys.u) isa System

    @mtkmodel NestedInput begin
        @components begin
            in = Input()
        end
        @variables begin
            x(t)
        end
        @equations begin
            D(x) ~ in.u
        end
    end
    @named nested_input_sys = NestedInput()
    nested_input_sys = complete(nested_input_sys; flatten = false)
    @test change_independent_variable(nested_input_sys, nested_input_sys.x) isa System
end

@testset "Change of variables, connections" begin
    @mtkmodel ConnectSys begin
        @components begin
            in = RealInput()
            out = RealOutput()
        end
        @variables begin
            x(t)
            y(t)
        end
        @equations begin
            connect(in, out)
            in.u ~ x
            D(x) ~ -out.u
            D(y) ~ 1
        end
    end
    @named sys = ConnectSys()
    sys = complete(sys; flatten = false)
    new_sys = change_independent_variable(sys, sys.y; add_old_diff = true)
    ss = mtkcompile(new_sys; allow_symbolic = true)
    prob = ODEProblem(ss, [ss.t => 0.0, ss.x => 1.0], (0.0, 1.0))
    sol = solve(prob, Tsit5(); reltol = 1e-5)
    @test all(isapprox.(sol[ss.t], sol[ss.y]; atol = 1e-10))
    @test all(sol[ss.x][2:end] .< sol[ss.x][1])
end

@testset "Change independent variable with array variables" begin
    @variables x(t) y(t) z(t)[1:2]
    eqs = [
        D(x) ~ 2,
        z ~ ModelingToolkit.scalarize.([sin(y), cos(y)]),
        D(y) ~ z[1]^2 + z[2]^2
    ]
    @named sys = System(eqs, t)
    sys = complete(sys)
    new_sys = change_independent_variable(sys, sys.x; add_old_diff = true)
    ss_new_sys = mtkcompile(new_sys; allow_symbolic = true)
    u0 = [new_sys.y => 0.5, new_sys.t => 0.0]
    prob = ODEProblem(ss_new_sys, u0, (0.0, 0.5))
    sol = solve(prob, Tsit5(); reltol = 1e-5)
    @test sol[new_sys.y][end] ≈ 0.75
end

@testset "`add_accumulations`" begin
    @parameters a
    @variables x(t) y(t) z(t)
    @named sys = System([D(x) ~ y, 0 ~ x + z, 0 ~ x - y], t, [z, y, x], [])
    asys = add_accumulations(sys)
    @variables accumulation_x(t) accumulation_y(t) accumulation_z(t)
    eqs = [0 ~ x + z
           0 ~ x - y
           D(accumulation_x) ~ x
           D(accumulation_y) ~ y
           D(accumulation_z) ~ z
           D(x) ~ y]
    @test issetequal(equations(asys), eqs)
    @variables ac(t)
    asys = add_accumulations(sys, [ac => (x + y)^2])
    eqs = [0 ~ x + z
           0 ~ x - y
           D(ac) ~ (x + y)^2
           D(x) ~ y]
    @test issetequal(equations(asys), eqs)
end

abstract type AbstractFoo end

struct Bar <: AbstractFoo end
struct Baz <: AbstractFoo end

foofn(x) = 4
@register_symbolic foofn(x::AbstractFoo)

@testset "`respecialize`" begin
    @parameters p::AbstractFoo p2(t)::AbstractFoo = p q[1:2]::AbstractFoo r
    rp = only(let p = nothing
        @parameters p::Bar
    end)
    rp2 = only(let p2 = nothing
        @parameters p2(t)::Baz
    end)
    @variables x(t) = 1.0
    @named sys1 = System([D(x) ~ foofn(p) + foofn(p2) + x], t, [x], [p, p2, q, r])

    @test_throws ["completed systems"] respecialize(sys1)
    @test_throws ["completed systems"] respecialize(sys1, [])
    @test_throws ["split systems"] respecialize(complete(sys1; split = false))
    @test_throws ["split systems"] respecialize(complete(sys1; split = false), [])

    sys = complete(sys1)

    @test_throws ["Parameter p", "associated value"] respecialize(sys)
    @test_throws ["Parameter p", "associated value"] respecialize(sys, [p])

    @test_throws ["Parameter p2", "symbolic default"] respecialize(sys, [p2])

    sys2 = respecialize(sys, [p => Bar()])
    @test ModelingToolkit.iscomplete(sys2)
    @test ModelingToolkit.is_split(sys2)
    ps = ModelingToolkit.get_ps(sys2)
    idx = findfirst(isequal(rp), ps)
    @test defaults(sys2)[rp] == Bar()
    @test symtype(ps[idx]) <: Bar
    ic = ModelingToolkit.get_index_cache(sys2)
    @test any(x -> x.type == Bar && x.length == 1, ic.nonnumeric_buffer_sizes)
    prob = ODEProblem(sys2, [p2 => Bar(), q => [Bar(), Bar()], r => 1], (0.0, 1.0))
    @test any(x -> x isa Vector{Bar} && length(x) == 1, prob.p.nonnumeric)

    defaults(sys)[p2] = Baz()
    sys2 = respecialize(sys, [p => Bar()]; all = true)
    @test ModelingToolkit.iscomplete(sys2)
    @test ModelingToolkit.is_split(sys2)
    ps = ModelingToolkit.get_ps(sys2)
    idx = findfirst(isequal(rp2), ps)
    @test defaults(sys2)[rp2] == Baz()
    @test symtype(ps[idx]) <: Baz
    ic = ModelingToolkit.get_index_cache(sys2)
    @test any(x -> x.type == Baz && x.length == 1, ic.nonnumeric_buffer_sizes)
    delete!(defaults(sys), p2)
    prob = ODEProblem(sys2, [q => [Bar(), Bar()], r => 1], (0.0, 1.0))
    @test any(x -> x isa Vector{Bar} && length(x) == 1, prob.p.nonnumeric)
    @test any(x -> x isa Vector{Baz} && length(x) == 1, prob.p.nonnumeric)

    @test_throws ["Numeric types cannot be respecialized"] respecialize(sys, [r => 1])
    @test_throws ["array symbolics"] respecialize(sys, [q => Bar[Bar(), Bar()]])
    @test_throws ["scalarized array"] respecialize(sys, [q[1] => Bar()])

    @parameters foo::AbstractFoo
    @test_throws ["does not exist"] respecialize(sys, [foo => Bar()])
end
