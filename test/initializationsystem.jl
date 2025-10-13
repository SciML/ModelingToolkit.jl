using ModelingToolkit, OrdinaryDiffEq, NonlinearSolve, Test
using StochasticDiffEq, DelayDiffEq, StochasticDelayDiffEq, JumpProcesses
using ForwardDiff, StaticArrays
using SymbolicIndexingInterface, SciMLStructures
using SciMLStructures: Tunable
using ModelingToolkit: t_nounits as t, D_nounits as D, observed
using DynamicQuantities
using DiffEqBase: BrownFullBasicInit

@parameters g
@variables x(t) y(t) [state_priority = 10] λ(t)
eqs = [D(D(x)) ~ λ * x
       D(D(y)) ~ λ * y - g
       x^2 + y^2 ~ 1]
@mtkcompile pend = System(eqs, t)

initprob = ModelingToolkit.InitializationProblem(pend, 0.0, [g => 1];
    guesses = [ModelingToolkit.missing_variable_defaults(pend); x => 1; y => 0.2])
conditions = getfield.(equations(initprob.f.sys), :rhs)

@test initprob isa NonlinearLeastSquaresProblem
sol = solve(initprob)
@test SciMLBase.successful_retcode(sol)
@test maximum(abs.(sol[conditions])) < 1e-14

@test_throws ModelingToolkit.ExtraVariablesSystemException ModelingToolkit.InitializationProblem(
    pend, 0.0, [g => 1];
    guesses = [ModelingToolkit.missing_variable_defaults(pend); x => 1; y => 0.2],
    fully_determined = true)

initprob = ModelingToolkit.InitializationProblem(pend, 0.0, [x => 1, y => 0, g => 1];
    guesses = ModelingToolkit.missing_variable_defaults(pend))
@test initprob isa NonlinearLeastSquaresProblem
sol = solve(initprob)
@test SciMLBase.successful_retcode(sol)
@test all(iszero, sol.u)
@test maximum(abs.(sol[conditions])) < 1e-14

initprob = ModelingToolkit.InitializationProblem(
    pend, 0.0, [g => 1]; guesses = ModelingToolkit.missing_variable_defaults(pend))
@test initprob isa NonlinearLeastSquaresProblem
sol = solve(initprob)
@test !SciMLBase.successful_retcode(sol) ||
      sol.retcode == SciMLBase.ReturnCode.StalledSuccess

@test_throws ModelingToolkit.ExtraVariablesSystemException ModelingToolkit.InitializationProblem(
    pend, 0.0, [g => 1]; guesses = ModelingToolkit.missing_variable_defaults(pend),
    fully_determined = true)

prob = ODEProblem(pend, [x => 1, y => 0, g => 1], (0.0, 1.5),
    guesses = ModelingToolkit.missing_variable_defaults(pend))
prob.f.initializeprob isa NonlinearProblem
sol = solve(prob.f.initializeprob)
@test maximum(abs.(sol[conditions])) < 1e-14
sol = solve(prob, Rodas5P())
@test maximum(abs.(sol[conditions][1])) < 1e-14

prob = ODEProblem(pend, [x => 1, g => 1], (0.0, 1.5),
    guesses = ModelingToolkit.missing_variable_defaults(pend))
prob.f.initializeprob isa NonlinearLeastSquaresProblem
sol = solve(prob.f.initializeprob)
@test maximum(abs.(sol[conditions])) < 1e-14
sol = solve(prob, Rodas5P())
@test maximum(abs.(sol[conditions][1])) < 1e-14

@test_throws ModelingToolkit.ExtraVariablesSystemException ODEProblem(
    pend, [x => 1, g => 1], (0.0, 1.5),
    guesses = ModelingToolkit.missing_variable_defaults(pend),
    fully_determined = true)

@connector Port begin
    p(t)
    dm(t) = 0, [connect = Flow]
end

@connector Flange begin
    dx(t) = 0
    f(t), [connect = Flow]
end

# Components ----
@mtkmodel Orifice begin
    @parameters begin
        Cₒ = 2.7
        Aₒ = 0.00094
        ρ₀ = 1000
        p′ = 0
    end
    @variables begin
        dm(t) = 0
        p₁(t) = p′
        p₂(t) = p′
    end
    @components begin
        port₁ = Port(p = p′)
        port₂ = Port(p = p′)
    end
    begin
        u = dm / (ρ₀ * Aₒ)
    end
    @equations begin
        dm ~ +port₁.dm
        dm ~ -port₂.dm
        p₁ ~ port₁.p
        p₂ ~ port₂.p

        p₁ - p₂ ~ (1 / 2) * ρ₀ * u^2 * Cₒ
    end
end

@mtkmodel Volume begin
    @parameters begin
        A = 0.1
        ρ₀ = 1000
        β = 2e9
        direction = +1
        p′
        x′
    end
    @variables begin
        p(t)
        x(t) = x′
        dm(t) = 0
        f(t) = p′ * A
        dx(t) = 0
        r(t), [guess = 1000]
        dr(t), [guess = 1000]
    end
    @components begin
        port = Port(p = p′)
        flange = Flange(f = -p′ * A * direction)
    end
    @equations begin
        D(x) ~ dx
        D(r) ~ dr

        p ~ +port.p
        dm ~ +port.dm # mass is entering
        f ~ -flange.f * direction # force is leaving
        dx ~ flange.dx * direction

        r ~ ρ₀ * (1 + p / β)
        dm ~ (r * dx * A) + (dr * x * A)
        f ~ p * A
    end
end

@mtkmodel Mass begin
    @parameters begin
        m = 100
        f′
    end
    @variables begin
        f(t) = f′
        x(t) = 0
        dx(t) = 0
        ẍ(t) = f′ / m
    end
    @components begin
        flange = Flange(f = f′)
    end
    @equations begin
        D(x) ~ dx
        D(dx) ~ ẍ

        f ~ flange.f
        dx ~ flange.dx

        m * ẍ ~ f
    end
end

@mtkmodel Actuator begin
    @parameters begin
        p₁′
        p₂′
    end
    begin #constants
        x′ = 0.5
        A = 0.1
    end
    @components begin
        port₁ = Port(p = p₁′)
        port₂ = Port(p = p₂′)
        vol₁ = Volume(p′ = p₁′, x′ = x′, direction = -1)
        vol₂ = Volume(p′ = p₂′, x′ = x′, direction = +1)
        mass = Mass(f′ = (p₂′ - p₁′) * A)
        flange = Flange(f = 0)
    end
    @equations begin
        connect(port₁, vol₁.port)
        connect(port₂, vol₂.port)
        connect(vol₁.flange, vol₂.flange, mass.flange, flange)
    end
end

@mtkmodel Source begin
    @parameters begin
        p′
    end
    @components begin
        port = Port(p = p′)
    end
    @equations begin
        port.p ~ p′
    end
end

@mtkmodel Damper begin
    @parameters begin
        c = 1000
    end
    @components begin
        flange = Flange(f = 0)
    end
    @equations begin
        flange.f ~ c * flange.dx
    end
end

@mtkmodel HydraulicSystem begin
    @components begin
        res₁ = Orifice(p′ = 300e5)
        res₂ = Orifice(p′ = 0)
        act = Actuator(p₁′ = 300e5, p₂′ = 0)
        src = Source(p′ = 300e5)
        snk = Source(p′ = 0)
        dmp = Damper()
    end
    @equations begin
        connect(src.port, res₁.port₁)
        connect(res₁.port₂, act.port₁)
        connect(act.port₂, res₂.port₁)
        connect(res₂.port₂, snk.port)
        connect(dmp.flange, act.flange)
    end
end

@mtkcompile sys = HydraulicSystem()
initprob = ModelingToolkit.InitializationProblem(sys, 0.0)
conditions = getfield.(equations(initprob.f.sys), :rhs)

@test initprob isa NonlinearLeastSquaresProblem
@test length(initprob.u0) == 2
initsol = solve(initprob, reltol = 1e-12, abstol = 1e-12)
@test SciMLBase.successful_retcode(initsol)
@test maximum(abs.(initsol[conditions])) < 1e-14

@test_throws ModelingToolkit.ExtraEquationsSystemException ModelingToolkit.InitializationProblem(
    sys, 0.0, fully_determined = true)

allinit = unknowns(sys) .=> initsol[unknowns(sys)]
prob = ODEProblem(sys, allinit, (0, 0.1))
sol = solve(prob, Rodas5P(), initializealg = BrownFullBasicInit())
# If initialized incorrectly, then it would be InitialFailure
@test sol.retcode == SciMLBase.ReturnCode.Unstable
@test maximum(abs.(initsol[conditions][1])) < 1e-14

prob = ODEProblem(sys, allinit, (0, 0.1))
prob = ODEProblem(sys, [], (0, 0.1), check = false)

@test_throws ModelingToolkit.ExtraEquationsSystemException ODEProblem(
    sys, [], (0, 0.1), fully_determined = true)

sol = solve(prob, Rodas5P())
# If initialized incorrectly, then it would be InitialFailure
@test sol.retcode == SciMLBase.ReturnCode.Unstable
@test maximum(abs.(initsol[conditions][1])) < 1e-14

@connector Flange begin
    dx(t), [guess = 0]
    f(t), [guess = 0, connect = Flow]
end

@mtkmodel Mass begin
    @parameters begin
        m = 100
    end
    @variables begin
        dx(t)
        f(t), [guess = 0]
    end
    @components begin
        flange = Flange()
    end
    @equations begin
        # connectors
        flange.dx ~ dx
        flange.f ~ -f

        # physics
        f ~ m * D(dx)
    end
end

@mtkmodel Damper begin
    @parameters begin
        d = 1
    end
    @variables begin
        dx(t), [guess = 0]
        f(t), [guess = 0]
    end
    @components begin
        flange = Flange()
    end
    @equations begin
        # connectors
        flange.dx ~ dx
        flange.f ~ -f

        # physics
        f ~ d * dx
    end
end

@mtkmodel MassDamperSystem begin
    @components begin
        mass = Mass(; dx = 100, m = 10)
        damper = Damper(; d = 1)
    end
    @equations begin
        connect(mass.flange, damper.flange)
    end
end

@mtkcompile sys = MassDamperSystem()
initprob = ModelingToolkit.InitializationProblem(sys, 0.0)
@test initprob isa NonlinearProblem
initsol = solve(initprob, reltol = 1e-12, abstol = 1e-12)
@test SciMLBase.successful_retcode(initsol)

allinit = unknowns(sys) .=> initsol[unknowns(sys)]
prob = ODEProblem(sys, allinit, (0, 0.1))
sol = solve(prob, Rodas5P())
# If initialized incorrectly, then it would be InitialFailure
@test sol.retcode == SciMLBase.ReturnCode.Success

prob = ODEProblem(sys, [], (0, 0.1))
sol = solve(prob, Rodas5P())
@test sol.retcode == SciMLBase.ReturnCode.Success

### Ensure that non-DAEs still throw for missing variables without the initialize system

@parameters σ ρ β
@variables x(t) y(t) z(t)

eqs = [D(D(x)) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@mtkcompile sys = System(eqs, t)

u0 = [D(x) => 2.0,
    y => 0.0,
    z => 0.0]

p = [σ => 28.0,
    ρ => 10.0,
    β => 8 / 3]

tspan = (0.0, 100.0)
@test_throws ModelingToolkit.IncompleteInitializationError prob=ODEProblem(
    sys, [u0; p], tspan, jac = true)

u0 = [y => 0.0,
    z => 0.0]
@test_throws "Differential(t)(x(t))" prob=ODEProblem(
    sys, [u0; p], tspan, jac = true)

# DAE Initialization on ODE with nonlinear system for initial conditions
# https://github.com/SciML/ModelingToolkit.jl/issues/2508

using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

function System2(; name)
    vars = @variables begin
        dx(t), [guess = 0]
        ddx(t), [guess = 0]
    end
    eqs = [D(dx) ~ ddx
           0 ~ ddx + dx + 1]
    return System(eqs, t, vars, []; name)
end

@mtkcompile sys = System2()
prob = ODEProblem(sys, [sys.dx => 1], (0, 1)) # OK
prob = ODEProblem(sys, [sys.ddx => -2], (0, 1), guesses = [sys.dx => 1])
sol = solve(prob, Tsit5())
@test SciMLBase.successful_retcode(sol)
@test sol.u[1] == [1.0]

## Late binding initialization_eqs

function System3(; name)
    vars = @variables begin
        dx(t), [guess = 0]
        ddx(t), [guess = 0]
    end
    eqs = [D(dx) ~ ddx
           0 ~ ddx + dx + 1]
    initialization_eqs = [
        ddx ~ -2
    ]
    return System(eqs, t, vars, []; name, initialization_eqs)
end

@mtkcompile sys = System3()
prob = ODEProblem(sys, [], (0, 1), guesses = [sys.dx => 1])
sol = solve(prob, Tsit5())
@test SciMLBase.successful_retcode(sol)
@test sol.u[1] == [1.0]

# Steady state initialization
@testset "Steady state initialization" begin
    @parameters σ ρ β
    @variables x(t) y(t) z(t)

    eqs = [D(D(x)) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z]

    @named sys = System(eqs, t)
    sys = mtkcompile(sys)

    u0 = [D(x) => 2.0,
        x => 1.0,
        D(y) => 0.0,
        z => 0.0]

    p = [σ => 28.0,
        ρ => 10.0,
        β => 8 / 3]

    tspan = (0.0, 0.2)
    prob_mtk = ODEProblem(sys, [u0; p], tspan)
    sol = solve(prob_mtk, Tsit5())
    @test sol[x * (ρ - z) - y][1] == 0.0

    prob_mtk.ps[Initial(D(y))] = 1.0
    sol = solve(prob_mtk, Tsit5())
    @test sol[x * (ρ - z) - y][1] == 1.0
end

@variables x(t) y(t) z(t)
@parameters α=1.5 β=1.0 γ=3.0 δ=1.0

eqs = [D(x) ~ α * x - β * x * y
       D(y) ~ -γ * y + δ * x * y
       z ~ x + y]

@named sys = System(eqs, t)
simpsys = mtkcompile(sys)
tspan = (0.0, 10.0)

prob = ODEProblem(simpsys, [D(x) => 0.0, y => 0.0], tspan, guesses = [x => 0.0])
sol = solve(prob, Tsit5())
@test sol.u[1] == [0.0, 0.0]

# Initialize with an observed variable
prob = ODEProblem(simpsys, [z => 0.0], tspan, guesses = [x => 2.0, y => 4.0])
sol = solve(prob, Tsit5())
@test sol[z, 1] == 0.0

prob = ODEProblem(simpsys, [z => 1.0, y => 1.0], tspan, guesses = [x => 2.0])
sol = solve(prob, Tsit5())
@test sol[[x, y], 1] == [0.0, 1.0]

@test_warn "underdetermined" prob = ODEProblem(
    simpsys, [], tspan, guesses = [x => 2.0, y => 1.0])

# Late Binding initialization_eqs
# https://github.com/SciML/ModelingToolkit.jl/issues/2787

@parameters g
@variables x(t) y(t) [state_priority = 10] λ(t)
eqs = [D(D(x)) ~ λ * x
       D(D(y)) ~ λ * y - g
       x^2 + y^2 ~ 1]
@mtkcompile pend = System(eqs, t)

prob = ODEProblem(pend, [x => 1, g => 1], (0.0, 1.5),
    guesses = [λ => 0, y => 1], initialization_eqs = [y ~ 1])

unsimp = generate_initializesystem(pend; op = [x => 1], initialization_eqs = [y ~ 1])
sys = mtkcompile(unsimp; fully_determined = false)
@test length(equations(sys)) in (3, 4) # could be either depending on tearing

# Extend two systems with initialization equations and guesses
# https://github.com/SciML/ModelingToolkit.jl/issues/2845
@variables x(t) y(t)
@named sysx = System([D(x) ~ 0], t; initialization_eqs = [x ~ 1])
@named sysy = System([D(y) ~ 0], t; initialization_eqs = [y^2 ~ 2], guesses = [y => 1])
sys = complete(extend(sysx, sysy))
@test length(equations(generate_initializesystem(sys))) == 2
@test length(ModelingToolkit.guesses(sys)) == 1

# https://github.com/SciML/ModelingToolkit.jl/issues/2873
@testset "Error on missing defaults" begin
    @variables x(t) y(t)
    @named sys = System([x^2 + y^2 ~ 25, D(x) ~ 1], t)
    ssys = mtkcompile(sys)
    @test_throws ModelingToolkit.MissingGuessError ODEProblem(
        ssys, [x => 3], (0, 1)) # y should have a guess
end

# https://github.com/SciML/ModelingToolkit.jl/issues/3025
@testset "Override defaults when setting initial conditions with unknowns(sys) or similar" begin
    @variables x(t) y(t)

    # system 1 should solve to x = 1
    ics1 = [x => 1]
    sys1 = System([D(x) ~ 0], t; defaults = ics1, name = :sys1) |> mtkcompile
    prob1 = ODEProblem(sys1, [], (0.0, 1.0))
    sol1 = solve(prob1, Tsit5())
    @test all(sol1[x] .== 1)

    # system 2 should solve to x = y = 2
    sys2 = extend(
        sys1,
        System([D(y) ~ 0], t; initialization_eqs = [y ~ 2], name = :sys2)
    ) |> mtkcompile
    ics2 = unknowns(sys1) .=> 2 # should be equivalent to "ics2 = [x => 2]"
    prob2 = ODEProblem(sys2, ics2, (0.0, 1.0); fully_determined = true)
    sol2 = solve(prob2, Tsit5())
    @test all(sol2[x] .== 2) && all(sol2[y] .== 2)
end

# https://github.com/SciML/ModelingToolkit.jl/issues/3029
@testset "Derivatives in initialization equations" begin
    @variables x(t)
    sys = System(
        [D(D(x)) ~ 0], t; initialization_eqs = [x ~ 0, D(x) ~ 1], name = :sys) |>
          mtkcompile
    @test_nowarn ODEProblem(sys, [], (0.0, 1.0))

    sys = System(
        [D(D(x)) ~ 0], t; initialization_eqs = [x ~ 0, D(D(x)) ~ 0], name = :sys) |>
          mtkcompile
    @test_nowarn ODEProblem(sys, [D(x) => 1.0], (0.0, 1.0))
end

# https://github.com/SciML/ModelingToolkit.jl/issues/3049
@testset "Derivatives in initialization guesses" begin
    for sign in [-1.0, +1.0]
        @variables x(t)
        sys = System(
            [D(D(x)) ~ 0], t;
            initialization_eqs = [D(x)^2 ~ 1, x ~ 0], guesses = [D(x) => sign], name = :sys
        ) |> mtkcompile
        prob = ODEProblem(sys, [], (0.0, 1.0))
        sol = solve(prob, Tsit5())
        @test sol(1.0, idxs = sys.x) ≈ sign # system with D(x(0)) = ±1 should solve to x(1) = ±1
    end
end

# https://github.com/SciML/ModelingToolkit.jl/issues/2619
@parameters k1 k2 ω
@variables X(t) Y(t)
eqs_1st_order = [D(Y) + Y - ω ~ 0,
    X + k1 ~ Y + k2]
eqs_2nd_order = [D(D(Y)) + 2ω * D(Y) + (ω^2) * Y ~ 0,
    X + k1 ~ Y + k2]
@mtkcompile sys_1st_order = System(eqs_1st_order, t)
@mtkcompile sys_2nd_order = System(eqs_2nd_order, t)

u0_1st_order_1 = [X => 1.0, Y => 2.0]
u0_1st_order_2 = [Y => 2.0]
u0_2nd_order_1 = [X => 1.0, Y => 2.0, D(Y) => 0.5]
u0_2nd_order_2 = [Y => 2.0, D(Y) => 0.5]
tspan = (0.0, 10.0)
ps = [ω => 0.5, k1 => 2.0, k2 => 3.0]

oprob_1st_order_1 = ODEProblem(sys_1st_order, [u0_1st_order_1; ps], tspan)
oprob_1st_order_2 = ODEProblem(sys_1st_order, [u0_1st_order_2; ps], tspan)
oprob_2nd_order_1 = ODEProblem(sys_2nd_order, [u0_2nd_order_1; ps], tspan) # gives sys_2nd_order
oprob_2nd_order_2 = ODEProblem(sys_2nd_order, [u0_2nd_order_2; ps], tspan)

@test solve(oprob_1st_order_1, Rosenbrock23()).retcode ==
      SciMLBase.ReturnCode.InitialFailure
@test solve(oprob_1st_order_2, Rosenbrock23())[Y][1] == 2.0
@test solve(oprob_2nd_order_1, Rosenbrock23()).retcode ==
      SciMLBase.ReturnCode.InitialFailure
sol = solve(oprob_2nd_order_2, Rosenbrock23()) # retcode: Success
@test sol[Y][1] == 2.0
@test sol[D(Y)][1] == 0.5

@testset "Vector in initial conditions" begin
    @variables x(t)[1:5] y(t)[1:5]
    @named sys = System([D(x) ~ x, D(y) ~ y], t; initialization_eqs = [y ~ -x])
    sys = mtkcompile(sys)
    prob = ODEProblem(sys, [sys.x => ones(5)], (0.0, 1.0))
    sol = solve(prob, Tsit5(), reltol = 1e-4)
    @test all(sol(1.0, idxs = sys.x) .≈ +exp(1)) && all(sol(1.0, idxs = sys.y) .≈ -exp(1))
end

@testset "Initialization of parameters" begin
    @variables _x(..) y(t)
    @parameters p q
    @brownians a b
    x = _x(t)
    sarray_ctor = splat(SVector)
    # `System` constructor creates appropriate type with mtkcompile
    # `Problem` and `alg` create the problem to test and allow calling `init` with
    # the correct solver.
    # `rhss` allows adding terms to the end of equations (only 2 equations allowed) to influence
    # the system type (brownian vars to turn it into an SDE).
    @testset "$Problem with $(SciMLBase.parameterless_type(alg)) and $ctor ctor" for (
        (Problem, alg, rhss), (ctor, expectedT)) in Iterators.product(
        [
            (ODEProblem, Tsit5(), zeros(2)),
            (SDEProblem, ImplicitEM(), [a, b]),
            (DDEProblem, MethodOfSteps(Tsit5()), [_x(t - 0.1), 0.0]),
            (SDDEProblem, ImplicitEM(), [_x(t - 0.1) + a, b])
        ],
        [(identity, Any), (sarray_ctor, SVector)])
        u0_constructor = p_constructor = ctor
        if ctor !== identity
            Problem = Problem{false}
        end
        function test_parameter(prob, sym, val)
            if prob.u0 !== nothing
                @test prob.u0 isa expectedT
                @test init(prob, alg).ps[sym] ≈ val
            end
            @test prob.p.tunable isa expectedT
            initprob = prob.f.initialization_data.initializeprob
            if state_values(initprob) !== nothing
                @test state_values(initprob) isa expectedT
            end
            @test parameter_values(initprob).tunable isa expectedT
            @test solve(prob, alg).ps[sym] ≈ val
        end
        function test_initializesystem(prob, p, equation)
            isys = prob.f.initialization_data.initializeprob.f.sys
            @test is_variable(isys, p) || ModelingToolkit.has_observed_with_lhs(isys, p)
            @test equation in [equations(isys); observed(isys)]
        end

        u0map = Dict(x => 1.0, y => 1.0)
        pmap = Dict()
        pmap[q] = 1.0
        # `missing` default, equation from Problem
        @mtkcompile sys = System(
            [D(x) ~ x * q + rhss[1], D(y) ~ y * p + rhss[2]], t; defaults = [p => missing], guesses = [p => 1.0])
        pmap[p] = 2q
        prob = Problem(sys, merge(u0map, pmap), (0.0, 1.0); u0_constructor, p_constructor)
        test_parameter(prob, p, 2.0)
        prob2 = remake(prob; u0 = u0map, p = pmap)
        prob2 = remake(prob2; p = setp_oop(prob2, p)(prob2, 0.0))
        test_parameter(prob2, p, 2.0)
        # `missing` default, provided guess
        @mtkcompile sys = System(
            [D(x) ~ x + rhss[1], p ~ x + y + rhss[2]], t; defaults = [p => missing], guesses = [p => 0.0])
        prob = Problem(sys, u0map, (0.0, 1.0); u0_constructor, p_constructor)
        test_parameter(prob, p, 2.0)
        test_initializesystem(prob, p, p ~ x + y)
        prob2 = remake(prob; u0 = u0map)
        prob2 = remake(prob2; p = setp_oop(prob2, p)(prob2, 0.0))
        test_parameter(prob2, p, 2.0)

        # `missing` to Problem, equation from default
        @mtkcompile sys = System(
            [D(x) ~ x * q + rhss[1], D(y) ~ y * p + rhss[2]], t; defaults = [p => 2q], guesses = [p => 1.0])
        pmap[p] = missing
        prob = Problem(sys, merge(u0map, pmap), (0.0, 1.0); u0_constructor, p_constructor)
        test_parameter(prob, p, 2.0)
        test_initializesystem(prob, p, p ~ 2q)
        prob2 = remake(prob; u0 = u0map, p = pmap)
        prob2 = remake(prob2; p = setp_oop(prob2, p)(prob2, 0.0))
        test_parameter(prob2, p, 2.0)
        # `missing` to Problem, provided guess
        @mtkcompile sys = System(
            [D(x) ~ x + rhss[1], p ~ x + y + rhss[2]], t; guesses = [p => 0.0])
        prob = Problem(sys, merge(u0map, pmap), (0.0, 1.0); u0_constructor, p_constructor)
        test_parameter(prob, p, 2.0)
        test_initializesystem(prob, p, p ~ x + y)
        prob2 = remake(prob; u0 = u0map, p = pmap)
        prob2 = remake(prob2; p = setp_oop(prob2, p)(prob2, 0.0))
        test_parameter(prob2, p, 2.0)

        # No `missing`, default and guess
        @mtkcompile sys = System(
            [D(x) ~ x * q + rhss[1], D(y) ~ y * p + rhss[2]], t; defaults = [p => 2q], guesses = [p => 0.0])
        delete!(pmap, p)
        prob = Problem(sys, merge(u0map, pmap), (0.0, 1.0); u0_constructor, p_constructor)
        test_parameter(prob, p, 2.0)
        test_initializesystem(prob, p, p ~ 2q)
        prob2 = remake(prob; u0 = u0map, p = pmap)
        prob2 = remake(prob2; p = setp_oop(prob2, p)(prob2, 0.0))
        test_parameter(prob2, p, 2.0)

        # Default overridden by Problem, guess provided
        @mtkcompile sys = System(
            [D(x) ~ q * x + rhss[1], D(y) ~ y * p + rhss[2]], t; defaults = [p => 2q], guesses = [p => 1.0])
        _pmap = merge(pmap, Dict(p => q))
        prob = Problem(sys, merge(u0map, _pmap), (0.0, 1.0); u0_constructor, p_constructor)
        test_parameter(prob, p, _pmap[q])
        test_initializesystem(prob, p, p ~ q)
        # Problem dependent value with guess, no `missing`
        @mtkcompile sys = System(
            [D(x) ~ y * q + p + rhss[1], D(y) ~ x * p + q + rhss[2]], t; guesses = [p => 0.0])
        _pmap = merge(pmap, Dict(p => 3q))
        prob = Problem(sys, merge(u0map, _pmap), (0.0, 1.0); u0_constructor, p_constructor)
        test_parameter(prob, p, 3pmap[q])

        # Should not be solved for:
        # Override dependent default with direct value
        @mtkcompile sys = System(
            [D(x) ~ q * x + rhss[1], D(y) ~ y * p + rhss[2]], t; defaults = [p => 2q], guesses = [p => 1.0])
        _pmap = merge(pmap, Dict(p => 1.0))
        prob = Problem(sys, merge(u0map, _pmap), (0.0, 1.0); u0_constructor, p_constructor)
        @test prob.ps[p] ≈ 1.0
        initsys = prob.f.initialization_data.initializeprob.f.sys
        @test is_parameter(initsys, p)

        # Non-floating point
        @parameters r::Int s::Int
        @mtkcompile sys = System(
            [D(x) ~ s * x + rhss[1], D(y) ~ y * r + rhss[2]], t; defaults = [s => 2r], guesses = [s => 1.0])
        prob = Problem(
            sys, merge(u0map, Dict(r => 1)), (0.0, 1.0); u0_constructor, p_constructor)
        @test prob.ps[r] == 1
        @test prob.ps[s] == 2
        initsys = prob.f.initialization_data.initializeprob.f.sys
        @test is_parameter(initsys, r)
        @test is_parameter(initsys, s)

        @mtkcompile sys = System(
            [D(x) ~ x + rhss[1], p ~ x + y + rhss[2]], t; guesses = [p => 0.0])
        @test_throws ModelingToolkit.MissingParametersError Problem(
            sys, [x => 1.0, y => 1.0], (0.0, 1.0))

        # Unsatisfiable initialization
        prob = Problem(sys, [x => 1.0, y => 1.0, p => 2.0], (0.0, 1.0);
            initialization_eqs = [x^2 + y^2 ~ 3], u0_constructor, p_constructor)
        @test prob.f.initialization_data !== nothing
        @test solve(prob, alg).retcode == ReturnCode.InitialFailure
        cache = init(prob, alg)
        @test solve!(cache).retcode == ReturnCode.InitialFailure
    end

    @testset "Null system" begin
        @variables x(t) y(t) s(t)
        @parameters x0 y0
        @mtkcompile sys = System([x ~ x0, y ~ y0, s ~ x + y], t; guesses = [y0 => 0.0])
        prob = ODEProblem(sys, [s => 1.0, x0 => 0.3, y0 => missing], (0.0, 1.0))
        # trivial initialization run immediately
        @test prob.ps[y0] ≈ 0.7
        @test init(prob, Tsit5()).ps[y0] ≈ 0.7
        @test solve(prob, Tsit5()).ps[y0] ≈ 0.7
    end

    using ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica: Fixed, Mass,
                                                                           Spring, Force,
                                                                           Damper
    using ModelingToolkitStandardLibrary.Mechanical: TranslationalModelica as TM
    using ModelingToolkitStandardLibrary.Blocks: Constant

    @named mass = TM.Mass(; m = 1.0, s = 1.0, v = 0.0, a = 0.0)
    @named fixed = Fixed(; s0 = 0.0)
    @named spring = Spring(; c = 2.0, s_rel0 = nothing)
    @named gravity = Force()
    @named constant = Constant(; k = 9.81)
    @named damper = TM.Damper(; d = 0.1)
    @mtkcompile sys = System(
        [connect(fixed.flange, spring.flange_a), connect(spring.flange_b, mass.flange_a),
            connect(mass.flange_a, gravity.flange), connect(constant.output, gravity.f),
            connect(fixed.flange, damper.flange_a), connect(damper.flange_b, mass.flange_a)],
        t;
        systems = [fixed, spring, mass, gravity, constant, damper],
        guesses = [spring.s_rel0 => 1.0])
    prob = ODEProblem(sys, [spring.s_rel0 => missing], (0.0, 1.0))
    # trivial initialization run immediately
    @test prob.ps[spring.s_rel0] ≈ -3.905
    @test init(prob, Tsit5()).ps[spring.s_rel0] ≈ -3.905
    @test solve(prob, Tsit5()).ps[spring.s_rel0] ≈ -3.905
end

@testset "NonlinearSystem initialization" begin
    nl_algs = [FastShortcutNonlinearPolyalg(), NewtonRaphson(),
        Klement(), SimpleNewtonRaphson(), DFSane()]
    nlls_algs = [FastShortcutNLLSPolyalg(), LevenbergMarquardt(), SimpleGaussNewton()]

    @testset "No initialization for variables" begin
        @variables x=1.0 y=0.0 z=0.0
        @parameters σ=10.0 ρ=26.0 β=8/3

        eqs = [0 ~ σ * (y - x),
            0 ~ x * (ρ - z) - y,
            0 ~ x * y - β * z]
        @mtkcompile ns = System(eqs, [x, y, z], [σ, ρ, β])

        prob = NonlinearProblem(ns, [])
        @test prob.f.initialization_data.update_initializeprob! === nothing
        @test prob.f.initialization_data.initializeprobmap === nothing
        @test prob.f.initialization_data.initializeprobpmap === nothing
        for alg in nl_algs
            @test SciMLBase.successful_retcode(solve(prob, alg))
        end

        prob = NonlinearLeastSquaresProblem(ns, [])
        @test prob.f.initialization_data.update_initializeprob! === nothing
        @test prob.f.initialization_data.initializeprobmap === nothing
        @test prob.f.initialization_data.initializeprobpmap === nothing
        for alg in nlls_algs
            @test SciMLBase.successful_retcode(solve(prob, alg))
        end
    end

    prob_alg_combinations = zip(
        [NonlinearProblem, NonlinearLeastSquaresProblem], [nl_algs, nlls_algs])
    sarray_ctor = splat(SVector)
    @testset "Parameter initialization with ctor $ctor" for (ctor, expectedT) in [
        (identity, Any),
        (sarray_ctor, SVector)
    ]
        u0_constructor = p_constructor = ctor
        function test_parameter(prob, alg, param, val)
            if prob.u0 !== nothing
                @test prob.u0 isa expectedT
            end
            @test prob.p.tunable isa expectedT
            integ = init(prob, alg)
            @test integ.ps[param]≈val rtol=1e-5
            # some algorithms are a little temperamental
            sol = solve(prob, alg)
            @test sol.ps[param]≈val rtol=1e-5 broken=(alg===SimpleNewtonRaphson())
            @test SciMLBase.successful_retcode(sol)
        end

        @parameters p=2.0 q=missing [guess=1.0] c=1.0
        @variables x=1.0 z=3.0

        # eqs = [0 ~ p * (y - x),
        #     0 ~ x * (q - z) - y,
        #     0 ~ x * y - c * z]
        # specifically written this way due to
        # https://github.com/SciML/NonlinearSolve.jl/issues/586
        eqs = [0 ~ -c * z + (q - z) * (x^2)
               0 ~ p * (-x + (q - z) * x)]
        @named sys = System(eqs; initialization_eqs = [p^2 + q^2 + 2p * q ~ 0])
        sys = complete(sys)
        # @mtkcompile sys = NonlinearSystem(
        #     [p * x^2 + q * y^3 ~ 0, x - q ~ 0]; defaults = [q => missing],
        #     guesses = [q => 1.0], initialization_eqs = [p^2 + q^2 + 2p * q ~ 0])

        for (probT, algs) in prob_alg_combinations
            if ctor != identity
                probT = probT{false}
            end
            prob = probT(sys, []; u0_constructor, p_constructor)
            @test prob.f.initialization_data !== nothing
            @test prob.f.initialization_data.initializeprobmap === nothing
            for alg in algs
                test_parameter(prob, alg, q, -2.0)
            end

            # `update_initializeprob!` works
            prob = remake(prob; p = setp_oop(prob, p)(prob, -2.0))
            for alg in algs
                test_parameter(prob, alg, q, 2.0)
            end
            prob = remake(prob; p = setp_oop(prob, p)(prob, 2.0))

            # `remake` works
            prob2 = remake(prob; p = [p => -2.0])
            @test prob2.f.initialization_data !== nothing
            @test prob2.f.initialization_data.initializeprobmap === nothing
            for alg in algs
                test_parameter(prob2, alg, q, 2.0)
            end

            # changing types works
            ps = parameter_values(prob)
            newps = SciMLStructures.replace(Tunable(), ps, ForwardDiff.Dual.(ps.tunable))
            prob3 = remake(prob; p = newps)
            @test prob3.f.initialization_data !== nothing
            @test eltype(state_values(prob3.f.initialization_data.initializeprob)) <:
                  ForwardDiff.Dual
            @test eltype(prob3.f.initialization_data.initializeprob.p.tunable) <:
                  ForwardDiff.Dual
        end
    end
end

@testset "Update initializeprob parameters" begin
    @variables _x(..) y(t)
    @parameters p q
    @brownians a b
    x = _x(t)

    @testset "$Problem with $(SciMLBase.parameterless_type(typeof(alg)))" for (
        System, Problem, alg, rhss) in [
        (ModelingToolkit.System, ODEProblem, Tsit5(), zeros(2)),
        (ModelingToolkit.System, SDEProblem, ImplicitEM(), [a, b]),
        (ModelingToolkit.System, DDEProblem, MethodOfSteps(Tsit5()), [_x(t - 0.1), 0.0]),
        (ModelingToolkit.System, SDDEProblem, ImplicitEM(), [_x(t - 0.1) + a, b])
    ]
        @mtkcompile sys = System(
            [D(x) ~ x + rhss[1], p ~ x + y + rhss[2]], t; guesses = [x => 0.0, p => 0.0])
        prob = Problem(sys, [y => 1.0, p => 3.0], (0.0, 1.0))
        @test prob.f.initialization_data.initializeprob.ps[p] ≈ 3.0
        @test init(prob, alg)[x] ≈ 2.0
        prob.ps[p] = 2.0
        @test prob.f.initialization_data.initializeprob.ps[p] ≈ 3.0
        @test init(prob, alg)[x] ≈ 1.0
        ModelingToolkit.defaults(prob.f.sys)[p] = missing
        prob2 = remake(prob; u0 = [y => 1.0], p = [p => 3x])
        @test !is_variable(prob2.f.initialization_data.initializeprob, p) &&
              !is_parameter(prob2.f.initialization_data.initializeprob, p)
        @test init(prob2, alg)[x] ≈ 0.5
        @test_nowarn solve(prob2, alg)
    end
end

@testset "Equations for dependent parameters" begin
    @variables _x(..)
    @parameters p q=5 r
    @brownians a
    x = _x(t)

    @testset "$Problem with $(SciMLBase.parameterless_type(typeof(alg)))" for (
        System, Problem, alg, rhss) in [
        (ModelingToolkit.System, ODEProblem, Tsit5(), 0),
        (ModelingToolkit.System, SDEProblem, ImplicitEM(), a),
        (ModelingToolkit.System, DDEProblem, MethodOfSteps(Tsit5()), _x(t - 0.1)),
        (ModelingToolkit.System, SDDEProblem, ImplicitEM(), _x(t - 0.1) + a)
    ]
        @mtkcompile sys = System(
            [D(x) ~ 2x + r + rhss, r ~ p + 2q, q ~ p + 3], t;
            guesses = [p => 1.0])
        prob = Problem(sys, [x => 1.0, p => missing], (0.0, 1.0))
        parent_isys = ModelingToolkit.get_parent(prob.f.initialization_data.initializeprob.f.sys)
        @test length(equations(parent_isys)) == 4
        integ = init(prob, alg)
        @test integ.ps[p] ≈ 2
    end
end

@testset "Re-creating initialization problem on remake" begin
    @variables _x(..) y(t)
    @parameters p q
    @brownians a b
    x = _x(t)

    @testset "$Problem with $(SciMLBase.parameterless_type(typeof(alg)))" for (
        Problem, alg, rhss) in [
        (ODEProblem, Tsit5(), zeros(2)),
        (SDEProblem, ImplicitEM(), [a, b]),
        (DDEProblem, MethodOfSteps(Tsit5()), [_x(t - 0.1), 0.0]),
        (SDDEProblem, ImplicitEM(), [_x(t - 0.1) + a, b])
    ]
        @mtkcompile sys = System(
            [D(x) ~ x + rhss[1], p ~ x + y + rhss[2]], t; defaults = [p => missing], guesses = [
                x => 0.0, p => 0.0])
        prob = Problem(sys, [x => 1.0, y => 1.0], (0.0, 1.0))
        @test init(prob, alg).ps[p] ≈ 2.0
        # nonsensical value for y just to test that equations work
        prob2 = remake(prob; u0 = [x => 1.0, y => 2x + exp(x)])
        @test init(prob2, alg).ps[p] ≈ 3 + exp(1)
        # solve for `x` given `p` and `y`
        prob3 = remake(prob; u0 = [x => nothing, y => 1.0], p = [p => 2x + exp(y)])
        @test init(prob3, alg)[x] ≈ 1 - exp(1)
        @test_logs (:warn, r"overdetermined") remake(
            prob; u0 = [x => 1.0, y => 2.0], p = [p => 4.0])
        prob4 = remake(prob; u0 = [x => 1.0, y => 2.0], p = [p => 4.0])
        @test solve(prob4, alg).retcode == ReturnCode.InitialFailure
        prob5 = remake(prob)
        @test init(prob, alg).ps[p] ≈ 2.0
    end
end

@testset "`remake` changes initialization problem types" begin
    @variables _x(..) y(t) z(t)
    @parameters p q
    @brownians a
    x = _x(t)

    @testset "$Problem with $(SciMLBase.parameterless_type(typeof(alg)))" for (
        System, Problem, alg, rhss) in [
        (ModelingToolkit.System, ODEProblem, Tsit5(), 0),
        (ModelingToolkit.System, SDEProblem, ImplicitEM(), a),
        (ModelingToolkit.System, DDEProblem, MethodOfSteps(Tsit5()), _x(t - 0.1)),
        (ModelingToolkit.System, SDDEProblem, ImplicitEM(), _x(t - 0.1) + a)
    ]
        alge_eqs = [y^2 * q + q^2 * x ~ 0, z * p - p^2 * x * z ~ 0]

        @mtkcompile sys = System(
            [D(x) ~ x * p + y^2 * q + rhss; alge_eqs],
            t; guesses = [x => 0.0, y => 0.0, z => 0.0, p => 0.0, q => 0.0])
        prob = Problem(sys, [x => 1.0, p => 1.0, q => missing], (0.0, 1.0))
        @test is_variable(prob.f.initialization_data.initializeprob, q)
        ps = prob.p
        newps = SciMLStructures.replace(Tunable(), ps, ForwardDiff.Dual.(ps.tunable))
        prob2 = remake(prob; p = newps)
        @test eltype(state_values(prob2.f.initialization_data.initializeprob)) <:
              ForwardDiff.Dual
        @test eltype(prob2.f.initialization_data.initializeprob.p.tunable) <:
              ForwardDiff.Dual
        @test state_values(prob2.f.initialization_data.initializeprob) ≈
              state_values(prob.f.initialization_data.initializeprob)

        prob2 = remake(prob; u0 = ForwardDiff.Dual.(prob.u0))
        @test eltype(state_values(prob2.f.initialization_data.initializeprob)) <:
              ForwardDiff.Dual
        @test eltype(prob2.f.initialization_data.initializeprob.p.tunable) <:
              ForwardDiff.Dual
        @test state_values(prob2.f.initialization_data.initializeprob) ≈
              state_values(prob.f.initialization_data.initializeprob)

        prob2 = remake(prob; u0 = ForwardDiff.Dual.(prob.u0), p = newps)
        @test eltype(state_values(prob2.f.initialization_data.initializeprob)) <:
              ForwardDiff.Dual
        @test eltype(prob2.f.initialization_data.initializeprob.p.tunable) <:
              ForwardDiff.Dual
        @test state_values(prob2.f.initialization_data.initializeprob) ≈
              state_values(prob.f.initialization_data.initializeprob)

        prob2 = remake(prob; u0 = [x => ForwardDiff.Dual(1.0)],
            p = [p => ForwardDiff.Dual(1.0), q => missing])
        @test eltype(state_values(prob2.f.initialization_data.initializeprob)) <:
              ForwardDiff.Dual
        @test eltype(prob2.f.initialization_data.initializeprob.p.tunable) <:
              ForwardDiff.Dual
        @test state_values(prob2.f.initialization_data.initializeprob) ≈
              state_values(prob.f.initialization_data.initializeprob)
        @test eltype(prob2.p.initials) <: ForwardDiff.Dual
    end
end

@testset "`remake` preserves old u0map and pmap" begin
    @variables _x(..) y(t)
    @parameters p
    @brownians a
    x = _x(t)

    @testset "$Problem with $(SciMLBase.parameterless_type(typeof(alg)))" for (
        System, Problem, alg, rhss) in [
        (ModelingToolkit.System, ODEProblem, Tsit5(), 0),
        (ModelingToolkit.System, SDEProblem, ImplicitEM(), a),
        (ModelingToolkit.System, DDEProblem, MethodOfSteps(Tsit5()), _x(t - 0.1)),
        (ModelingToolkit.System, SDDEProblem, ImplicitEM(), _x(t - 0.1) + a)
    ]
        alge_eqs = [y^2 + 4y * p^2 ~ x^3]
        @mtkcompile sys = System(
            [D(x) ~ x + p * y^2 + rhss; alge_eqs], t; guesses = [
                y => 1.0, p => 1.0])
        prob = Problem(sys, [x => 1.0, p => 1.0], (0.0, 1.0))
        @test is_variable(prob.f.initialization_data.initializeprob, y)
        prob2 = @test_nowarn remake(prob; p = [p => 3.0]) # ensure no over/under-determined warning
        @test is_variable(prob.f.initialization_data.initializeprob, y)

        prob = Problem(sys, [y => 1.0, x => 2.0, p => missing], (0.0, 1.0))
        @test is_variable(prob.f.initialization_data.initializeprob, p)
        prob2 = @test_nowarn remake(prob; u0 = [y => 0.5])
        @test is_variable(prob.f.initialization_data.initializeprob, p)
    end
end

struct Multiplier{T}
    a::T
    b::T
end

function (m::Multiplier)(x, y)
    m.a * x + m.b * y
end

@register_symbolic Multiplier(x::Real, y::Real)

@testset "Nonnumeric parameter dependencies are retained" begin
    @variables x(t) y(t)
    @parameters foo(::Real, ::Real) p
    @mtkcompile sys = System([D(x) ~ t, 0 ~ foo(x, y), foo ~ Multiplier(p, 2p)], t;
        guesses = [y => -1.0])
    prob = ODEProblem(sys, [x => 1.0, p => 1.0], (0.0, 1.0))
    integ = init(prob, Rosenbrock23())
    @test integ[y] ≈ -0.5
end

@testset "Use observed equations for guesses of observed variables" begin
    @variables x(t) y(t) [state_priority = 100]
    @mtkcompile sys = System(
        [D(x) ~ x + t, y ~ 2x + 1], t; initialization_eqs = [x^3 + y^3 ~ 1])
    isys = ModelingToolkit.generate_initializesystem(sys)
    @test isequal(defaults(isys)[y], 2x + 1)
end

@testset "Create initializeprob when unknown has dependent value" begin
    @variables x(t) y(t)
    @mtkcompile sys = System([D(x) ~ x, D(y) ~ t * y], t; defaults = [x => 2y])
    prob = ODEProblem(sys, [y => 1.0], (0.0, 1.0))
    @test prob.f.initializeprob !== nothing
    integ = init(prob)
    @test integ[x] ≈ 2.0

    @variables x(t)[1:2] y(t)
    @mtkcompile sys = System([D(x) ~ x, D(y) ~ t], t; defaults = [x => [y, 3.0]])
    prob = ODEProblem(sys, [y => 1.0], (0.0, 1.0))
    @test prob.f.initializeprob !== nothing
    integ = init(prob)
    @test integ[x] ≈ [1.0, 3.0]
end

@testset "units" begin
    t = ModelingToolkit.t
    D = ModelingToolkit.D
    @parameters g [unit = u"m/s^2"] L=1 [unit = u"m^2"]
    @variables x(t) [unit = u"m"] y(t) [unit = u"m" state_priority = 10] λ(t) [unit = u"s^-2"]
    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ L]
    @mtkcompile pend = System(eqs, t)

    prob = ODEProblem(pend, [x => 1, y => 0, g => 1], (0.0, 1.5),
        guesses = ModelingToolkit.missing_variable_defaults(pend))
    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)

    prob2 = remake(prob, u0 = [x => 0.5, y=>nothing])
    sol2 = solve(prob2, Rodas5P())
    @test SciMLBase.successful_retcode(sol2)
end

@testset "Issue#3205" begin
    using ModelingToolkitStandardLibrary.Electrical
    import ModelingToolkitStandardLibrary.Mechanical.Rotational as MR
    using ModelingToolkitStandardLibrary.Blocks
    using SciMLBase

    function dc_motor(R1 = 0.5)
        R = R1 # [Ohm] armature resistance
        L = 4.5e-3 # [H] armature inductance
        k = 0.5 # [N.m/A] motor constant
        J = 0.02 # [kg.m²] inertia
        f = 0.01 # [N.m.s/rad] friction factor
        tau_L_step = -0.3 # [N.m] amplitude of the load torque step

        @named ground = Ground()
        @named source = Voltage()
        @named ref = Blocks.Step(height = 0.2, start_time = 0)
        @named pi_controller = Blocks.LimPI(k = 1.1, T = 0.035, u_max = 10, Ta = 0.035)
        @named feedback = Blocks.Feedback()
        @named R1 = Resistor(R = R)
        @named L1 = Inductor(L = L)
        @named emf = EMF(k = k)
        @named fixed = MR.Fixed()
        @named load = MR.Torque()
        @named load_step = Blocks.Step(height = tau_L_step, start_time = 3)
        @named inertia = MR.Inertia(J = J)
        @named friction = MR.Damper(d = f)
        @named speed_sensor = MR.SpeedSensor()

        connections = [connect(fixed.flange, emf.support, friction.flange_b)
                       connect(emf.flange, friction.flange_a, inertia.flange_a)
                       connect(inertia.flange_b, load.flange)
                       connect(inertia.flange_b, speed_sensor.flange)
                       connect(load_step.output, load.tau)
                       connect(ref.output, feedback.input1)
                       connect(speed_sensor.w, :y, feedback.input2)
                       connect(feedback.output, pi_controller.err_input)
                       connect(pi_controller.ctr_output, :u, source.V)
                       connect(source.p, R1.p)
                       connect(R1.n, L1.p)
                       connect(L1.n, emf.p)
                       connect(emf.n, source.n, ground.g)]

        @named model = System(connections, t,
            systems = [
                ground,
                ref,
                pi_controller,
                feedback,
                source,
                R1,
                L1,
                emf,
                fixed,
                load,
                load_step,
                inertia,
                friction,
                speed_sensor
            ])
    end

    model = dc_motor()
    sys = mtkcompile(model)

    prob = ODEProblem(sys, [sys.L1.i => 0.0], (0, 6.0))

    @test_nowarn remake(prob, p = prob.p)
end

@testset "Singular initialization prints a warning" begin
    @parameters g
    @variables x(t) y(t) [state_priority = 10] λ(t)
    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ 1]
    @mtkcompile pend = System(eqs, t)
    @test_warn ["structurally singular", "initialization", "Guess", "heuristic"] ODEProblem(
        pend, [x => 1, y => 0, g => 1], (0.0, 1.5), guesses = [λ => 1])
end

@testset "DAEProblem initialization" begin
    @variables x(t) [guess = 1.0] y(t) [guess = 1.0]
    @parameters p=missing [guess=1.0] q=missing [guess=1.0]
    @mtkcompile sys = System(
        [D(x) ~ p * y + q, x^3 + y^3 ~ 5], t; initialization_eqs = [p^2 + q^3 ~ 3])

    # FIXME: solve for du0
    prob = DAEProblem(
        sys, [D(x) => cbrt(4) + cbrt(2), D(y) => -1 / cbrt(4), x => 1.0, p => 1.0], (
            0.0, 1.0))

    integ = init(prob, DImplicitEuler())
    @test integ[x] ≈ 1.0
    @test integ[y]≈cbrt(4) rtol=1e-6
    @test integ.ps[p] ≈ 1.0
    @test integ.ps[q]≈cbrt(2) rtol=1e-6
end

@testset "Guesses provided to `ODEProblem` are used in `remake`" begin
    @variables x(t) y(t)=2x
    @parameters p q=3x
    @mtkcompile sys = System([D(x) ~ x * p + q, x^3 + y^3 ~ 3], t)
    prob = ODEProblem(
        sys, [p => 1.0], (0.0, 1.0); guesses = [x => 1.0, y => 1.0, q => 1.0])
    @test prob[x] == 1.0
    @test prob[y] == 2.0
    @test prob.ps[p] == 1.0
    @test prob.ps[q] == 3.0
    integ = init(prob)
    @test integ[x] ≈ 1 / cbrt(3)
    @test integ[y] ≈ 2 / cbrt(3)
    @test integ.ps[p] == 1.0
    @test integ.ps[q]≈3 / cbrt(3) atol=1e-5
    prob2 = remake(prob; u0 = [y => 3x], p = [q => 2x])
    integ2 = init(prob2)
    @test integ2[x]≈cbrt(3 / 28) atol=1e-5
    @test integ2[y]≈3cbrt(3 / 28) atol=1e-5
    @test integ2.ps[p] == 1.0
    @test integ2.ps[q]≈2cbrt(3 / 28) atol=1e-5
end

function test_dummy_initialization_equation(prob, var)
    initsys = prob.f.initialization_data.initializeprob.f.sys
    @test isempty(equations(initsys))
    idx = findfirst(eq -> isequal(var, eq.lhs), observed(initsys))
    @test idx !== nothing && is_parameter(initsys, observed(initsys)[idx].rhs)
end

@testset "Remake problem with no initializeprob" begin
    @variables x(t) [guess = 1.0] y(t) [guess = 1.0]
    @parameters p [guess = 1.0] q [guess = 1.0]
    @mtkcompile sys = System(
        [D(x) ~ p * x + q * y, y ~ 2x, q ~ 2p], t)
    prob = ODEProblem(sys, [x => 1.0, p => 1.0], (0.0, 1.0))
    test_dummy_initialization_equation(prob, x)
    prob2 = remake(prob; u0 = [x => 2.0])
    @test prob2[x] == 2.0
    test_dummy_initialization_equation(prob2, x)
    # otherwise we have `x ~ 2, y ~ 2` which is unsatisfiable
    prob3 = remake(prob; u0 = [x => nothing, y => 2.0])
    @test prob3.f.initialization_data !== nothing
    @test init(prob3)[x] ≈ 1.0
    prob4 = remake(prob; p = [p => 1.0])
    test_dummy_initialization_equation(prob4, x)
    prob5 = remake(prob; p = [p => missing, q => 4.0])
    @test prob5.f.initialization_data !== nothing
    @test init(prob5).ps[p] ≈ 2.0
end

@testset "Variables provided as symbols" begin
    @variables x(t) [guess = 1.0] y(t) [guess = 1.0]
    @parameters p [guess = 1.0] q [guess = 1.0]
    @mtkcompile sys = System(
        [D(x) ~ p * x + q * y, y ~ 2x, q ~ 2p], t)
    prob = ODEProblem(sys, [:x => 1.0, p => 1.0], (0.0, 1.0))
    test_dummy_initialization_equation(prob, x)
    prob2 = remake(prob; u0 = [:x => 2.0])
    test_dummy_initialization_equation(prob2, x)
    prob3 = remake(prob; u0 = [:y => 1.0, :x => nothing])
    @test init(prob3)[x] ≈ 0.5
    @test SciMLBase.successful_retcode(solve(prob3))
end

@testset "Issue#3246: type promotion with parameter dependent initialization_eqs" begin
    @variables x(t)=1 y(t)=1
    @parameters a = 1
    @named sys = System([D(x) ~ 0, D(y) ~ x + a], t; initialization_eqs = [y ~ a])

    ssys = mtkcompile(sys)
    prob = ODEProblem(ssys, [], (0, 1))

    @test SciMLBase.successful_retcode(solve(prob))

    seta = setsym_oop(prob, [a])
    (newu0, newp) = seta(prob, ForwardDiff.Dual{ForwardDiff.Tag{:tag, Float64}}.([1.0], 0))
    newprob = remake(prob, u0 = newu0, p = newp)

    @test SciMLBase.successful_retcode(solve(newprob))
end

@testset "Issue#3295: Incomplete initialization of pure-ODE systems" begin
    @variables X(t) Y(t)
    @parameters p d
    eqs = [
        D(X) ~ p - d * X,
        D(Y) ~ p - d * Y
    ]
    @mtkcompile osys = System(eqs, t)

    # Make problem.
    u0_vals = [X => 4, Y => 5.0]
    tspan = (0.0, 10.0)
    p_vals = [p => 1.0, d => 0.1]
    oprob = ODEProblem(osys, [u0_vals; p_vals], tspan)
    integ = init(oprob)
    @test integ[X] ≈ 4.0
    @test integ[Y] ≈ 5.0
    # Attempt to `remake`.
    rp = remake(oprob; u0 = [Y => 7])
    integ = init(rp)
    @test integ[X] ≈ 4.0
    @test integ[Y] ≈ 7.0
end

@testset "Issue#3297: `generate_initializesystem(::JumpSystem)`" begin
    @parameters β γ S0
    @variables S(t)=S0 I(t) R(t)
    rate₁ = β * S * I
    affect₁ = [S ~ Pre(S) - 1, I ~ Pre(I) + 1]
    rate₂ = γ * I
    affect₂ = [I ~ Pre(I) - 1, R ~ Pre(R) + 1]
    j₁ = ConstantRateJump(rate₁, affect₁)
    j₂ = ConstantRateJump(rate₂, affect₂)
    j₃ = MassActionJump(2 * β + γ, [R => 1], [S => 1, R => -1])
    @mtkcompile js = JumpSystem([j₁, j₂, j₃], t, [S, I, R], [β, γ, S0])

    u0s = [I => 1, R => 0]
    ps = [S0 => 999, β => 0.01, γ => 0.001]

    jprob = JumpProblem(js, [u0s; ps], (0.0, 10.0))
    sol = solve(jprob, SSAStepper())
    @test sol[S, 1] ≈ 999
    @test SciMLBase.successful_retcode(sol)
end

@testset "Solvable array parameters with scalarized guesses" begin
    @variables x(t)
    @parameters p[1:2] q
    @mtkcompile sys = System(
        D(x) ~ p[1] + p[2] + q, t; defaults = [p[1] => q, p[2] => 2q],
        guesses = [p[1] => q, p[2] => 2q])
    @test ModelingToolkit.is_parameter_solvable(p, Dict(), defaults(sys), guesses(sys))
    prob = ODEProblem(sys, [x => 1.0, q => 2.0], (0.0, 1.0))
    initsys = prob.f.initialization_data.initializeprob.f.sys
    @test length(ModelingToolkit.observed(initsys)) == 4
    sol = solve(prob, Tsit5())
    @test sol.ps[p] ≈ [2.0, 4.0]
end

@testset "Issue#3318: Mutating `Initial` parameters works" begin
    @variables x(t) y(t)[1:2] [guess = ones(2)]
    @parameters p[1:2, 1:2]
    @mtkcompile sys = System(
        [D(x) ~ x, D(y) ~ p * y], t; initialization_eqs = [x^2 + y[1]^2 + y[2]^2 ~ 4])
    prob = ODEProblem(sys, [x => 1.0, y[1] => 1, p => 2ones(2, 2)], (0.0, 1.0))
    integ = init(prob, Tsit5())
    @test integ[x] ≈ 1.0
    @test integ[y] ≈ [1.0, sqrt(2.0)]
    prob.ps[Initial(x)] = 0.5
    integ = init(prob, Tsit5())
    @test integ[x] ≈ 0.5
    @test integ[y] ≈ [1.0, sqrt(2.75)]
    prob.ps[Initial(y[1])] = 0.5
    integ = init(prob, Tsit5())
    @test integ[x] ≈ 0.5
    @test integ[y]≈[0.5, sqrt(3.5)] atol=1e-6
end

@testset "Issue#3342" begin
    @variables x(t) y(t)
    stop!(mod, obs, ctx, integrator) = (terminate!(integrator); return (;))
    @named sys = System([D(x) ~ 1.0
                         D(y) ~ 1.0], t; initialization_eqs = [
            y ~ 0.0
        ],
        continuous_events = [
            [y ~ 0.5] => (; f = stop!)
        ])
    sys = mtkcompile(sys)
    prob0 = ODEProblem(sys, [x => NaN], (0.0, 1.0))

    # final_x(x0) is equivalent to x0 + 0.5
    function final_x(x0)
        prob = remake(prob0; u0 = [x => x0])
        sol = solve(prob)
        return sol[x][end]
    end
    @test final_x(0.3) ≈ 0.8 # should be 0.8
    @test ForwardDiff.derivative(final_x, 0.3) ≈ 1.0
end

@testset "Issue#3330: Initialization for unsimplified systems" begin
    @variables x(t) [guess = 1.0]
    @mtkcompile sys = System(D(x) ~ x, t; initialization_eqs = [x^2 ~ 4])
    prob = ODEProblem(sys, [], (0.0, 1.0))
    @test prob.f.initialization_data !== nothing
end

@testset "`ReconstructInitializeprob` with `nothing` state" begin
    @parameters p
    @variables x(t)
    @mtkcompile sys = System(x ~ p * t, t)
    prob = @test_nowarn ODEProblem(sys, [p => 1.0], (0.0, 1.0))
    @test_nowarn remake(prob, p = [p => 1.0])
    @test_nowarn remake(prob, p = [p => ForwardDiff.Dual(1.0)])
end

@testset "`late_binding_update_u0_p` copies `newp`" begin
    @parameters k1 k2
    @variables X1(t) X2(t)
    @parameters Γ[1:1]=missing [guess = [1.0]]
    eqs = [
        D(X1) ~ k1 * (Γ[1] - X1) - k2 * X1
    ]
    obs = [X2 ~ Γ[1] - X1]
    @mtkcompile osys = System(eqs, t, [X1, X2], [k1, k2, Γ]; observed = obs)
    u0 = [X1 => 1.0, X2 => 2.0]
    ps = [k1 => 0.1, k2 => 0.2]

    oprob1 = ODEProblem(osys, [u0; ps], 1.0)
    oprob2 = remake(oprob1, u0 = [X1 => 10.0])
    integ1 = init(oprob1)
    @test integ1[X1] ≈ 1.0
end

@testset "Trivial initialization is run on problem construction" begin
    @variables _x(..) y(t)
    @brownians a
    @parameters tot
    x = _x(t)
    @testset "$Problem" for (Problem, lhs, rhs) in [
        (ODEProblem, D, 0.0),
        (SDEProblem, D, a),
        (DDEProblem, D, _x(t - 0.1)),
        (SDDEProblem, D, _x(t - 0.1) + a)
    ]
        @mtkcompile sys = ModelingToolkit.System([lhs(x) ~ x + rhs, x + y ~ tot], t;
            guesses = [tot => 1.0], defaults = [tot => missing])
        prob = Problem(sys, [x => 1.0, y => 1.0], (0.0, 1.0))
        @test prob.ps[tot] ≈ 2.0
    end
    @testset "$Problem" for Problem in [NonlinearProblem, NonlinearLeastSquaresProblem]
        @parameters p1 p2
        @mtkcompile sys = System([x^2 + y^2 ~ p1, (x - 1)^2 + (y - 1)^2 ~ p2, p2 ~ 2p1];
            guesses = [p1 => 0.0], defaults = [p1 => missing])
        prob = Problem(sys, [x => 1.0, y => 1.0, p2 => 6.0])
        @test prob.ps[p1] ≈ 3.0
    end
end

@testset "`Initial(X)` in time-independent systems: $Problem" for Problem in [
    NonlinearProblem, NonlinearLeastSquaresProblem]
    @parameters k1 k2
    @variables X1(t) X2(t)
    @parameters Γ[1:1]=missing [guess = [1.0]]
    eqs = [
        0 ~ k1 * (Γ[1] - X1) - k2 * X1
    ]
    initialization_eqs = [
        X2 ~ Γ[1] - X1
    ]
    @mtkcompile nlsys = System(eqs, [X1, X2], [k1, k2, Γ]; initialization_eqs)

    @testset "throws if initialization_eqs contain unknowns" begin
        u0 = [X1 => 1.0, X2 => 2.0]
        ps = [k1 => 0.1, k2 => 0.2]
        @test_throws ArgumentError Problem(nlsys, [u0; ps])
    end

    eqs = [0 ~ k1 * (Γ[1] - X1) - k2 * X1
           X2 ~ Γ[1] - X1]
    initialization_eqs = [
        Initial(X2) ~ Γ[1] - Initial(X1)
    ]
    @mtkcompile nlsys = System(eqs, [X1, X2], [k1, k2, Γ]; initialization_eqs)

    @testset "solves initialization" begin
        u0 = [X1 => 1.0, X2 => 2.0]
        ps = [k1 => 0.1, k2 => 0.2]
        prob = Problem(nlsys, [u0; ps])
        @test state_values(prob.f.initialization_data.initializeprob) === nothing
        @test prob.ps[Γ[1]] ≈ 3.0
    end

    @testset "respects explicitly provided value" begin
        ps = [k1 => 0.1, k2 => 0.2, Γ => [5.0]]
        prob = Problem(nlsys, ps)
        @test prob.ps[Γ[1]] ≈ 5.0
    end

    @testset "fails initialization if inconsistent explicit value" begin
        u0 = [X1 => 1.0, X2 => 2.0]
        ps = [k1 => 0.1, k2 => 0.2, Γ => [5.0]]
        prob = Problem(nlsys, [u0; ps])
        sol = solve(prob)
        @test sol.retcode == SciMLBase.ReturnCode.InitialFailure
    end

    @testset "Ignores initial equation if given insufficient u0" begin
        u0 = [X2 => 2.0]
        ps = [k1 => 0.1, k2 => 0.2, Γ => [5.0]]
        prob = Problem(nlsys, [u0; ps])
        sol = solve(prob)
        @test SciMLBase.successful_retcode(sol)
        @test sol.ps[Γ[1]] ≈ 5.0
    end
end

@testset "Issue#3504: Update initials when `remake` called with non-symbolic `u0`" begin
    @variables x(t) y(t)
    @parameters c1 c2
    @mtkcompile sys = System([D(x) ~ -c1 * x + c2 * y, D(y) ~ c1 * x - c2 * y], t)
    prob1 = ODEProblem(sys, [x => 1.0, y => 2.0, c1 => 1.0, c2 => 2.0], (0.0, 1.0))
    prob2 = remake(prob1, u0 = [2.0, 3.0])
    prob3 = remake(prob1, u0 = [2.0, 3.0], p = [c1 => 2.0])
    integ1 = init(prob1, Tsit5())
    integ2 = init(prob2, Tsit5())
    integ3 = init(prob3, Tsit5())
    @test integ2.u ≈ [2.0, 3.0]
    @test integ3.u ≈ [2.0, 3.0]
    @test integ3.ps[c1] ≈ 2.0
end

# https://github.com/SciML/SciMLBase.jl/issues/985
@testset "Type-stability of `remake`" begin
    @parameters α=1 β=1 γ=1 δ=1
    @variables x(t)=1 y(t)=1
    eqs = [D(x) ~ α * x - β * x * y, D(y) ~ -δ * y + γ * x * y]
    @named sys = System(eqs, t)
    prob = ODEProblem(complete(sys), [], (0.0, 1))
    @inferred remake(prob; u0 = 2 .* prob.u0, p = prob.p)
    @inferred solve(prob)
end

@testset "Issue#3570, #3552: `Initial`s/guesses are copied to `u0` during `solve`/`init`" begin
    @parameters g
    @variables x(t) [state_priority = 10] y(t) λ(t)
    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ 1]
    @mtkcompile pend = System(eqs, t)

    prob = ODEProblem(
        pend, [x => (√2 / 2), D(x) => 0.0, g => 1], (0.0, 1.5),
        guesses = [λ => 1, y => √2 / 2])
    sol = solve(prob)

    @testset "Guesses of initialization problem copied to algebraic variables" begin
        prob.f.initialization_data.initializeprob[λ] = 1.0
        prob2 = DiffEqBase.get_updated_symbolic_problem(
            pend, prob; u0 = prob.u0, p = prob.p)
        @test prob2[λ] ≈ 1.0
    end

    @testset "Initial values for algebraic variables are retained" begin
        prob2 = ODEProblem(
            pend, [x => (√2 / 2), D(y) => 0.0, g => 1], (0.0, 1.5),
            guesses = [λ => 1, y => √2 / 2])
        sol = solve(prob)
        @test SciMLBase.successful_retcode(sol)
        prob3 = DiffEqBase.get_updated_symbolic_problem(
            pend, prob2; u0 = prob2.u0, p = prob2.p)
        @test prob3[D(y)] ≈ 0.0
    end

    @testset "`setsym_oop`" begin
        setter = setsym_oop(prob, [Initial(x)])
        (u0, p) = setter(prob, [0.8])
        new_prob = remake(prob; u0, p, initializealg = BrownFullBasicInit())
        new_sol = solve(new_prob)
        @test new_sol[x, 1] ≈ 0.8
        integ = init(new_prob)
        @test integ[x] ≈ 0.8
    end

    @testset "`setsym`" begin
        @test prob.ps[Initial(x)] ≈ √2 / 2
        prob.ps[Initial(x)] = 0.8
        sol = solve(prob; initializealg = BrownFullBasicInit())
        @test sol[x, 1] ≈ 0.8
        integ = init(prob; initializealg = BrownFullBasicInit())
        @test integ[x] ≈ 0.8
    end
end

@testset "Initialization copies solved `u0` to `p`" begin
    @parameters σ ρ β A[1:3]
    @variables x(t) y(t) z(t) w(t) w2(t)
    eqs = [D(D(x)) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z,
        w ~ x + y + z + 2 * β,
        0 ~ x^2 + y^2 - w2^2
    ]

    @mtkcompile sys = System(eqs, t)

    u0 = [D(x) => 2.0,
        x => 1.0,
        y => 0.0,
        z => 0.0]

    p = [σ => 28.0,
        ρ => 10.0,
        β => 8 / 3]

    tspan = (0.0, 100.0)
    getter = getsym(sys, Initial.(unknowns(sys)))
    prob = ODEProblem(sys, [u0; p], tspan; guesses = [w2 => 3.0])
    new_u0, new_p,
    _ = SciMLBase.get_initial_values(
        prob, prob, prob.f, SciMLBase.OverrideInit(), Val(true);
        nlsolve_alg = NewtonRaphson(), abstol = 1e-6, reltol = 1e-3)
    @test getter(prob) != getter(new_p)
    @test getter(new_p) == new_u0
    _prob = remake(prob, u0 = new_u0, p = new_p)
    sol = solve(_prob; initializealg = CheckInit())
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[1] ≈ new_u0
end

@testset "Initialization system retains `split` kwarg of parent" begin
    @parameters g
    @variables x(t) y(t) [state_priority = 10] λ(t)
    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ 1]
    @mtkcompile pend=System(eqs, t) split=false
    prob = ODEProblem(
        pend, [x => 1.0, D(x) => 0.0, g => 1.0], (0.0, 1.0); guesses = [y => 1.0, λ => 1.0])
    @test !ModelingToolkit.is_split(prob.f.initialization_data.initializeprob.f.sys)
end

@testset "`InitializationProblem` retains `iip` of parent" begin
    @parameters g
    @variables x(t) y(t) [state_priority = 10] λ(t)
    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ 1]
    @mtkcompile pend = System(eqs, t)
    prob = ODEProblem(pend, SA[x => 1.0, D(x) => 0.0, g => 1.0],
        (0.0, 1.0); guesses = [y => 1.0, λ => 1.0])
    @test !SciMLBase.isinplace(prob)
    @test !SciMLBase.isinplace(prob.f.initialization_data.initializeprob)
end

@testset "Array unknowns occurring unscalarized in initializeprobpmap" begin
    @variables begin
        u(t)[1:2] = 0.9ones(2)
        x(t)[1:2], [guess = 0.01ones(2)]
        o(t)[1:2]
    end
    @parameters p[1:4] = [2.0, 1.875, 2.0, 1.875]

    eqs = [D(u[1]) ~ p[1] * u[1] - p[2] * u[1] * u[2] + x[1] + 0.1
           D(u[2]) ~ p[4] * u[1] * u[2] - p[3] * u[2] - x[2]
           o[1] ~ sum(p) * sum(u)
           o[2] ~ sum(p) * sum(x)
           x[1] ~ 0.01exp(-1)
           x[2] ~ 0.01cos(t)]

    @mtkcompile sys = System(eqs, t)
    prob = ODEProblem(sys, [], (0.0, 1.0))
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
end

@testset "Defaults removed with ` => nothing` aren't retained" begin
    @variables x(t)[1:2]
    @mtkcompile sys = System([D(x[1]) ~ -x[1], x[1] + x[2] ~ 3], t; defaults = [x[1] => 1])
    prob = ODEProblem(sys, [x[1] => nothing, x[2] => 1], (0.0, 1.0))
    @test SciMLBase.initialization_status(prob) == SciMLBase.FULLY_DETERMINED
end

@testset "Nonnumerics aren't narrowed" begin
    @mtkmodel Foo begin
        @variables begin
            x(t) = 1.0
        end
        @parameters begin
            p::AbstractString
            r = 1.0
        end
        @equations begin
            D(x) ~ r * x
        end
    end
    @mtkcompile sys = Foo(p = "a")
    prob = ODEProblem(sys, [], (0.0, 1.0))
    @test prob.p.nonnumeric[1] isa Vector{AbstractString}
    integ = init(prob)
    @test integ.p.nonnumeric[1] isa Vector{AbstractString}
end
