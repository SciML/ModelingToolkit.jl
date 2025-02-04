using ModelingToolkit, OrdinaryDiffEq, NonlinearSolve, Test
using StochasticDiffEq, DelayDiffEq, StochasticDelayDiffEq, JumpProcesses
using ForwardDiff
using SymbolicIndexingInterface, SciMLStructures
using SciMLStructures: Tunable
using ModelingToolkit: t_nounits as t, D_nounits as D
using DynamicQuantities

@parameters g
@variables x(t) y(t) [state_priority = 10] λ(t)
eqs = [D(D(x)) ~ λ * x
       D(D(y)) ~ λ * y - g
       x^2 + y^2 ~ 1]
@mtkbuild pend = ODESystem(eqs, t)

initprob = ModelingToolkit.InitializationProblem(pend, 0.0, [], [g => 1];
    guesses = [ModelingToolkit.missing_variable_defaults(pend); x => 1; y => 0.2])
conditions = getfield.(equations(initprob.f.sys), :rhs)

@test initprob isa NonlinearLeastSquaresProblem
sol = solve(initprob)
@test SciMLBase.successful_retcode(sol)
@test maximum(abs.(sol[conditions])) < 1e-14

@test_throws ModelingToolkit.ExtraVariablesSystemException ModelingToolkit.InitializationProblem(
    pend, 0.0, [], [g => 1];
    guesses = [ModelingToolkit.missing_variable_defaults(pend); x => 1; y => 0.2],
    fully_determined = true)

initprob = ModelingToolkit.InitializationProblem(pend, 0.0, [x => 1, y => 0], [g => 1];
    guesses = ModelingToolkit.missing_variable_defaults(pend))
@test initprob isa NonlinearLeastSquaresProblem
sol = solve(initprob)
@test SciMLBase.successful_retcode(sol)
@test all(iszero, sol.u)
@test maximum(abs.(sol[conditions])) < 1e-14

initprob = ModelingToolkit.InitializationProblem(
    pend, 0.0, [], [g => 1]; guesses = ModelingToolkit.missing_variable_defaults(pend))
@test initprob isa NonlinearLeastSquaresProblem
sol = solve(initprob)
@test !SciMLBase.successful_retcode(sol)

@test_throws ModelingToolkit.ExtraVariablesSystemException ModelingToolkit.InitializationProblem(
    pend, 0.0, [], [g => 1]; guesses = ModelingToolkit.missing_variable_defaults(pend),
    fully_determined = true)

prob = ODEProblem(pend, [x => 1, y => 0], (0.0, 1.5), [g => 1],
    guesses = ModelingToolkit.missing_variable_defaults(pend))
prob.f.initializeprob isa NonlinearProblem
sol = solve(prob.f.initializeprob)
@test maximum(abs.(sol[conditions])) < 1e-14
sol = solve(prob, Rodas5P())
@test maximum(abs.(sol[conditions][1])) < 1e-14

prob = ODEProblem(pend, [x => 1], (0.0, 1.5), [g => 1],
    guesses = ModelingToolkit.missing_variable_defaults(pend))
prob.f.initializeprob isa NonlinearLeastSquaresProblem
sol = solve(prob.f.initializeprob)
@test maximum(abs.(sol[conditions])) < 1e-14
sol = solve(prob, Rodas5P())
@test maximum(abs.(sol[conditions][1])) < 1e-14

@test_throws ModelingToolkit.ExtraVariablesSystemException ODEProblem(
    pend, [x => 1], (0.0, 1.5), [g => 1],
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

@mtkmodel System begin
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

@mtkbuild sys = System()
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

@mtkbuild sys = MassDamperSystem()
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

@mtkbuild sys = ODESystem(eqs, t)

u0 = [D(x) => 2.0,
    y => 0.0,
    z => 0.0]

p = [σ => 28.0,
    ρ => 10.0,
    β => 8 / 3]

tspan = (0.0, 100.0)
@test_throws ModelingToolkit.IncompleteInitializationError prob=ODEProblem(
    sys, u0, tspan, p, jac = true)

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
    return ODESystem(eqs, t, vars, []; name)
end

@mtkbuild sys = System2()
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
    return ODESystem(eqs, t, vars, []; name, initialization_eqs)
end

@mtkbuild sys = System3()
prob = ODEProblem(sys, [], (0, 1), guesses = [sys.dx => 1])
sol = solve(prob, Tsit5())
@test SciMLBase.successful_retcode(sol)
@test sol.u[1] == [1.0]

# Steady state initialization

@parameters σ ρ β
@variables x(t) y(t) z(t)

eqs = [D(D(x)) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@named sys = ODESystem(eqs, t)
sys = structural_simplify(sys)

u0 = [D(x) => 2.0,
    x => 1.0,
    D(y) => 0.0,
    z => 0.0]

p = [σ => 28.0,
    ρ => 10.0,
    β => 8 / 3]

tspan = (0.0, 0.2)
prob_mtk = ODEProblem(sys, u0, tspan, p)
sol = solve(prob_mtk, Tsit5())
@test sol[x * (ρ - z) - y][1] == 0.0

@variables x(t) y(t) z(t)
@parameters α=1.5 β=1.0 γ=3.0 δ=1.0

eqs = [D(x) ~ α * x - β * x * y
       D(y) ~ -γ * y + δ * x * y
       z ~ x + y]

@named sys = ODESystem(eqs, t)
simpsys = structural_simplify(sys)
tspan = (0.0, 10.0)

prob = ODEProblem(simpsys, [D(x) => 0.0, y => 0.0], tspan, guesses = [x => 0.0])
sol = solve(prob, Tsit5())
@test sol.u[1] == [0.0, 0.0]

# Initialize with an observed variable
prob = ODEProblem(simpsys, [z => 0.0], tspan, guesses = [x => 2.0, y => 4.0])
sol = solve(prob, Tsit5())
@test sol.u[1] == [0.0, 0.0]

prob = ODEProblem(simpsys, [z => 1.0, y => 1.0], tspan, guesses = [x => 2.0])
sol = solve(prob, Tsit5())
@test sol.u[1] == [0.0, 1.0]

# This should warn, but logging tests can't be marked as broken
@test_logs prob = ODEProblem(simpsys, [], tspan, guesses = [x => 2.0])

# Late Binding initialization_eqs
# https://github.com/SciML/ModelingToolkit.jl/issues/2787

@parameters g
@variables x(t) y(t) [state_priority = 10] λ(t)
eqs = [D(D(x)) ~ λ * x
       D(D(y)) ~ λ * y - g
       x^2 + y^2 ~ 1]
@mtkbuild pend = ODESystem(eqs, t)

prob = ODEProblem(pend, [x => 1], (0.0, 1.5), [g => 1],
    guesses = [λ => 0, y => 1], initialization_eqs = [y ~ 1])

unsimp = generate_initializesystem(pend; u0map = [x => 1], initialization_eqs = [y ~ 1])
sys = structural_simplify(unsimp; fully_determined = false)
@test length(equations(sys)) == 3

# Extend two systems with initialization equations and guesses
# https://github.com/SciML/ModelingToolkit.jl/issues/2845
@variables x(t) y(t)
@named sysx = ODESystem([D(x) ~ 0], t; initialization_eqs = [x ~ 1])
@named sysy = ODESystem([D(y) ~ 0], t; initialization_eqs = [y^2 ~ 2], guesses = [y => 1])
sys = extend(sysx, sysy)
@test length(equations(generate_initializesystem(sys))) == 2
@test length(ModelingToolkit.guesses(sys)) == 1

# https://github.com/SciML/ModelingToolkit.jl/issues/2873
@testset "Error on missing defaults" begin
    @variables x(t) y(t)
    @named sys = ODESystem([x^2 + y^2 ~ 25, D(x) ~ 1], t)
    ssys = structural_simplify(sys)
    @test_throws ModelingToolkit.MissingVariablesError ODEProblem(
        ssys, [x => 3], (0, 1), []) # y should have a guess
end

# https://github.com/SciML/ModelingToolkit.jl/issues/3025
@testset "Override defaults when setting initial conditions with unknowns(sys) or similar" begin
    @variables x(t) y(t)

    # system 1 should solve to x = 1
    ics1 = [x => 1]
    sys1 = ODESystem([D(x) ~ 0], t; defaults = ics1, name = :sys1) |> structural_simplify
    prob1 = ODEProblem(sys1, [], (0.0, 1.0), [])
    sol1 = solve(prob1, Tsit5())
    @test all(sol1[x] .== 1)

    # system 2 should solve to x = y = 2
    sys2 = extend(
        sys1,
        ODESystem([D(y) ~ 0], t; initialization_eqs = [y ~ 2], name = :sys2)
    ) |> structural_simplify
    ics2 = unknowns(sys1) .=> 2 # should be equivalent to "ics2 = [x => 2]"
    prob2 = ODEProblem(sys2, ics2, (0.0, 1.0), []; fully_determined = true)
    sol2 = solve(prob2, Tsit5())
    @test all(sol2[x] .== 2) && all(sol2[y] .== 2)
end

# https://github.com/SciML/ModelingToolkit.jl/issues/3029
@testset "Derivatives in initialization equations" begin
    @variables x(t)
    sys = ODESystem(
        [D(D(x)) ~ 0], t; initialization_eqs = [x ~ 0, D(x) ~ 1], name = :sys) |>
          structural_simplify
    @test_nowarn ODEProblem(sys, [], (0.0, 1.0), [])

    sys = ODESystem(
        [D(D(x)) ~ 0], t; initialization_eqs = [x ~ 0, D(D(x)) ~ 0], name = :sys) |>
          structural_simplify
    @test_nowarn ODEProblem(sys, [D(x) => 1.0], (0.0, 1.0), [])
end

# https://github.com/SciML/ModelingToolkit.jl/issues/3049
@testset "Derivatives in initialization guesses" begin
    for sign in [-1.0, +1.0]
        @variables x(t)
        sys = ODESystem(
            [D(D(x)) ~ 0], t;
            initialization_eqs = [D(x)^2 ~ 1, x ~ 0], guesses = [D(x) => sign], name = :sys
        ) |> structural_simplify
        prob = ODEProblem(sys, [], (0.0, 1.0), [])
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
@mtkbuild sys_1st_order = ODESystem(eqs_1st_order, t)
@mtkbuild sys_2nd_order = ODESystem(eqs_2nd_order, t)

u0_1st_order_1 = [X => 1.0, Y => 2.0]
u0_1st_order_2 = [Y => 2.0]
u0_2nd_order_1 = [X => 1.0, Y => 2.0, D(Y) => 0.5]
u0_2nd_order_2 = [Y => 2.0, D(Y) => 0.5]
tspan = (0.0, 10.0)
ps = [ω => 0.5, k1 => 2.0, k2 => 3.0]

oprob_1st_order_1 = ODEProblem(sys_1st_order, u0_1st_order_1, tspan, ps)
oprob_1st_order_2 = ODEProblem(sys_1st_order, u0_1st_order_2, tspan, ps)
oprob_2nd_order_1 = ODEProblem(sys_2nd_order, u0_2nd_order_1, tspan, ps) # gives sys_2nd_order
oprob_2nd_order_2 = ODEProblem(sys_2nd_order, u0_2nd_order_2, tspan, ps)

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
    @named sys = ODESystem([D(x) ~ x, D(y) ~ y], t; initialization_eqs = [y ~ -x])
    sys = structural_simplify(sys)
    prob = ODEProblem(sys, [sys.x => ones(5)], (0.0, 1.0), [])
    sol = solve(prob, Tsit5(), reltol = 1e-4)
    @test all(sol(1.0, idxs = sys.x) .≈ +exp(1)) && all(sol(1.0, idxs = sys.y) .≈ -exp(1))
end

@testset "Initialization of parameters" begin
    @variables _x(..) y(t)
    @parameters p q
    @brownian a b
    x = _x(t)

    # `System` constructor creates appropriate type with mtkbuild
    # `Problem` and `alg` create the problem to test and allow calling `init` with
    # the correct solver.
    # `rhss` allows adding terms to the end of equations (only 2 equations allowed) to influence
    # the system type (brownian vars to turn it into an SDE).
    @testset "$Problem with $(SciMLBase.parameterless_type(alg))" for (System, Problem, alg, rhss) in [
        (ModelingToolkit.System, ODEProblem, Tsit5(), zeros(2)),
        (ModelingToolkit.System, SDEProblem, ImplicitEM(), [a, b]),
        (ModelingToolkit.System, DDEProblem, MethodOfSteps(Tsit5()), [_x(t - 0.1), 0.0]),
        (ModelingToolkit.System, SDDEProblem, ImplicitEM(), [_x(t - 0.1) + a, b])
    ]
        function test_parameter(prob, sym, val, initialval = zero(val))
            @test prob.ps[sym] ≈ initialval
            if prob.u0 !== nothing
                @test init(prob, alg).ps[sym] ≈ val
            end
            @test solve(prob, alg).ps[sym] ≈ val
        end
        function test_initializesystem(sys, u0map, pmap, p, equation)
            isys = ModelingToolkit.generate_initializesystem(
                sys; u0map, pmap, guesses = ModelingToolkit.guesses(sys))
            @test is_variable(isys, p)
            @test equation in equations(isys) || (0 ~ -equation.rhs) in equations(isys)
        end

        u0map = Dict(x => 1.0, y => 1.0)
        pmap = Dict()
        pmap[q] = 1.0
        # `missing` default, equation from Problem
        @mtkbuild sys = System(
            [D(x) ~ x * q + rhss[1], D(y) ~ y * p + rhss[2]], t; defaults = [p => missing], guesses = [p => 1.0])
        pmap[p] = 2q
        prob = Problem(sys, u0map, (0.0, 1.0), pmap)
        test_parameter(prob, p, 2.0)
        prob2 = remake(prob; u0 = u0map, p = pmap)
        prob2.ps[p] = 0.0
        test_parameter(prob2, p, 2.0)
        # `missing` default, provided guess
        @mtkbuild sys = System(
            [D(x) ~ x + rhss[1], p ~ x + y + rhss[2]], t; defaults = [p => missing], guesses = [p => 0.0])
        prob = Problem(sys, u0map, (0.0, 1.0))
        test_parameter(prob, p, 2.0)
        test_initializesystem(sys, u0map, pmap, p, 0 ~ p - x - y)
        prob2 = remake(prob; u0 = u0map)
        prob2.ps[p] = 0.0
        test_parameter(prob2, p, 2.0)

        # `missing` to Problem, equation from default
        @mtkbuild sys = System(
            [D(x) ~ x * q + rhss[1], D(y) ~ y * p + rhss[2]], t; defaults = [p => 2q], guesses = [p => 1.0])
        pmap[p] = missing
        prob = Problem(sys, u0map, (0.0, 1.0), pmap)
        test_parameter(prob, p, 2.0)
        test_initializesystem(sys, u0map, pmap, p, 0 ~ 2q - p)
        prob2 = remake(prob; u0 = u0map, p = pmap)
        prob2.ps[p] = 0.0
        test_parameter(prob2, p, 2.0)
        # `missing` to Problem, provided guess
        @mtkbuild sys = System(
            [D(x) ~ x + rhss[1], p ~ x + y + rhss[2]], t; guesses = [p => 0.0])
        prob = Problem(sys, u0map, (0.0, 1.0), pmap)
        test_parameter(prob, p, 2.0)
        test_initializesystem(sys, u0map, pmap, p, 0 ~ x + y - p)
        prob2 = remake(prob; u0 = u0map, p = pmap)
        prob2.ps[p] = 0.0
        test_parameter(prob2, p, 2.0)

        # No `missing`, default and guess
        @mtkbuild sys = System(
            [D(x) ~ x * q + rhss[1], D(y) ~ y * p + rhss[2]], t; defaults = [p => 2q], guesses = [p => 0.0])
        delete!(pmap, p)
        prob = Problem(sys, u0map, (0.0, 1.0), pmap)
        test_parameter(prob, p, 2.0)
        test_initializesystem(sys, u0map, pmap, p, 0 ~ 2q - p)
        prob2 = remake(prob; u0 = u0map, p = pmap)
        prob2.ps[p] = 0.0
        test_parameter(prob2, p, 2.0)

        # Default overridden by Problem, guess provided
        @mtkbuild sys = System(
            [D(x) ~ q * x + rhss[1], D(y) ~ y * p + rhss[2]], t; defaults = [p => 2q], guesses = [p => 1.0])
        _pmap = merge(pmap, Dict(p => q))
        prob = Problem(sys, u0map, (0.0, 1.0), _pmap)
        test_parameter(prob, p, _pmap[q])
        test_initializesystem(sys, u0map, _pmap, p, 0 ~ q - p)
        # Problem dependent value with guess, no `missing`
        @mtkbuild sys = System(
            [D(x) ~ y * q + p + rhss[1], D(y) ~ x * p + q + rhss[2]], t; guesses = [p => 0.0])
        _pmap = merge(pmap, Dict(p => 3q))
        prob = Problem(sys, u0map, (0.0, 1.0), _pmap)
        test_parameter(prob, p, 3pmap[q])

        # Should not be solved for:
        # Override dependent default with direct value
        @mtkbuild sys = System(
            [D(x) ~ q * x + rhss[1], D(y) ~ y * p + rhss[2]], t; defaults = [p => 2q], guesses = [p => 1.0])
        _pmap = merge(pmap, Dict(p => 1.0))
        prob = Problem(sys, u0map, (0.0, 1.0), _pmap)
        @test prob.ps[p] ≈ 1.0
        @test prob.f.initialization_data === nothing

        # Non-floating point
        @parameters r::Int s::Int
        @mtkbuild sys = System(
            [D(x) ~ s * x + rhss[1], D(y) ~ y * r + rhss[2]], t; defaults = [s => 2r], guesses = [s => 1.0])
        prob = Problem(sys, u0map, (0.0, 1.0), [r => 1])
        @test prob.ps[r] == 1
        @test prob.ps[s] == 2
        @test prob.f.initialization_data === nothing

        @mtkbuild sys = System(
            [D(x) ~ x + rhss[1], p ~ x + y + rhss[2]], t; guesses = [p => 0.0])
        @test_throws ModelingToolkit.MissingParametersError Problem(
            sys, [x => 1.0, y => 1.0], (0.0, 1.0))

        # Unsatisfiable initialization
        prob = Problem(sys, [x => 1.0, y => 1.0], (0.0, 1.0),
            [p => 2.0]; initialization_eqs = [x^2 + y^2 ~ 3])
        @test prob.f.initialization_data !== nothing
        @test solve(prob, alg).retcode == ReturnCode.InitialFailure
        cache = init(prob, alg)
        @test solve!(cache).retcode == ReturnCode.InitialFailure
    end

    @testset "Null system" begin
        @variables x(t) y(t) s(t)
        @parameters x0 y0
        @mtkbuild sys = ODESystem([x ~ x0, y ~ y0, s ~ x + y], t; guesses = [y0 => 0.0])
        prob = ODEProblem(sys, [s => 1.0], (0.0, 1.0), [x0 => 0.3, y0 => missing])
        @test prob.ps[y0] ≈ 0.0
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
    @mtkbuild sys = ODESystem(
        [connect(fixed.flange, spring.flange_a), connect(spring.flange_b, mass.flange_a),
            connect(mass.flange_a, gravity.flange), connect(constant.output, gravity.f),
            connect(fixed.flange, damper.flange_a), connect(damper.flange_b, mass.flange_a)],
        t;
        systems = [fixed, spring, mass, gravity, constant, damper],
        guesses = [spring.s_rel0 => 1.0])
    prob = ODEProblem(sys, [], (0.0, 1.0), [spring.s_rel0 => missing])
    @test prob.ps[spring.s_rel0] ≈ 0.0
    @test init(prob, Tsit5()).ps[spring.s_rel0] ≈ -3.905
    @test solve(prob, Tsit5()).ps[spring.s_rel0] ≈ -3.905
end

@testset "NonlinearSystem initialization" begin
    nl_algs = [FastShortcutNonlinearPolyalg(), NewtonRaphson(),
        Klement(), SimpleNewtonRaphson(), DFSane()]
    nlls_algs = [FastShortcutNLLSPolyalg(), LevenbergMarquardt(), SimpleGaussNewton()]

    @testset "No initialization for variables" begin
        @variables x=1.0 y=0.0 z=0.0
        @parameters σ=10.0 ρ=26.0 β=8 / 3

        eqs = [0 ~ σ * (y - x),
            0 ~ x * (ρ - z) - y,
            0 ~ x * y - β * z]
        @mtkbuild ns = NonlinearSystem(eqs, [x, y, z], [σ, ρ, β])

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
    @testset "Parameter initialization" begin
        function test_parameter(prob, alg, param, val)
            integ = init(prob, alg)
            @test integ.ps[param]≈val rtol=1e-6
            sol = solve(prob, alg)
            @test sol.ps[param]≈val rtol=1e-6
            @test SciMLBase.successful_retcode(sol)
        end
        @variables x=1.0 y=3.0
        @parameters p q

        @named sys = NonlinearSystem(
            [x + y - p ~ 0, x - q ~ 0]; defaults = [q => missing],
            guesses = [q => 1.0], initialization_eqs = [p^2 + q^2 + 2p * q ~ 0])
        sys = complete(sys)
        for (probT, algs) in prob_alg_combinations
            prob = probT(sys, [], [p => 2.0])
            @test prob.f.initialization_data !== nothing
            @test prob.f.initialization_data.initializeprobmap === nothing
            for alg in algs
                test_parameter(prob, alg, q, -2.0)
            end

            # `update_initializeprob!` works
            prob.ps[p] = -2.0
            for alg in algs
                test_parameter(prob, alg, q, 2.0)
            end
            prob.ps[p] = 2.0

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
    @brownian a b
    x = _x(t)

    @testset "$Problem with $(SciMLBase.parameterless_type(typeof(alg)))" for (System, Problem, alg, rhss) in [
        (ModelingToolkit.System, ODEProblem, Tsit5(), zeros(2)),
        (ModelingToolkit.System, SDEProblem, ImplicitEM(), [a, b]),
        (ModelingToolkit.System, DDEProblem, MethodOfSteps(Tsit5()), [_x(t - 0.1), 0.0]),
        (ModelingToolkit.System, SDDEProblem, ImplicitEM(), [_x(t - 0.1) + a, b])
    ]
        @mtkbuild sys = System(
            [D(x) ~ x + rhss[1], p ~ x + y + rhss[2]], t; guesses = [x => 0.0, p => 0.0])
        prob = Problem(sys, [y => 1.0], (0.0, 1.0), [p => 3.0])
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
    @brownian a
    x = _x(t)

    @testset "$Problem with $(SciMLBase.parameterless_type(typeof(alg)))" for (System, Problem, alg, rhss) in [
        (ModelingToolkit.System, ODEProblem, Tsit5(), 0),
        (ModelingToolkit.System, SDEProblem, ImplicitEM(), a),
        (ModelingToolkit.System, DDEProblem, MethodOfSteps(Tsit5()), _x(t - 0.1)),
        (ModelingToolkit.System, SDDEProblem, ImplicitEM(), _x(t - 0.1) + a)
    ]
        @mtkbuild sys = System(
            [D(x) ~ 2x + r + rhss], t; parameter_dependencies = [r ~ p + 2q, q ~ p + 3],
            guesses = [p => 1.0])
        prob = Problem(sys, [x => 1.0], (0.0, 1.0), [p => missing])
        @test length(equations(ModelingToolkit.get_parent(prob.f.initialization_data.initializeprob.f.sys))) ==
              4
        integ = init(prob, alg)
        @test integ.ps[p] ≈ 2
    end
end

@testset "Re-creating initialization problem on remake" begin
    @variables _x(..) y(t)
    @parameters p q
    @brownian a b
    x = _x(t)

    @testset "$Problem with $(SciMLBase.parameterless_type(typeof(alg)))" for (System, Problem, alg, rhss) in [
        (ModelingToolkit.System, ODEProblem, Tsit5(), zeros(2)),
        (ModelingToolkit.System, SDEProblem, ImplicitEM(), [a, b]),
        (ModelingToolkit.System, DDEProblem, MethodOfSteps(Tsit5()), [_x(t - 0.1), 0.0]),
        (ModelingToolkit.System, SDDEProblem, ImplicitEM(), [_x(t - 0.1) + a, b])
    ]
        @mtkbuild sys = System(
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
    @brownian a
    x = _x(t)

    @testset "$Problem with $(SciMLBase.parameterless_type(typeof(alg)))" for (System, Problem, alg, rhss) in [
        (ModelingToolkit.System, ODEProblem, Tsit5(), 0),
        (ModelingToolkit.System, SDEProblem, ImplicitEM(), a),
        (ModelingToolkit.System, DDEProblem, MethodOfSteps(Tsit5()), _x(t - 0.1)),
        (ModelingToolkit.System, SDDEProblem, ImplicitEM(), _x(t - 0.1) + a)
    ]
        alge_eqs = [y^2 * q + q^2 * x ~ 0, z * p - p^2 * x * z ~ 0]

        @mtkbuild sys = System(
            [D(x) ~ x * p + y^2 * q + rhss; alge_eqs],
            t; guesses = [x => 0.0, y => 0.0, z => 0.0, p => 0.0, q => 0.0])
        prob = Problem(sys, [x => 1.0], (0.0, 1.0), [p => 1.0, q => missing])
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
        @test eltype(prob2.f.initialization_data.initializeprob.p.tunable) <: Float64
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
    end
end

@testset "`remake` preserves old u0map and pmap" begin
    @variables _x(..) y(t)
    @parameters p
    @brownian a
    x = _x(t)

    @testset "$Problem with $(SciMLBase.parameterless_type(typeof(alg)))" for (System, Problem, alg, rhss) in [
        (ModelingToolkit.System, ODEProblem, Tsit5(), 0),
        (ModelingToolkit.System, SDEProblem, ImplicitEM(), a),
        (ModelingToolkit.System, DDEProblem, MethodOfSteps(Tsit5()), _x(t - 0.1)),
        (ModelingToolkit.System, SDDEProblem, ImplicitEM(), _x(t - 0.1) + a)
    ]
        alge_eqs = [y^2 + 4y * p^2 ~ x^3]
        @mtkbuild sys = System(
            [D(x) ~ x + p * y^2 + rhss; alge_eqs], t; guesses = [
                y => 1.0, p => 1.0])
        prob = Problem(sys, [x => 1.0], (0.0, 1.0), [p => 1.0])
        @test is_variable(prob.f.initialization_data.initializeprob, y)
        prob2 = @test_nowarn remake(prob; p = [p => 3.0]) # ensure no over/under-determined warning
        @test is_variable(prob.f.initialization_data.initializeprob, y)

        prob = Problem(sys, [y => 1.0, x => 2.0], (0.0, 1.0), [p => missing])
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
    @mtkbuild sys = ODESystem([D(x) ~ t, 0 ~ foo(x, y)], t;
        parameter_dependencies = [foo ~ Multiplier(p, 2p)], guesses = [y => -1.0])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0), [p => 1.0])
    integ = init(prob, Rosenbrock23())
    @test integ[y] ≈ -0.5
end

@testset "Use observed equations for guesses of observed variables" begin
    @variables x(t) y(t) [state_priority = 100]
    @mtkbuild sys = ODESystem(
        [D(x) ~ x + t, y ~ 2x + 1], t; initialization_eqs = [x^3 + y^3 ~ 1])
    isys = ModelingToolkit.generate_initializesystem(sys)
    @test isequal(defaults(isys)[y], 2x + 1)
end

@testset "Create initializeprob when unknown has dependent value" begin
    @variables x(t) y(t)
    @mtkbuild sys = ODESystem([D(x) ~ x, D(y) ~ t * y], t; defaults = [x => 2y])
    prob = ODEProblem(sys, [y => 1.0], (0.0, 1.0))
    @test prob.f.initializeprob !== nothing
    integ = init(prob)
    @test integ[x] ≈ 2.0

    @variables x(t)[1:2] y(t)
    @mtkbuild sys = ODESystem([D(x) ~ x, D(y) ~ t], t; defaults = [x => [y, 3.0]])
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
    @mtkbuild pend = ODESystem(eqs, t)

    prob = ODEProblem(pend, [x => 1, y => 0], (0.0, 1.5), [g => 1],
        guesses = ModelingToolkit.missing_variable_defaults(pend))
    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)
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

        @named model = ODESystem(connections, t,
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
    sys = structural_simplify(model)

    prob = ODEProblem(sys, [sys.L1.i => 0.0], (0, 6.0))

    @test_nowarn remake(prob, p = prob.p)
end

@testset "Singular initialization prints a warning" begin
    @parameters g
    @variables x(t) y(t) [state_priority = 10] λ(t)
    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ 1]
    @mtkbuild pend = ODESystem(eqs, t)
    @test_warn ["structurally singular", "initialization", "Guess", "heuristic"] ODEProblem(
        pend, [x => 1, y => 0], (0.0, 1.5), [g => 1], guesses = [λ => 1])
end

@testset "DAEProblem initialization" begin
    @variables x(t) [guess = 1.0] y(t) [guess = 1.0]
    @parameters p=missing [guess = 1.0] q=missing [guess = 1.0]
    @mtkbuild sys = ODESystem(
        [D(x) ~ p * y + q * t, x^3 + y^3 ~ 5], t; initialization_eqs = [p^2 + q^3 ~ 3])

    # FIXME: solve for du0
    prob = DAEProblem(
        sys, [D(x) => cbrt(4), D(y) => -1 / cbrt(4)], [x => 1.0], (0.0, 1.0), [p => 1.0])

    integ = init(prob, DImplicitEuler())
    @test integ[x] ≈ 1.0
    @test integ[y]≈cbrt(4) rtol=1e-6
    @test integ.ps[p] ≈ 1.0
    @test integ.ps[q]≈cbrt(2) rtol=1e-6
end

@testset "Guesses provided to `ODEProblem` are used in `remake`" begin
    @variables x(t) y(t)=2x
    @parameters p q=3x
    @mtkbuild sys = ODESystem([D(x) ~ x * p + q, x^3 + y^3 ~ 3], t)
    prob = ODEProblem(
        sys, [], (0.0, 1.0), [p => 1.0]; guesses = [x => 1.0, y => 1.0, q => 1.0])
    @test prob[x] == 0.0
    @test prob[y] == 0.0
    @test prob.ps[p] == 1.0
    @test prob.ps[q] == 0.0
    integ = init(prob)
    @test integ[x] ≈ 1 / cbrt(3)
    @test integ[y] ≈ 2 / cbrt(3)
    @test integ.ps[p] == 1.0
    @test integ.ps[q] ≈ 3 / cbrt(3)
    prob2 = remake(prob; u0 = [y => 3x], p = [q => 2x])
    integ2 = init(prob2)
    @test integ2[x] ≈ cbrt(3 / 28)
    @test integ2[y] ≈ 3cbrt(3 / 28)
    @test integ2.ps[p] == 1.0
    @test integ2.ps[q] ≈ 2cbrt(3 / 28)
end

@testset "Remake problem with no initializeprob" begin
    @variables x(t) [guess = 1.0] y(t) [guess = 1.0]
    @parameters p [guess = 1.0] q [guess = 1.0]
    @mtkbuild sys = ODESystem(
        [D(x) ~ p * x + q * y, y ~ 2x], t; parameter_dependencies = [q ~ 2p])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0), [p => 1.0])
    @test prob.f.initialization_data === nothing
    prob2 = remake(prob; u0 = [x => 2.0])
    @test prob2[x] == 2.0
    @test prob2.f.initialization_data === nothing
    prob3 = remake(prob; u0 = [y => 2.0])
    @test prob3.f.initialization_data !== nothing
    @test init(prob3)[x] ≈ 1.0
    prob4 = remake(prob; p = [p => 1.0])
    @test prob4.f.initialization_data === nothing
    prob5 = remake(prob; p = [p => missing, q => 2.0])
    @test prob5.f.initialization_data !== nothing
    @test init(prob5).ps[p] ≈ 1.0
end

@testset "Variables provided as symbols" begin
    @variables x(t) [guess = 1.0] y(t) [guess = 1.0]
    @parameters p [guess = 1.0] q [guess = 1.0]
    @mtkbuild sys = ODESystem(
        [D(x) ~ p * x + q * y, y ~ 2x], t; parameter_dependencies = [q ~ 2p])
    prob = ODEProblem(sys, [:x => 1.0], (0.0, 1.0), [p => 1.0])
    @test prob.f.initialization_data === nothing
    prob2 = remake(prob; u0 = [:x => 2.0])
    @test prob2.f.initialization_data === nothing
    prob3 = remake(prob; u0 = [:y => 1.0])
    @test prob3.f.initialization_data !== nothing
    @test init(prob3)[x] ≈ 0.5
end

@testset "Issue#3246: type promotion with parameter dependent initialization_eqs" begin
    @variables x(t)=1 y(t)=1
    @parameters a = 1
    @named sys = ODESystem([D(x) ~ 0, D(y) ~ x + a], t; initialization_eqs = [y ~ a])

    ssys = structural_simplify(sys)
    prob = ODEProblem(ssys, [], (0, 1), [])

    @test SciMLBase.successful_retcode(solve(prob))

    seta = setsym_oop(prob, [a])
    (newu0, newp) = seta(prob, ForwardDiff.Dual{ForwardDiff.Tag{:tag, Float64}}.([1.0], 1))
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
    @mtkbuild osys = ODESystem(eqs, t)

    # Make problem.
    u0_vals = [X => 4, Y => 5.0]
    tspan = (0.0, 10.0)
    p_vals = [p => 1.0, d => 0.1]
    oprob = ODEProblem(osys, u0_vals, tspan, p_vals)
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
    affect₁ = [S ~ S - 1, I ~ I + 1]
    rate₂ = γ * I
    affect₂ = [I ~ I - 1, R ~ R + 1]
    j₁ = ConstantRateJump(rate₁, affect₁)
    j₂ = ConstantRateJump(rate₂, affect₂)
    j₃ = MassActionJump(2 * β + γ, [R => 1], [S => 1, R => -1])
    @mtkbuild js = JumpSystem([j₁, j₂, j₃], t, [S, I, R], [β, γ, S0])

    u0s = [I => 1, R => 0]
    ps = [S0 => 999, β => 0.01, γ => 0.001]
    dprob = DiscreteProblem(js, u0s, (0.0, 10.0), ps)
    @test dprob.f.initialization_data !== nothing
    sol = solve(dprob, FunctionMap())
    @test sol[S, 1] ≈ 999
    @test SciMLBase.successful_retcode(sol)

    jprob = JumpProblem(js, dprob)
    sol = solve(jprob, SSAStepper())
    @test sol[S, 1] ≈ 999
    @test SciMLBase.successful_retcode(sol)
end

@testset "Solvable array parameters with scalarized guesses" begin
    @variables x(t)
    @parameters p[1:2] q
    @mtkbuild sys = ODESystem(
        D(x) ~ p[1] + p[2] + q, t; defaults = [p[1] => q, p[2] => 2q],
        guesses = [p[1] => q, p[2] => 2q])
    @test ModelingToolkit.is_parameter_solvable(p, Dict(), defaults(sys), guesses(sys))
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0), [q => 2.0])
    @test length(ModelingToolkit.observed(prob.f.initialization_data.initializeprob.f.sys)) ==
          3
    sol = solve(prob, Tsit5())
    @test sol.ps[p] ≈ [2.0, 4.0]
end
