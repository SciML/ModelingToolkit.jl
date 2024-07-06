using ModelingToolkit, OrdinaryDiffEq, NonlinearSolve, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

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
@test initprob isa NonlinearProblem
sol = solve(initprob)
@test SciMLBase.successful_retcode(sol)
@test sol.u == [0.0, 0.0, 0.0, 0.0]
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
        p(t) = p′
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
    @test_throws ArgumentError ODEProblem(ssys, [x => 3], (0, 1), []) # y should have a guess
end

# https://github.com/SciML/ModelingToolkit.jl/issues/2619

@parameters k1 k2 ω
@variables X(t) Y(t)
eqs_1st_order = [D(Y) + Y - ω ~ 0,
       X + k1 ~ Y + k2]
eqs_2nd_order = [D(D(Y)) + 2ω*D(Y) + (ω^2)*Y ~ 0,
        X + k1 ~ Y + k2]
@mtkbuild sys_1st_order = ODESystem(eqs_1st_order, t)
@mtkbuild sys_2nd_order = ODESystem(eqs_2nd_order, t)

u0_1st_order_1 = [X => 1.0, Y => 2.0]
u0_1st_order_2 = [Y => 2.0]
u0_2nd_order_1 = [X => 1.0, Y => 2.0, D(Y) => 0.5]
u0_2nd_order_2 = [Y => 2.0, D(Y) => 0.5]
tspan = (0.0, 10.)
ps = [ω => 0.5, k1 => 2.0, k2 => 3.0]

oprob_1st_order_1 = ODEProblem(sys_1st_order, u0_1st_order_1, tspan, ps)
oprob_1st_order_2 = ODEProblem(sys_1st_order, u0_1st_order_2, tspan, ps)
oprob_2nd_order_1 = ODEProblem(sys_2nd_order, u0_2nd_order_1, tspan, ps) # gives sys_2nd_order
oprob_2nd_order_2 = ODEProblem(sys_2nd_order, u0_2nd_order_2, tspan, ps)

@test solve(oprob_1st_order_1, Rosenbrock23()).retcode == SciMLBase.ReturnCode.InitialFailure
@test solve(oprob_1st_order_2, Rosenbrock23())[Y][1] == 2.0
@test solve(oprob_2nd_order_1, Rosenbrock23()).retcode == SciMLBase.ReturnCode.InitialFailure
sol = solve(oprob_2nd_order_2, Rosenbrock23()) # retcode: Success
@test sol[Y][1] == 2.0
@test sol[1][2] == 0.5
