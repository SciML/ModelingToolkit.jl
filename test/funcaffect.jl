using ModelingToolkit, Test, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

@constants h=1 zr=0
@variables u(t)

eqs = [D(u) ~ -u]

affect1!(integ, u, p, ctx) = integ.u[u.u] += 10

@named sys = ODESystem(eqs, t, [u], [],
    discrete_events = [[4.0] => (affect1!, [u], [], [], nothing)])
prob = ODEProblem(complete(sys), [u => 10.0], (0, 10.0))
sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4 + 1][1] > 10.0

# callback
cb = ModelingToolkit.SymbolicDiscreteCallback(t == zr,
    (f = affect1!, sts = [], pars = [], discretes = [],
        ctx = [1]))
cb1 = ModelingToolkit.SymbolicDiscreteCallback(t == zr, (affect1!, [], [], [], [1]))
@test ModelingToolkit.affects(cb) isa ModelingToolkit.FunctionalAffect
@test cb == cb1
@test ModelingToolkit.SymbolicDiscreteCallback(cb) === cb # passthrough
@test hash(cb) == hash(cb1)
ModelingToolkit.generate_discrete_callback(cb, sys, ModelingToolkit.get_variables(sys),
    ModelingToolkit.get_ps(sys));

cb = ModelingToolkit.SymbolicContinuousCallback([t ~ zr],
    (f = affect1!, sts = [], pars = [], discretes = [],
        ctx = [1]))
cb1 = ModelingToolkit.SymbolicContinuousCallback([t ~ zr], (affect1!, [], [], [], [1]))
@test cb == cb1
@test ModelingToolkit.SymbolicContinuousCallback(cb) === cb # passthrough
@test hash(cb) == hash(cb1)

# named tuple
sys1 = ODESystem(eqs, t, [u], [], name = :sys,
    discrete_events = [
        [4.0] => (f = affect1!, sts = [u], pars = [], discretes = [], ctx = nothing)
    ])
@test sys == sys1

# has_functional_affect
de = ModelingToolkit.get_discrete_events(sys1)
@test length(de) == 1
de = de[1]
@test ModelingToolkit.condition(de) == [4.0]
@test ModelingToolkit.has_functional_affect(de)

sys2 = ODESystem(eqs, t, [u], [], name = :sys,
    discrete_events = [[4.0] => [u ~ -u * h]])
@test !ModelingToolkit.has_functional_affect(ModelingToolkit.get_discrete_events(sys2)[1])

# context
function affect2!(integ, u, p, ctx)
    integ.u[u.u] += ctx[1]
    ctx[1] *= 2
end
ctx1 = [10.0]
@named sys = ODESystem(eqs, t, [u], [],
    discrete_events = [[4.0, 8.0] => (affect2!, [u], [], [], ctx1)])
prob = ODEProblem(complete(sys), [u => 10.0], (0, 10.0))
sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4 + 1][1] > 10.0
i8 = findfirst(==(8.0), sol[:t])
@test sol.u[i8 + 1][1] > 20.0
@test ctx1[1] == 40.0

# parameter
function affect3!(integ, u, p, ctx)
    integ.u[u.u] += integ.ps[p.a]
    integ.ps[p.a] *= 2
end

@parameters a = 10.0
@named sys = ODESystem(eqs, t, [u], [a],
    discrete_events = [[4.0, 8.0] => (affect3!, [u], [a], [a], nothing)])
prob = ODEProblem(complete(sys), [u => 10.0], (0, 10.0))

sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4 + 1][1] > 10.0
i8 = findfirst(==(8.0), sol[:t])
@test sol.u[i8 + 1][1] > 20.0

# rename parameter
function affect3!(integ, u, p, ctx)
    integ.u[u.u] += integ.ps[p.b]
    integ.ps[p.b] *= 2
end

@named sys = ODESystem(eqs, t, [u], [a],
    discrete_events = [
        [4.0, 8.0] => (affect3!, [u], [a => :b], [a], nothing)
    ])
prob = ODEProblem(complete(sys), [u => 10.0], (0, 10.0))

sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4 + 1][1] > 10.0
i8 = findfirst(==(8.0), sol[:t])
@test sol.u[i8 + 1][1] > 20.0

# same name
@variables v(t)
@test_throws ErrorException ODESystem(eqs, t, [u], [a],
    discrete_events = [
        [4.0, 8.0] => (affect3!, [u, v => :u], [a], [a],
        nothing)
    ]; name = :sys)

@test_nowarn ODESystem(eqs, t, [u], [a],
    discrete_events = [
        [4.0, 8.0] => (affect3!, [u], [a => :u], [a], nothing)
    ]; name = :sys)

@named resistor = ODESystem(D(v) ~ v, t, [v], [])

# nested namespace
ctx = [0]
function affect4!(integ, u, p, ctx)
    ctx[1] += 1
    @test u.resistor₊v == 1
end
s1 = compose(
    ODESystem(Equation[], t, [], [], name = :s1,
        discrete_events = 1.0 => (affect4!, [resistor.v], [], [], ctx)),
    resistor)
s2 = structural_simplify(s1)
prob = ODEProblem(s2, [resistor.v => 10.0], (0, 2.01))
sol = solve(prob, Tsit5())
@test ctx[1] == 2

include("../examples/rc_model.jl")

function affect5!(integ, u, p, ctx)
    @test integ.u[u.capacitor₊v] ≈ 0.3
    integ.ps[p.C] *= 200
end

@named rc_model = ODESystem(rc_eqs, t,
    continuous_events = [
        [capacitor.v ~ 0.3] => (affect5!, [capacitor.v],
        [capacitor.C => :C], [capacitor.C], nothing)
    ])
rc_model = compose(rc_model, [resistor, capacitor, source, ground])

sys = structural_simplify(rc_model)
u0 = [capacitor.v => 0.0
      capacitor.p.i => 0.0
      resistor.v => 0.0]

prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())
@test all(sol[rc_model.capacitor.v] .< 0.4)

# hierarchical - result should be identical

function affect6!(integ, u, p, ctx)
    @test integ.u[u.v] ≈ 0.3
    integ.ps[p.C] *= 200
end

function Capacitor2(; name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C = C
    eqs = [
        D(v) ~ i / C
    ]
    extend(
        ODESystem(eqs, t, [], ps; name = name,
            continuous_events = [[v ~ 0.3] => (affect6!, [v], [C], [C], nothing)]),
        oneport)
end

@named capacitor2 = Capacitor2(C = C)

rc_eqs2 = [connect(source.p, resistor.p)
           connect(resistor.n, capacitor2.p)
           connect(capacitor2.n, source.n)
           connect(capacitor2.n, ground.g)]

@named rc_model2 = ODESystem(rc_eqs2, t)
rc_model2 = compose(rc_model2, [resistor, capacitor2, source, ground])

sys2 = structural_simplify(rc_model2)
u0 = [capacitor2.v => 0.0
      capacitor2.p.i => 0.0
      resistor.v => 0.0]

prob2 = ODEProblem(sys2, u0, (0, 10.0))
sol2 = solve(prob2, Rodas4())
@test all(sol2[rc_model2.capacitor2.v] .== sol[rc_model.capacitor.v])

# discrete events

a7_count = 0
function affect7!(integ, u, p, ctx)
    integ.ps[p.g] = 0
    ctx[1] += 1
    @test ctx[1] <= 2
    @test (ctx[1] == 1 && integ.t == 1.0) || (ctx[1] == 2 && integ.t == 2.0)
    global a7_count += 1
end

a7_ctx = [0]
function Ball(; name, g = 9.8, anti_gravity_time = 1.0)
    pars = @parameters g = g
    sts = @variables x(t), v(t)
    eqs = [D(x) ~ v, D(v) ~ g]
    ODESystem(eqs, t, sts, pars; name = name,
        discrete_events = [[anti_gravity_time] => (affect7!, [], [g], [g], a7_ctx)])
end

@named ball1 = Ball(anti_gravity_time = 1.0)
@named ball2 = Ball(anti_gravity_time = 2.0)

@named balls = ODESystem(Equation[], t)
balls = compose(balls, [ball1, ball2])

@test ModelingToolkit.has_discrete_events(balls)
@test length(ModelingToolkit.affects(ModelingToolkit.discrete_events(balls))) == 2

prob = ODEProblem(complete(balls),
    [ball1.x => 10.0, ball1.v => 0, ball2.x => 10.0, ball2.v => 0],
    (0, 3.0))
sol = solve(prob, Tsit5())

@test a7_count == 2
@test sol(0.99)[1] == sol(0.99)[3]
@test sol(1.01)[4] > sol(1.01)[2]
@test sol(1.99)[2] == sol(1.01)[2]
@test sol(1.99)[4] > sol(1.01)[4]
@test sol(2.5)[4] == sol(3.0)[4]

# bouncing ball

# DiffEq implementation
function f_(du, u, p, t)
    du[1] = u[2]
    du[2] = -p
end

function condition_(u, t, integrator) # Event when event_f(u,t) == 0
    u[1]
end

function affect_!(integrator)
    integrator.u[2] = -integrator.u[2]
end

cb_ = ContinuousCallback(condition_, affect_!)

u0 = [50.0, 0.0]
tspan = (0.0, 15.0)
p = 9.8
prob_ = ODEProblem(f_, u0, tspan, p)
sol_ = solve(prob_, Tsit5(), callback = cb_)

# same - with MTK
sts = @variables y(t), v(t)
par = @parameters g = 9.8
bb_eqs = [D(y) ~ v
          D(v) ~ -g]

function bb_affect!(integ, u, p, ctx)
    integ.u[u.v] = -integ.u[u.v]
end

@named bb_model = ODESystem(bb_eqs, t, sts, par,
    continuous_events = [
        [y ~ zr] => (bb_affect!, [v], [], [], nothing)
    ])

bb_sys = structural_simplify(bb_model)
@test only(ModelingToolkit.affects(ModelingToolkit.continuous_events(bb_sys))) isa
      ModelingToolkit.FunctionalAffect

u0 = [v => 0.0, y => 50.0]

bb_prob = ODEProblem(bb_sys, u0, (0, 15.0))
bb_sol = solve(bb_prob, Tsit5())

@test bb_sol[y] ≈ map(u -> u[1], sol_.u)
@test bb_sol[v] ≈ map(u -> u[2], sol_.u)
