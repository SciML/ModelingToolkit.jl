using ModelingToolkit, Test, OrdinaryDiffEq

@parameters t
@constants h=1 zr=0
@variables u(t)
D = Differential(t)

eqs = [D(u) ~ -u]

affect1!(integ, u, p, ctx) = integ.u[u.u] += 10

@named sys = ODESystem(eqs, t, [u], [],
                       discrete_events = [[4.0] => (affect1!, [u], [], nothing)])
prob = ODEProblem(sys, [u => 10.0], (0, 10.0))
sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4 + 1][1] > 10.0

# callback
cb = ModelingToolkit.SymbolicDiscreteCallback(t == zr,
                                              (f = affect1!, sts = [], pars = [],
                                               ctx = [1]))
cb1 = ModelingToolkit.SymbolicDiscreteCallback(t == zr, (affect1!, [], [], [1]))
@test ModelingToolkit.affects(cb) isa ModelingToolkit.FunctionalAffect
@test cb == cb1
@test ModelingToolkit.SymbolicDiscreteCallback(cb) === cb # passthrough
@test hash(cb) == hash(cb1)
ModelingToolkit.generate_discrete_callback(cb, sys, ModelingToolkit.get_variables(sys),
                                           ModelingToolkit.get_ps(sys));

cb = ModelingToolkit.SymbolicContinuousCallback([t ~ zr],
                                                (f = affect1!, sts = [], pars = [],
                                                 ctx = [1]))
cb1 = ModelingToolkit.SymbolicContinuousCallback([t ~ zr], (affect1!, [], [], [1]))
@test cb == cb1
@test ModelingToolkit.SymbolicContinuousCallback(cb) === cb # passthrough
@test hash(cb) == hash(cb1)

# named tuple
sys1 = ODESystem(eqs, t, [u], [], name = :sys,
                 discrete_events = [
                     [4.0] => (f = affect1!, sts = [u], pars = [], ctx = nothing),
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
                       discrete_events = [[4.0, 8.0] => (affect2!, [u], [], ctx1)])
prob = ODEProblem(sys, [u => 10.0], (0, 10.0))
sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4 + 1][1] > 10.0
i8 = findfirst(==(8.0), sol[:t])
@test sol.u[i8 + 1][1] > 20.0
@test ctx1[1] == 40.0

# parameter
function affect3!(integ, u, p, ctx)
    integ.u[u.u] += integ.p[p.a]
    integ.p[p.a] *= 2
end

@parameters a = 10.0
@named sys = ODESystem(eqs, t, [u], [a],
                       discrete_events = [[4.0, 8.0] => (affect3!, [u], [a], nothing)])
prob = ODEProblem(sys, [u => 10.0], (0, 10.0))

sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4 + 1][1] > 10.0
i8 = findfirst(==(8.0), sol[:t])
@test sol.u[i8 + 1][1] > 20.0

# rename parameter
function affect3!(integ, u, p, ctx)
    integ.u[u.u] += integ.p[p.b]
    integ.p[p.b] *= 2
end

@named sys = ODESystem(eqs, t, [u], [a],
                       discrete_events = [
                           [4.0, 8.0] => (affect3!, [u], [a => :b], nothing),
                       ])
prob = ODEProblem(sys, [u => 10.0], (0, 10.0))

sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4 + 1][1] > 10.0
i8 = findfirst(==(8.0), sol[:t])
@test sol.u[i8 + 1][1] > 20.0

# same name
@variables v(t)
@test_throws ErrorException ODESystem(eqs, t, [u], [a],
                                      discrete_events = [
                                          [4.0, 8.0] => (affect3!, [u, v => :u], [a],
                                                         nothing),
                                      ]; name = :sys)

@test_nowarn ODESystem(eqs, t, [u], [a],
                       discrete_events = [
                           [4.0, 8.0] => (affect3!, [u], [a => :u], nothing),
                       ]; name = :sys)

@named resistor = ODESystem(D(v) ~ v, t, [v], [])

# nested namespace
ctx = [0]
function affect4!(integ, u, p, ctx)
    ctx[1] += 1
    @test u.resistor₊v == 1
end
s1 = compose(ODESystem(Equation[], t, [], [], name = :s1,
                       discrete_events = 1.0 => (affect4!, [resistor.v], [], ctx)),
             resistor)
s2 = structural_simplify(s1)
prob = ODEProblem(s2, [resistor.v => 10.0], (0, 2.01))
sol = solve(prob, Tsit5())
@test ctx[1] == 2

include("../examples/rc_model.jl")

function affect5!(integ, u, p, ctx)
    @test integ.u[u.capacitor₊v] ≈ 0.3
    integ.p[p.C] *= 200
end

@named rc_model = ODESystem(rc_eqs, t,
                            continuous_events = [
                                [capacitor.v ~ 0.3] => (affect5!, [capacitor.v],
                                                        [capacitor.C => :C], nothing),
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
    integ.p[p.C] *= 200
end

function Capacitor2(; name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C = C
    D = Differential(t)
    eqs = [
        D(v) ~ i / C,
    ]
    extend(ODESystem(eqs, t, [], ps; name = name,
                     continuous_events = [[v ~ 0.3] => (affect6!, [v], [C], nothing)]),
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
    integ.p[p.g] = 0
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
              discrete_events = [[anti_gravity_time] => (affect7!, [], [g], a7_ctx)])
end

@named ball1 = Ball(anti_gravity_time = 1.0)
@named ball2 = Ball(anti_gravity_time = 2.0)

@named balls = ODESystem(Equation[], t)
balls = compose(balls, [ball1, ball2])

@test ModelingToolkit.has_discrete_events(balls)
@test length(ModelingToolkit.affects(ModelingToolkit.discrete_events(balls))) == 2

prob = ODEProblem(balls, [ball1.x => 10.0, ball1.v => 0, ball2.x => 10.0, ball2.v => 0],
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
                                [y ~ zr] => (bb_affect!, [v], [], nothing),
                            ])

bb_sys = structural_simplify(bb_model)
@test ModelingToolkit.affects(ModelingToolkit.continuous_events(bb_sys)) isa
      ModelingToolkit.FunctionalAffect

u0 = [v => 0.0, y => 50.0]

bb_prob = ODEProblem(bb_sys, u0, (0, 15.0))
bb_sol = solve(bb_prob, Tsit5())

@test bb_sol[y] ≈ map(u -> u[1], sol_.u)
@test bb_sol[v] ≈ map(u -> u[2], sol_.u)

@testset "SymbolicIterativeCallback" begin
    # This is loosely based on the Population example from the docs
    # `t` is the time
    # `N` is the population size
    # `m` is an additive quantity to add to the population size `N`
    # `M` is an additive quantity to add to the population size `N`
    # `α` is the limit value of `N` for large `t`
    @parameters m=-10 M=50 α=100.0
    @variables t N(t)

    Dt = Differential(t)
    eqs = [Dt(N) ~ α - N]
    #=
    NOTE
    The differential equation is 
    d/dt N(t) = α - N
    Suppose that N(0) >= α, then the solution is 
    N(t) = (N0 - α) * exp(-t) + α
    Otherwise, the solution is 
    N(t) = α - (α - N0) * exp(-t)
    =#

    u0 = [N => 0.0]
    tspan = (0.0, 20.0)

    # Until the first call to `user_affect!` hits, the population develops according to
    # N(t) = 100 - 100 * exp(-t)
    # At ``t₁ ≈ ln(2)`` the population will be ``N(t₁) ≈ 50``.
    t1 = 1.1 * log(2)
    N1 = 100 - 100 * exp(-t1)

    # From t1 onwards, the solution is 
    # N(t) = (N1+M-α) * exp(t1 - t) + α
    # N(t) = (N1+M-α) * exp(t1) * exp(-t) + α
    # At ``t₂ = t₁ + ln(2) = 2 ln(2)`` the population will roughly be
    # N(t) = 0.5 (N1+M-α) * exp(t1) * exp(-t) + α
    t2 = 2 * t1
    N2 = (N1 + 50 - 100) * exp(t1) * exp(-t2) + 100

    time_choice_history = Float64[]
    solution_history = Float64[]

    function time_choice(integ, sts, pars, ctx)
        local N = integ.u[sts.N]
        local t = integ.t

        tnext = if N <= 50
            t1
        elseif t < t2
            t2
        elseif t <= 15
            t + 1
        else
            nothing
        end
        if !isnothing(tnext)
            push!(time_choice_history, tnext)
        end
        return tnext
    end

    function user_affect!(integ, sts, pars, ctx)
        local N = integ.u[sts.N]
        local α = integ.p[pars.α]
        local m = integ.p[pars.m]
        local M = integ.p[pars.M]
        local t = integ.t

        push!(solution_history, N)

        if t == t1
            @test abs(N - N1) <= 1e-3
            integ.u[sts.N] += M # now the population is > 100 and begins to exponentially shrink
        elseif t == t2
            @test abs(N - N2) <= 1e-3
            integ.u[sts.N] += m
        else
            if N > α - m
                integ.u[sts.N] += m
            else
                integ.u[sts.N] += M
            end
        end
    end

    cb = ModelingToolkit.SymbolicIterativeCallback((time_choice, [N], [], nothing),
                                                   (user_affect!, [N], [α, m, M], nothing))

    @test cb.initial_affect == false

    @named osys = ODESystem(eqs, t, [N], [α, m, M]; discrete_events = [cb])

    @test ModelingToolkit.has_discrete_events(osys)
    @test length(ModelingToolkit.get_discrete_events(osys)) == 1
    @test only(ModelingToolkit.get_discrete_events(osys)).wrapped == cb

    oprob = ODEProblem(osys, u0, tspan)
    sol = solve(oprob, Tsit5())

    @test all(abs(first(sol(t)) - _N) <= 1e-3
              for (t, _N) in zip(time_choice_history, solution_history))

    # Now: `initial_affect = true`
    function time_choice2(integ, sts, pars, ctx)
        return nothing
    end

    function user_affect!2(integ, sts, pars, ctx)
        integ.u[sts.N] += integ.p[pars.M]
    end

    cb2 = ModelingToolkit.SymbolicIterativeCallback((time_choice2, [], [], nothing),
                                                    (user_affect!2, [N], [M], nothing),
                                                    true)
    @test cb2.initial_affect == true

    @named osys2 = ODESystem(eqs, t, [N], [α, m, M]; discrete_events = [cb2])

    @test ModelingToolkit.has_discrete_events(osys2)
    @test length(ModelingToolkit.get_discrete_events(osys2)) == 1
    @test only(ModelingToolkit.get_discrete_events(osys2)).wrapped == cb2

    oprob2 = ODEProblem(osys2, u0, tspan)
    sol2 = solve(oprob2, Tsit5())

    @test first(sol2(1e-3)) >= 50
end

@testset "SymbolicPeriodicCallback" begin
    @parameters α=100 m=20
    @variables t N(t)
    Dt = Differential(t)
    eqs = [Dt(N) ~ α - N]

    del_t = log(2) / 2

    u0 = [N => 0.0]
    tspan = (0.0, 20.0)

    T = Float64[]

    function user_affect!(integ, sts, pars, ctx)
        push!(T, integ.t)
        N = integ.u[sts.N]
        if N > 70
            integ.u[sts.N] -= integ.p[pars.m]
        end
    end

    # define discrete symbolic callback
    cb = ModelingToolkit.SymbolicPeriodicCallback(del_t, (user_affect!, [N], [m], nothing))

    @test cb.initial_affect == false
    @test cb.final_affect == false

    # setup system, problem, and solve:
    @named osys = ODESystem(eqs, t, [N], [α, m]; discrete_events = [cb])

    @test ModelingToolkit.has_discrete_events(osys)
    @test length(ModelingToolkit.get_discrete_events(osys)) == 1
    @test only(ModelingToolkit.get_discrete_events(osys)).wrapped == cb

    oprob = ODEProblem(osys, u0, tspan)
    sol = solve(oprob, Tsit5())

    for _t in T
        if first(sol(_t)) > 70
            @test first(sol(_t + 1e-5)) < 70
        end
    end

    # TODO initial_affect=true, final_affect=true
end

@testset "SymbolicPresetTimeCallback" begin
    @parameters α=100 m=20
    @variables t N(t)
    Dt = Differential(t)
    eqs = [Dt(N) ~ α - N]

    u0 = [N => 0.0]
    tspan = (0.0, 20.0)
    tstops = LinRange(tspan[1], tspan[2], 100)

    T = Float64[]

    function user_affect!(integ, sts, pars, ctx)
        push!(T, integ.t)
        N = integ.u[sts.N]
        if N > 70
            integ.u[sts.N] -= integ.p[pars.m]
        end
    end

    # define discrete symbolic callback
    cb = ModelingToolkit.SymbolicPresetTimeCallback(tstops,
                                                    (user_affect!, [N], [m], nothing))

    @test cb.filter_tstops == true

    # setup system, problem, and solve:
    @named osys = ODESystem(eqs, t, [N], [α, m]; discrete_events = [cb])

    @test ModelingToolkit.has_discrete_events(osys)
    @test length(ModelingToolkit.get_discrete_events(osys)) == 1
    @test only(ModelingToolkit.get_discrete_events(osys)).wrapped == cb

    oprob = ODEProblem(osys, u0, tspan)
    sol = solve(oprob, Tsit5())

    for _t in T
        if first(sol(_t)) > 70
            @test first(sol(_t + 1e-5)) < 70
        end
    end

    # TODO filter_tstops=false
end
