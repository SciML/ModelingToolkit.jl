using ModelingToolkit, OrdinaryDiffEq, StochasticDiffEq, JumpProcesses, Test
using SciMLStructures: canonicalize, Discrete
using ModelingToolkit: SymbolicContinuousCallback,
                       SymbolicDiscreteCallback,
                       t_nounits as t,
                       D_nounits as D,
                       affects, affect_negs, system, observed, AffectSystem
using StableRNGs
import SciMLBase
using SymbolicIndexingInterface
using Setfield
rng = StableRNG(12345)

function get_callback(prob)
    prob.kwargs[:callback]
end

@variables x(t) = 0

eqs = [D(x) ~ 1]
affect = [x ~ 0]
affect_neg = [x ~ 1]

@testset "SymbolicContinuousCallback constructors" begin
    e = SymbolicContinuousCallback(eqs[])
    @test e isa SymbolicContinuousCallback
    @test isequal(equations(e), eqs)
    @test e.affect === nothing
    @test e.affect_neg === nothing
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs)
    @test e isa SymbolicContinuousCallback
    @test isequal(equations(e), eqs)
    @test e.affect === nothing
    @test e.affect_neg === nothing
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs, nothing)
    @test e isa SymbolicContinuousCallback
    @test isequal(equations(e), eqs)
    @test e.affect === nothing
    @test e.affect_neg === nothing
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs[], nothing)
    @test e isa SymbolicContinuousCallback
    @test isequal(equations(e), eqs)
    @test e.affect === nothing
    @test e.affect_neg === nothing
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs => nothing)
    @test e isa SymbolicContinuousCallback
    @test isequal(equations(e), eqs)
    @test e.affect === nothing
    @test e.affect_neg === nothing
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs[] => nothing)
    @test e isa SymbolicContinuousCallback
    @test isequal(equations(e), eqs)
    @test e.affect === nothing
    @test e.affect_neg === nothing
    @test e.rootfind == SciMLBase.LeftRootFind

    ## With affect
    e = SymbolicContinuousCallback(eqs[], affect)
    @test e isa SymbolicContinuousCallback
    @test isequal(equations(e), eqs)
    @test e.rootfind == SciMLBase.LeftRootFind

    # with only positive edge affect
    e = SymbolicContinuousCallback(eqs[], affect, affect_neg = nothing)
    @test e isa SymbolicContinuousCallback
    @test isequal(equations(e), eqs)
    @test isnothing(e.affect_neg)
    @test e.rootfind == SciMLBase.LeftRootFind

    # with explicit edge affects
    e = SymbolicContinuousCallback(eqs[], affect, affect_neg = affect_neg)
    @test e isa SymbolicContinuousCallback
    @test isequal(equations(e), eqs)
    @test e.rootfind == SciMLBase.LeftRootFind

    # with different root finding ops
    e = SymbolicContinuousCallback(
        eqs[], affect, affect_neg = affect_neg, rootfind = SciMLBase.LeftRootFind)
    @test e isa SymbolicContinuousCallback
    @test isequal(equations(e), eqs)
    @test e.rootfind == SciMLBase.LeftRootFind
end

@testset "ImperativeAffect constructors" begin
    fmfa(o, x, i, c) = nothing
    m = ModelingToolkit.ImperativeAffect(fmfa)
    @test m isa ModelingToolkit.ImperativeAffect
    @test m.f == fmfa
    @test m.obs == []
    @test m.obs_syms == []
    @test m.modified == []
    @test m.mod_syms == []
    @test m.ctx === nothing

    m = ModelingToolkit.ImperativeAffect(fmfa, (;))
    @test m isa ModelingToolkit.ImperativeAffect
    @test m.f == fmfa
    @test m.obs == []
    @test m.obs_syms == []
    @test m.modified == []
    @test m.mod_syms == []
    @test m.ctx === nothing

    m = ModelingToolkit.ImperativeAffect(fmfa, (; x))
    @test m isa ModelingToolkit.ImperativeAffect
    @test m.f == fmfa
    @test isequal(m.obs, [])
    @test m.obs_syms == []
    @test isequal(m.modified, [x])
    @test m.mod_syms == [:x]
    @test m.ctx === nothing

    m = ModelingToolkit.ImperativeAffect(fmfa, (; y = x))
    @test m isa ModelingToolkit.ImperativeAffect
    @test m.f == fmfa
    @test isequal(m.obs, [])
    @test m.obs_syms == []
    @test isequal(m.modified, [x])
    @test m.mod_syms == [:y]
    @test m.ctx === nothing

    m = ModelingToolkit.ImperativeAffect(fmfa; observed = (; y = x))
    @test m isa ModelingToolkit.ImperativeAffect
    @test m.f == fmfa
    @test isequal(m.obs, [x])
    @test m.obs_syms == [:y]
    @test m.modified == []
    @test m.mod_syms == []
    @test m.ctx === nothing

    m = ModelingToolkit.ImperativeAffect(fmfa; modified = (; x))
    @test m isa ModelingToolkit.ImperativeAffect
    @test m.f == fmfa
    @test isequal(m.obs, [])
    @test m.obs_syms == []
    @test isequal(m.modified, [x])
    @test m.mod_syms == [:x]
    @test m.ctx === nothing

    m = ModelingToolkit.ImperativeAffect(fmfa; modified = (; y = x))
    @test m isa ModelingToolkit.ImperativeAffect
    @test m.f == fmfa
    @test isequal(m.obs, [])
    @test m.obs_syms == []
    @test isequal(m.modified, [x])
    @test m.mod_syms == [:y]
    @test m.ctx === nothing

    m = ModelingToolkit.ImperativeAffect(fmfa, (; x), (; x))
    @test m isa ModelingToolkit.ImperativeAffect
    @test m.f == fmfa
    @test isequal(m.obs, [x])
    @test m.obs_syms == [:x]
    @test isequal(m.modified, [x])
    @test m.mod_syms == [:x]
    @test m.ctx === nothing

    m = ModelingToolkit.ImperativeAffect(fmfa, (; y = x), (; y = x))
    @test m isa ModelingToolkit.ImperativeAffect
    @test m.f == fmfa
    @test isequal(m.obs, [x])
    @test m.obs_syms == [:y]
    @test isequal(m.modified, [x])
    @test m.mod_syms == [:y]
    @test m.ctx === nothing

    m = ModelingToolkit.ImperativeAffect(
        fmfa; modified = (; y = x), observed = (; y = x))
    @test m isa ModelingToolkit.ImperativeAffect
    @test m.f == fmfa
    @test isequal(m.obs, [x])
    @test m.obs_syms == [:y]
    @test isequal(m.modified, [x])
    @test m.mod_syms == [:y]
    @test m.ctx === nothing

    m = ModelingToolkit.ImperativeAffect(
        fmfa; modified = (; y = x), observed = (; y = x), ctx = 3)
    @test m isa ModelingToolkit.ImperativeAffect
    @test m.f == fmfa
    @test isequal(m.obs, [x])
    @test m.obs_syms == [:y]
    @test isequal(m.modified, [x])
    @test m.mod_syms == [:y]
    @test m.ctx === 3

    m = ModelingToolkit.ImperativeAffect(fmfa, (; x), (; x), 3)
    @test m isa ModelingToolkit.ImperativeAffect
    @test m.f == fmfa
    @test isequal(m.obs, [x])
    @test m.obs_syms == [:x]
    @test isequal(m.modified, [x])
    @test m.mod_syms == [:x]
    @test m.ctx === 3
end

@testset "Condition Compilation" begin
    @named sys = System(eqs, t, continuous_events = [x ~ 1])
    cevt1 = getfield(sys, :continuous_events)[]
    @test getfield(sys, :continuous_events)[] ==
          SymbolicContinuousCallback(Equation[x ~ 1], nothing; zero_crossing_id = cevt1.zero_crossing_id)
    @test isequal(equations(getfield(sys, :continuous_events))[], x ~ 1)
    fsys = flatten(sys)
    @test isequal(equations(getfield(fsys, :continuous_events))[], x ~ 1)

    @named sys2 = System([D(x) ~ 1], t, continuous_events = [x ~ 2], systems = [sys])
    cevt2 = getfield(sys2, :continuous_events)[]
    @test getfield(sys2, :continuous_events)[] ==
          SymbolicContinuousCallback(Equation[x ~ 2], nothing; zero_crossing_id = cevt2.zero_crossing_id)
    @test all(ModelingToolkit.continuous_events(sys2) .== [
        SymbolicContinuousCallback(Equation[x ~ 2], nothing; zero_crossing_id = cevt2.zero_crossing_id),
        SymbolicContinuousCallback(Equation[sys.x ~ 1], nothing; zero_crossing_id = cevt1.zero_crossing_id)
    ])

    @test isequal(equations(getfield(sys2, :continuous_events))[1], x ~ 2)
    @test length(ModelingToolkit.continuous_events(sys2)) == 2
    @test isequal(equations(ModelingToolkit.continuous_events(sys2)[1])[], x ~ 2)
    @test isequal(equations(ModelingToolkit.continuous_events(sys2)[2])[], sys.x ~ 1)

    sys = complete(sys)
    sys_nosplit = complete(sys; split = false)
    sys2 = complete(sys2)

    # Test proper rootfinding
    prob = ODEProblem(sys, Pair[], (0.0, 2.0))
    p0 = 0
    t0 = 0
    @test get_callback(prob) isa ModelingToolkit.DiffEqCallbacks.ContinuousCallback
    cb = ModelingToolkit.generate_continuous_callbacks(sys)
    cond = cb.condition
    out = [0.0]
    cond.f(out, [0], p0, t0)
    @test out[] ≈ -1 # signature is u,p,t
    cond.f(out, [1], p0, t0)
    @test out[] ≈ 0  # signature is u,p,t
    cond.f(out, [2], p0, t0)
    @test out[] ≈ 1  # signature is u,p,t

    prob = ODEProblem(sys, Pair[], (0.0, 2.0))
    prob_nosplit = ODEProblem(sys_nosplit, Pair[], (0.0, 2.0))
    sol = solve(prob, Tsit5())
    sol_nosplit = solve(prob_nosplit, Tsit5())
    @test minimum(t -> abs(t - 1), sol.t) < 1e-10 # test that the solver stepped at the root
    @test minimum(t -> abs(t - 1), sol_nosplit.t) < 1e-10 # test that the solver stepped at the root

    # Test user-provided callback is respected
    test_callback = DiscreteCallback(x -> x, x -> x)
    prob = ODEProblem(sys, Pair[], (0.0, 2.0), callback = test_callback)
    prob_nosplit = ODEProblem(sys_nosplit, Pair[], (0.0, 2.0), callback = test_callback)
    cbs = get_callback(prob)
    cbs_nosplit = get_callback(prob_nosplit)
    @test cbs isa CallbackSet
    @test cbs.discrete_callbacks[1] == test_callback
    @test cbs_nosplit isa CallbackSet
    @test cbs_nosplit.discrete_callbacks[1] == test_callback

    prob = ODEProblem(sys2, Pair[], (0.0, 3.0))
    cb = get_callback(prob)
    @test cb isa ModelingToolkit.DiffEqCallbacks.VectorContinuousCallback

    cond = cb.condition
    out = [0.0, 0.0]
    # the root to find is 2
    cond.f(out, [0, 0], p0, t0)
    @test out[1] ≈ -2 # signature is u,p,t
    cond.f(out, [1, 0], p0, t0)
    @test out[1] ≈ -1  # signature is u,p,t
    cond.f(out, [2, 0], p0, t0) # this should return 0
    @test out[1] ≈ 0  # signature is u,p,t

    # the root to find is 1
    out = [0.0, 0.0]
    cond.f(out, [0, 0], p0, t0)
    @test out[2] ≈ -1 # signature is u,p,t
    cond.f(out, [0, 1], p0, t0) # this should return 0
    @test out[2] ≈ 0  # signature is u,p,t
    cond.f(out, [0, 2], p0, t0)
    @test out[2] ≈ 1  # signature is u,p,t

    sol = solve(prob, Tsit5())
    @test minimum(t -> abs(t - 1), sol.t) < 1e-9 # test that the solver stepped at the first root
    @test minimum(t -> abs(t - 2), sol.t) < 1e-9 # test that the solver stepped at the second root

    @named sys = System(eqs, t, continuous_events = [x ~ 1, x ~ 2]) # two root eqs using the same unknown
    sys = complete(sys)
    prob = ODEProblem(sys, Pair[], (0.0, 3.0))
    @test get_callback(prob) isa ModelingToolkit.DiffEqCallbacks.VectorContinuousCallback
    sol = solve(prob, Tsit5())
    @test minimum(t -> abs(t - 1), sol.t) < 1e-9 # test that the solver stepped at the first root
    @test minimum(t -> abs(t - 2), sol.t) < 1e-9 # test that the solver stepped at the second root
end

@testset "Bouncing Ball" begin
    ###### 1D Bounce
    @variables x(t)=1 v(t)=0

    root_eqs = [x ~ 0]
    affect = [v ~ -Pre(v)]

    @named ball = System(
        [D(x) ~ v
         D(v) ~ -9.8], t, continuous_events = root_eqs => affect)

    cev = only(continuous_events(ball))
    @test isequal(only(equations(cev)), x ~ 0)
    @test isequal(only(cev.affect.affect), v ~ -Pre(v))
    ball = mtkcompile(ball)

    @test length(ModelingToolkit.continuous_events(ball)) == 1
    cev = only(continuous_events(ball))
    @test isequal(only(equations(cev)), x ~ 0)
    @test isequal(only(observed(cev.affect.system)), v ~ -Pre(v))

    tspan = (0.0, 5.0)
    prob = ODEProblem(ball, Pair[], tspan)
    sol = solve(prob, Tsit5())
    @test 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close

    ###### 2D bouncing ball
    @variables x(t)=1 y(t)=0 vx(t)=0 vy(t)=1

    events = [[x ~ 0] => [vx ~ -Pre(vx)]
              [y ~ -1.5, y ~ 1.5] => [vy ~ -Pre(vy)]]

    @named ball = System(
        [D(x) ~ vx
         D(y) ~ vy
         D(vx) ~ -9.8
         D(vy) ~ -0.01vy], t; continuous_events = events)

    _ball = ball
    ball = mtkcompile(_ball)
    ball_nosplit = mtkcompile(_ball; split = false)

    tspan = (0.0, 5.0)
    prob = ODEProblem(ball, Pair[], tspan)
    prob_nosplit = ODEProblem(ball_nosplit, Pair[], tspan)

    cb = get_callback(prob)
    @test cb isa ModelingToolkit.DiffEqCallbacks.VectorContinuousCallback
    _cevs = getfield(ball, :continuous_events)
    @test isequal(only(equations(_cevs[1])), x ~ 0)
    @test isequal(only(observed(_cevs[1].affect.system)), vx ~ -Pre(vx))
    @test issetequal(equations(_cevs[2]), [y ~ -1.5, y ~ 1.5])
    @test isequal(only(observed(_cevs[2].affect.system)), vy ~ -Pre(vy))
    cond = cb.condition
    out = [0.0, 0.0, 0.0]
    p0 = 0.0
    t0 = 0.0
    cond.f(out, [0, 0, 0, 0], p0, t0)
    @test out ≈ [0, 1.5, -1.5]

    sol = solve(prob, Tsit5())
    sol_nosplit = solve(prob_nosplit, Tsit5())
    @test 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close
    @test minimum(sol[y]) ≈ -1.5 # check wall conditions
    @test maximum(sol[y]) ≈ 1.5  # check wall conditions
    @test 0 <= minimum(sol_nosplit[x]) <= 1e-10 # the ball never went through the floor but got very close
    @test minimum(sol_nosplit[y]) ≈ -1.5 # check wall conditions
    @test maximum(sol_nosplit[y]) ≈ 1.5  # check wall conditions

    ## Test multi-variable affect
    # in this test, there are two variables affected by a single event.
    events = [[x ~ 0] => [vx ~ -Pre(vx), vy ~ -Pre(vy)]]

    @named ball = System(
        [D(x) ~ vx
         D(y) ~ vy
         D(vx) ~ -1
         D(vy) ~ 0], t; continuous_events = events)

    ball_nosplit = mtkcompile(ball)
    ball = mtkcompile(ball)

    tspan = (0.0, 5.0)
    prob = ODEProblem(ball, Pair[], tspan)
    prob_nosplit = ODEProblem(ball_nosplit, Pair[], tspan)
    sol = solve(prob, Tsit5())
    sol_nosplit = solve(prob_nosplit, Tsit5())
    @test 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close
    @test -minimum(sol[y]) ≈ maximum(sol[y]) ≈ sqrt(2)  # the ball will never go further than √2 in either direction (gravity was changed to 1 to get this particular number)
    @test 0 <= minimum(sol_nosplit[x]) <= 1e-10 # the ball never went through the floor but got very close
    @test -minimum(sol_nosplit[y]) ≈ maximum(sol_nosplit[y]) ≈ sqrt(2)  # the ball will never go further than √2 in either direction (gravity was changed to 1 to get this particular number)
end

# issue https://github.com/SciML/ModelingToolkit.jl/issues/1386
# tests that it works for ODAESystem
@testset "ODAESystem" begin
    @variables vs(t) v(t) vmeasured(t)
    eq = [vs ~ sin(2pi * t)
          D(v) ~ vs - v
          D(vmeasured) ~ 0.0]
    ev = [sin(20pi * t) ~ 0.0] => [vmeasured ~ Pre(v)]
    @named sys = System(eq, t, continuous_events = ev)
    sys = mtkcompile(sys)
    prob = ODEProblem(sys, zeros(2), (0.0, 5.1))
    sol = solve(prob, Tsit5())
    @test all(minimum((0:0.1:5) .- sol.t', dims = 2) .< 0.0001) # test that the solver stepped every 0.1s as dictated by event
    @test sol([0.25 - eps()])[vmeasured][] == sol([0.23])[vmeasured][] # test the hold property
end

##  https://github.com/SciML/ModelingToolkit.jl/issues/1528
@testset "Handle Empty Events" begin
    Dₜ = D

    @parameters u(t) [input = true]  # Indicate that this is a controlled input
    @parameters y(t) [output = true] # Indicate that this is a measured output

    function Mass(; name, m = 1.0, p = 0, v = 0)
        ps = @parameters m = m
        sts = @variables pos(t)=p vel(t)=v
        eqs = Dₜ(pos) ~ vel
        System(eqs, t, [pos, vel], ps; name)
    end
    function Spring(; name, k = 1e4)
        ps = @parameters k = k
        @variables x(t) = 0 # Spring deflection
        System(Equation[], t, [x], ps; name)
    end
    function Damper(; name, c = 10)
        ps = @parameters c = c
        @variables vel(t) = 0
        System(Equation[], t, [vel], ps; name)
    end
    function SpringDamper(; name, k = false, c = false)
        spring = Spring(; name = :spring, k)
        damper = Damper(; name = :damper, c)
        compose(System(Equation[], t; name),
            spring, damper)
    end
    connect_sd(
        sd, m1, m2) = [
        sd.spring.x ~ m1.pos - m2.pos, sd.damper.vel ~ m1.vel - m2.vel]
    sd_force(sd) = -sd.spring.k * sd.spring.x - sd.damper.c * sd.damper.vel
    @named mass1 = Mass(; m = 1)
    @named mass2 = Mass(; m = 1)
    @named sd = SpringDamper(; k = 1000, c = 10)
    function Model(u, d = 0)
        eqs = [connect_sd(sd, mass1, mass2)
               Dₜ(mass1.vel) ~ (sd_force(sd) + u) / mass1.m
               Dₜ(mass2.vel) ~ (-sd_force(sd) + d) / mass2.m]
        @named _model = System(eqs, t; observed = [y ~ mass2.pos])
        @named model = compose(_model, mass1, mass2, sd)
    end
    model = Model(sin(30t))
    sys = mtkcompile(model)
    @test isempty(ModelingToolkit.continuous_events(sys))
end

@testset "SDE/ODESystem Discrete Callbacks" begin
    function testsol(
            sys, probtype, solver, u0, p, tspan; tstops = Float64[], paramtotest = nothing,
            kwargs...)
        prob = probtype(complete(sys), [u0; p], tspan; kwargs...)
        sol = solve(prob, solver(); tstops = tstops, abstol = 1e-10, reltol = 1e-10)
        @test isapprox(sol(1.0000000001)[1] - sol(0.999999999)[1], 1.0; rtol = 1e-6)
        paramtotest === nothing || (@test sol.ps[paramtotest] == [0.0, 1.0])
        @test isapprox(sol(4.0)[1], 2 * exp(-2.0); rtol = 1e-6)
        sol
    end

    @parameters k(t) t1 t2
    @variables A(t) B(t)

    cond1 = (t == t1)
    affect1 = [A ~ Pre(A) + 1]
    cb1 = cond1 => affect1
    cond2 = (t == t2)
    affect2 = [k ~ 1.0]
    cb2 = cond2 => affect2
    cb2 = SymbolicDiscreteCallback(cb2, discrete_parameters = [k], iv = t)

    ∂ₜ = D
    eqs = [∂ₜ(A) ~ -k * A]
    @named osys = System(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2])
    @named ssys = SDESystem(eqs, [0.0], t, [A], [k, t1, t2], discrete_events = [cb1, cb2])
    u0 = [A => 1.0]
    p = [k => 0.0, t1 => 1.0, t2 => 2.0]
    tspan = (0.0, 4.0)
    testsol(osys, ODEProblem, Tsit5, u0, p, tspan; tstops = [1.0, 2.0], paramtotest = k)
    testsol(ssys, SDEProblem, RI5, u0, p, tspan; tstops = [1.0, 2.0], paramtotest = k)

    cond1a = (t == t1)
    affect1a = [A ~ Pre(A) + 1, B ~ A]
    cb1a = cond1a => affect1a
    @named osys1 = System(eqs, t, [A, B], [k, t1, t2], discrete_events = [cb1a, cb2])
    @named ssys1 = SDESystem(eqs, [0.0], t, [A, B], [k, t1, t2],
        discrete_events = [cb1a, cb2])
    u0′ = [A => 1.0, B => 0.0]
    sol = testsol(osys1, ODEProblem, Tsit5, u0′, p, tspan;
        tstops = [1.0, 2.0], check_length = false, paramtotest = k)
    @test sol(1.0000001, idxs = B) == 2.0

    sol = testsol(ssys1, SDEProblem, RI5, u0′, p, tspan; tstops = [1.0, 2.0],
        check_length = false, paramtotest = k)
    @test sol(1.0000001, idxs = B) == 2.0

    # same as above - but with set-time event syntax
    cb1‵ = [1.0] => affect1 # needs to be a Vector for the event to happen only once
    cb2‵ = SymbolicDiscreteCallback([2.0] => affect2, discrete_parameters = [k], iv = t)
    @named osys‵ = System(eqs, t, [A], [k], discrete_events = [cb1‵, cb2‵])
    @named ssys‵ = SDESystem(eqs, [0.0], t, [A], [k], discrete_events = [cb1‵, cb2‵])
    testsol(osys‵, ODEProblem, Tsit5, u0, p, tspan; paramtotest = k)
    testsol(ssys‵, SDEProblem, RI5, u0, p, tspan; paramtotest = k)

    # mixing discrete affects
    @named osys3 = System(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵])
    @named ssys3 = SDESystem(eqs, [0.0], t, [A], [k, t1, t2],
        discrete_events = [cb1, cb2‵])
    testsol(osys3, ODEProblem, Tsit5, u0, p, tspan; tstops = [1.0], paramtotest = k)
    testsol(ssys3, SDEProblem, RI5, u0, p, tspan; tstops = [1.0], paramtotest = k)

    # mixing with a func affect
    function affect!(mod, obs, ctx, integ)
        return (; k = 1.0)
    end
    cb2‵‵ = [2.0] => (f = affect!, modified = (; k))
    @named osys4 = System(eqs, t, [A], [k, t1], discrete_events = [cb1, cb2‵‵])
    @named ssys4 = SDESystem(eqs, [0.0], t, [A], [k, t1],
        discrete_events = [cb1, cb2‵‵])
    oprob4 = ODEProblem(complete(osys4), [u0; p], tspan)
    testsol(osys4, ODEProblem, Tsit5, u0, p, tspan; tstops = [1.0], paramtotest = k)
    testsol(ssys4, SDEProblem, RI5, u0, p, tspan; tstops = [1.0], paramtotest = k)

    # mixing with symbolic condition in the func affect
    cb2‵‵‵ = (t == t2) => (f = affect!, modified = (; k))
    @named osys5 = System(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵‵‵])
    @named ssys5 = SDESystem(eqs, [0.0], t, [A], [k, t1, t2],
        discrete_events = [cb1, cb2‵‵‵])
    testsol(osys5, ODEProblem, Tsit5, u0, p, tspan; tstops = [1.0, 2.0])
    testsol(ssys5, SDEProblem, RI5, u0, p, tspan; tstops = [1.0, 2.0])
    @named osys6 = System(eqs, t, [A], [k, t1, t2], discrete_events = [cb2‵‵‵, cb1])
    @named ssys6 = SDESystem(eqs, [0.0], t, [A], [k, t1, t2],
        discrete_events = [cb2‵‵‵, cb1])
    testsol(osys6, ODEProblem, Tsit5, u0, p, tspan; tstops = [1.0, 2.0], paramtotest = k)
    testsol(ssys6, SDEProblem, RI5, u0, p, tspan; tstops = [1.0, 2.0], paramtotest = k)

    # mix a continuous event too
    cond3 = A ~ 0.1
    affect3 = [k ~ 0.0]
    cb3 = SymbolicContinuousCallback(cond3 => affect3, discrete_parameters = [k], iv = t)
    @named osys7 = System(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵‵‵],
        continuous_events = [cb3])
    @named ssys7 = SDESystem(eqs, [0.0], t, [A], [k, t1, t2],
        discrete_events = [cb1, cb2‵‵‵],
        continuous_events = [cb3])

    sol = testsol(osys7, ODEProblem, Tsit5, u0, p, (0.0, 10.0); tstops = [1.0, 2.0])
    @test isapprox(sol(10.0)[1], 0.1; atol = 1e-10, rtol = 1e-10)
    sol = testsol(ssys7, SDEProblem, RI5, u0, p, (0.0, 10.0); tstops = [1.0, 2.0])
    @test isapprox(sol(10.0)[1], 0.1; atol = 1e-10, rtol = 1e-10)
end

@testset "JumpSystem Discrete Callbacks" begin
    function testsol(jsys, u0, p, tspan; tstops = Float64[], paramtotest = nothing,
            N = 40000, kwargs...)
        jsys = complete(jsys)
        jprob = JumpProblem(jsys, [u0; p], tspan; aggregator = Direct(), kwargs...)
        sol = solve(jprob, SSAStepper(); tstops = tstops)
        @test (sol(1.000000000001)[1] - sol(0.99999999999)[1]) == 1
        paramtotest === nothing || (@test sol.ps[paramtotest] == [0.0, 1.0])
        @test sol(40.0)[1] == 0
        sol
    end

    @parameters k(t) t1 t2
    @variables A(t) B(t)

    eqs = [MassActionJump(k, [A => 1], [A => -1])]
    cond1 = (t == t1)
    affect1 = [A ~ Pre(A) + 1]
    cb1 = cond1 => affect1
    cond2 = (t == t2)
    affect2 = [k ~ 1.0]
    cb2 = cond2 => affect2
    cb2 = SymbolicDiscreteCallback(cb2, discrete_parameters = [k], iv = t)

    @named jsys = JumpSystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2])
    u0 = [A => 1]
    p = [k => 0.0, t1 => 1.0, t2 => 2.0]
    tspan = (0.0, 40.0)
    testsol(jsys, u0, p, tspan; tstops = [1.0, 2.0], rng, paramtotest = k)

    cond1a = (t == t1)
    affect1a = [A ~ Pre(A) + 1, B ~ A]
    cb1a = cond1a => affect1a
    @named jsys1 = JumpSystem(eqs, t, [A, B], [k, t1, t2], discrete_events = [cb1a, cb2])
    u0′ = [A => 1, B => 0]
    sol = testsol(jsys1, u0′, p, tspan; tstops = [1.0, 2.0],
        check_length = false, rng, paramtotest = k)
    @test sol(1.000000001, idxs = B) == 2

    # same as above - but with set-time event syntax
    cb1‵ = [1.0] => affect1 # needs to be a Vector for the event to happen only once
    cb2‵ = SymbolicDiscreteCallback([2.0] => affect2, discrete_parameters = [k], iv = t)
    @named jsys‵ = JumpSystem(eqs, t, [A], [k], discrete_events = [cb1‵, cb2‵])
    testsol(jsys‵, u0, [p[1]], tspan; rng, paramtotest = k)

    # mixing discrete affects
    @named jsys3 = JumpSystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵])
    testsol(jsys3, u0, p, tspan; tstops = [1.0], rng, paramtotest = k)

    # mixing with a func affect
    function affect!(mod, obs, ctx, integrator)
        return (; k = 1.0)
    end
    cb2‵‵ = [2.0] => (f = affect!, modified = (; k))
    @named jsys4 = JumpSystem(eqs, t, [A], [k, t1], discrete_events = [cb1, cb2‵‵])
    testsol(jsys4, u0, p, tspan; tstops = [1.0], rng, paramtotest = k)

    # mixing with symbolic condition in the func affect
    cb2‵‵‵ = (t == t2) => (f = affect!, modified = (; k))
    @named jsys5 = JumpSystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵‵‵])
    testsol(jsys5, u0, p, tspan; tstops = [1.0, 2.0], rng, paramtotest = k)
    @named jsys6 = JumpSystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb2‵‵‵, cb1])
    testsol(jsys6, u0, p, tspan; tstops = [1.0, 2.0], rng, paramtotest = k)
end

@testset "Namespacing" begin
    function oscillator_ce(k = 1.0; name)
        sts = @variables x(t)=1.0 v(t)=0.0 F(t)
        ps = @parameters k=k Θ=0.5
        eqs = [D(x) ~ v, D(v) ~ -k * x + F]
        ev = [x ~ Θ] => [x ~ 1.0, v ~ 0.0]
        System(eqs, t, sts, ps, continuous_events = [ev]; name)
    end

    @named oscce = oscillator_ce()
    eqs = [oscce.F ~ 0]
    @named eqs_sys = System(eqs, t)
    @named oneosc_ce = compose(eqs_sys, oscce)
    oneosc_ce_simpl = mtkcompile(oneosc_ce)

    prob = ODEProblem(oneosc_ce_simpl, [], (0.0, 2.0))
    sol = solve(prob, Tsit5(), saveat = 0.1)

    @test typeof(oneosc_ce_simpl) == System
    @test sol(0.5, idxs = oscce.x) < 1.0 # test whether x(t) decreases over time
    @test sol(1.5, idxs = oscce.x) > 0.5 # test whether event happened
end

@testset "Additional SymbolicContinuousCallback options" begin
    # baseline affect (pos + neg + left root find)
    @variables c1(t)=1.0 c2(t)=1.0 # c1 = cos(t), c2 = cos(3t)
    eqs = [D(c1) ~ -sin(t); D(c2) ~ -3 * sin(3 * t)]
    function record_crossings(mod, obs, ctx, integ)
        push!(ctx, integ.t => obs.v)
        return (;)
    end
    cr1 = []
    cr2 = []
    evt1 = ModelingToolkit.SymbolicContinuousCallback(
        [c1 ~ 0], (f = record_crossings, observed = (; v = c1), ctx = cr1))
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (f = record_crossings, observed = (; v = c2), ctx = cr2))
    @named trigsys = System(eqs, t; continuous_events = [evt1, evt2])
    trigsys_ss = mtkcompile(trigsys)
    prob = ODEProblem(trigsys_ss, [], (0.0, 2π))
    sol = solve(prob, Tsit5())
    required_crossings_c1 = [π / 2, 3 * π / 2]
    required_crossings_c2 = [π / 6, π / 2, 5 * π / 6, 7 * π / 6, 3 * π / 2, 11 * π / 6]
    @test maximum(abs.(first.(cr1) .- required_crossings_c1)) < 1e-4
    @test maximum(abs.(first.(cr2) .- required_crossings_c2)) < 1e-4
    @test sign.(cos.(required_crossings_c1 .- 1e-6)) == sign.(last.(cr1))
    @test sign.(cos.(3 * (required_crossings_c2 .- 1e-6))) == sign.(last.(cr2))

    # with neg affect (pos * neg + left root find)
    cr1p = []
    cr2p = []
    cr1n = []
    cr2n = []
    evt1 = ModelingToolkit.SymbolicContinuousCallback(
        [c1 ~ 0], (f = record_crossings, observed = (; v = c1), ctx = cr1p);
        affect_neg = (f = record_crossings, observed = (; v = c1), ctx = cr1n))
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (f = record_crossings, observed = (; v = c2), ctx = cr2p);
        affect_neg = (f = record_crossings, observed = (; v = c2), ctx = cr2n))
    @named trigsys = System(eqs, t; continuous_events = [evt1, evt2])
    trigsys_ss = mtkcompile(trigsys)
    prob = ODEProblem(trigsys_ss, [], (0.0, 2π))
    sol = solve(prob, Tsit5(); dtmax = 0.01)
    c1_pc = filter((<=)(0) ∘ sin, required_crossings_c1)
    c1_nc = filter((>=)(0) ∘ sin, required_crossings_c1)
    c2_pc = filter(c -> -sin(3c) > 0, required_crossings_c2)
    c2_nc = filter(c -> -sin(3c) < 0, required_crossings_c2)
    @test maximum(abs.(c1_pc .- first.(cr1p))) < 1e-5
    @test maximum(abs.(c1_nc .- first.(cr1n))) < 1e-5
    @test maximum(abs.(c2_pc .- first.(cr2p))) < 1e-5
    @test maximum(abs.(c2_nc .- first.(cr2n))) < 1e-5
    @test sign.(cos.(c1_pc .- 1e-6)) == sign.(last.(cr1p))
    @test sign.(cos.(c1_nc .- 1e-6)) == sign.(last.(cr1n))
    @test sign.(cos.(3 * (c2_pc .- 1e-6))) == sign.(last.(cr2p))
    @test sign.(cos.(3 * (c2_nc .- 1e-6))) == sign.(last.(cr2n))

    # with nothing neg affect (pos * neg + left root find)
    cr1p = []
    cr2p = []
    evt1 = ModelingToolkit.SymbolicContinuousCallback(
        [c1 ~ 0], (f = record_crossings, observed = (; v = c1), ctx = cr1p); affect_neg = nothing)
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (f = record_crossings, observed = (; v = c2), ctx = cr2p); affect_neg = nothing)
    @named trigsys = System(eqs, t; continuous_events = [evt1, evt2])
    trigsys_ss = mtkcompile(trigsys)
    prob = ODEProblem(trigsys_ss, [], (0.0, 2π))
    sol = solve(prob, Tsit5(); dtmax = 0.01)
    @test maximum(abs.(c1_pc .- first.(cr1p))) < 1e-5
    @test maximum(abs.(c2_pc .- first.(cr2p))) < 1e-5
    @test sign.(cos.(c1_pc .- 1e-6)) == sign.(last.(cr1p))
    @test sign.(cos.(3 * (c2_pc .- 1e-6))) == sign.(last.(cr2p))

    #mixed
    cr1p = []
    cr2p = []
    cr1n = []
    cr2n = []
    evt1 = ModelingToolkit.SymbolicContinuousCallback(
        [c1 ~ 0], (f = record_crossings, observed = (; v = c1), ctx = cr1p); affect_neg = nothing)
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (f = record_crossings, observed = (; v = c2), ctx = cr2p);
        affect_neg = (f = record_crossings, observed = (; v = c2), ctx = cr2n))
    @named trigsys = System(eqs, t; continuous_events = [evt1, evt2])
    trigsys_ss = mtkcompile(trigsys)
    prob = ODEProblem(trigsys_ss, [], (0.0, 2π))
    sol = solve(prob, Tsit5(); dtmax = 0.01)
    c1_pc = filter((<=)(0) ∘ sin, required_crossings_c1)
    c2_pc = filter(c -> -sin(3c) > 0, required_crossings_c2)
    c2_nc = filter(c -> -sin(3c) < 0, required_crossings_c2)
    @test maximum(abs.(c1_pc .- first.(cr1p))) < 1e-5
    @test maximum(abs.(c2_pc .- first.(cr2p))) < 1e-5
    @test maximum(abs.(c2_nc .- first.(cr2n))) < 1e-5
    @test sign.(cos.(c1_pc .- 1e-6)) == sign.(last.(cr1p))
    @test sign.(cos.(3 * (c2_pc .- 1e-6))) == sign.(last.(cr2p))
    @test sign.(cos.(3 * (c2_nc .- 1e-6))) == sign.(last.(cr2n))

    # baseline affect w/ right rootfind (pos + neg + right root find)
    @variables c1(t)=1.0 c2(t)=1.0 # c1 = cos(t), c2 = cos(3t)
    cr1 = []
    cr2 = []
    evt1 = ModelingToolkit.SymbolicContinuousCallback(
        [c1 ~ 0], (f = record_crossings, observed = (; v = c1), ctx = cr1);
        rootfind = SciMLBase.RightRootFind)
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (f = record_crossings, observed = (; v = c2), ctx = cr2);
        rootfind = SciMLBase.RightRootFind)
    @named trigsys = System(eqs, t; continuous_events = [evt1, evt2])
    trigsys_ss = mtkcompile(trigsys)
    prob = ODEProblem(trigsys_ss, [], (0.0, 2π))
    sol = solve(prob, Tsit5(); dtmax = 0.01)
    required_crossings_c1 = [π / 2, 3 * π / 2]
    required_crossings_c2 = [π / 6, π / 2, 5 * π / 6, 7 * π / 6, 3 * π / 2, 11 * π / 6]
    @test maximum(abs.(first.(cr1) .- required_crossings_c1)) < 1e-4
    @test maximum(abs.(first.(cr2) .- required_crossings_c2)) < 1e-4
    @test sign.(cos.(required_crossings_c1 .+ 1e-6)) == sign.(last.(cr1))
    @test sign.(cos.(3 * (required_crossings_c2 .+ 1e-6))) == sign.(last.(cr2))

    # baseline affect w/ mixed rootfind (pos + neg + right root find)
    cr1 = []
    cr2 = []
    evt1 = ModelingToolkit.SymbolicContinuousCallback(
        [c1 ~ 0], (f = record_crossings, observed = (; v = c1), ctx = cr1);
        rootfind = SciMLBase.LeftRootFind)
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (f = record_crossings, observed = (; v = c2), ctx = cr2);
        rootfind = SciMLBase.RightRootFind)
    @named trigsys = System(eqs, t; continuous_events = [evt1, evt2])
    trigsys_ss = mtkcompile(trigsys)
    prob = ODEProblem(trigsys_ss, [], (0.0, 2π))
    sol = solve(prob, Tsit5())
    @test maximum(abs.(first.(cr1) .- required_crossings_c1)) < 1e-4
    @test maximum(abs.(first.(cr2) .- required_crossings_c2)) < 1e-4
    @test sign.(cos.(required_crossings_c1 .- 1e-6)) == sign.(last.(cr1))
    @test sign.(cos.(3 * (required_crossings_c2 .+ 1e-6))) == sign.(last.(cr2))

    #flip order and ensure results are okay
    cr1 = []
    cr2 = []
    evt1 = ModelingToolkit.SymbolicContinuousCallback(
        [c1 ~ 0], (f = record_crossings, observed = (; v = c1), ctx = cr1);
        rootfind = SciMLBase.LeftRootFind)
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (f = record_crossings, observed = (; v = c2), ctx = cr2);
        rootfind = SciMLBase.RightRootFind)
    @named trigsys = System(eqs, t; continuous_events = [evt2, evt1])
    trigsys_ss = mtkcompile(trigsys)
    prob = ODEProblem(trigsys_ss, [], (0.0, 2π))
    sol = solve(prob, Tsit5())
    @test maximum(abs.(first.(cr1) .- required_crossings_c1)) < 1e-4
    @test maximum(abs.(first.(cr2) .- required_crossings_c2)) < 1e-4
    @test sign.(cos.(required_crossings_c1 .- 1e-6)) == sign.(last.(cr1))
    @test sign.(cos.(3 * (required_crossings_c2 .+ 1e-6))) == sign.(last.(cr2))
end

@testset "Discrete event reinitialization (#3142)" begin
    @connector LiquidPort begin
        p(t)::Float64, [description = "Set pressure in bar",
            guess = 1.01325]
        Vdot(t)::Float64,
        [description = "Volume flow rate in L/min",
            guess = 0.0,
            connect = Flow]
    end

    @mtkmodel PressureSource begin
        @components begin
            port = LiquidPort()
        end
        @parameters begin
            p_set::Float64 = 1.01325, [description = "Set pressure in bar"]
        end
        @equations begin
            port.p ~ p_set
        end
    end

    @mtkmodel BinaryValve begin
        @constants begin
            p_ref::Float64 = 1.0, [description = "Reference pressure drop in bar"]
            ρ_ref::Float64 = 1000.0, [description = "Reference density in kg/m^3"]
        end
        @components begin
            port_in = LiquidPort()
            port_out = LiquidPort()
        end
        @parameters begin
            k_V::Float64 = 1.0, [description = "Valve coefficient in L/min/bar"]
            k_leakage::Float64 = 1e-08, [description = "Leakage coefficient in L/min/bar"]
            ρ::Float64 = 1000.0, [description = "Density in kg/m^3"]
        end
        @variables begin
            S(t)::Float64, [description = "Valve state", guess = 1.0, irreducible = true]
            Δp(t)::Float64, [description = "Pressure difference in bar", guess = 1.0]
            Vdot(t)::Float64, [description = "Volume flow rate in L/min", guess = 1.0]
        end
        @equations begin
            # Port handling
            port_in.Vdot ~ -Vdot
            port_out.Vdot ~ Vdot
            Δp ~ port_in.p - port_out.p
            # System behavior
            D(S) ~ 0.0
            Vdot ~ S * k_V * sign(Δp) * sqrt(abs(Δp) / p_ref * ρ_ref / ρ) + k_leakage * Δp # softplus alpha function to avoid negative values under the sqrt
        end
    end

    # Test System
    @mtkmodel TestSystem begin
        @components begin
            pressure_source_1 = PressureSource(p_set = 2.0)
            binary_valve_1 = BinaryValve(S = 1.0, k_leakage = 0.0)
            binary_valve_2 = BinaryValve(S = 1.0, k_leakage = 0.0)
            pressure_source_2 = PressureSource(p_set = 1.0)
        end
        @equations begin
            connect(pressure_source_1.port, binary_valve_1.port_in)
            connect(binary_valve_1.port_out, binary_valve_2.port_in)
            connect(binary_valve_2.port_out, pressure_source_2.port)
        end
        @discrete_events begin
            [30] => [binary_valve_1.S ~ 0.0, binary_valve_2.Δp ~ 0.0]
            [60] => [binary_valve_1.S ~ 1.0, binary_valve_2.Δp ~ 1.0]
            [120] => [binary_valve_1.S ~ 0.0, binary_valve_2.Δp ~ 0.0]
        end
    end

    # Test Simulation
    @mtkcompile sys = TestSystem()

    # Test Simulation
    prob = ODEProblem(sys, [], (0.0, 150.0))
    sol = solve(prob)
    # This is singular at the second event, but the derivatives are zero so it's
    # constant after that point anyway. Just make sure it hits the last event and
    # had the correct `u`.
    @test_broken SciMLBase.successful_retcode(sol)
    @test sol.t[end] >= 120.0
    @test sol[end] == [0.0, 0.0, 0.0]
end

@testset "Discrete variable timeseries" begin
    @variables x(t)
    @parameters a(t) b(t) c(t)
    cb1 = SymbolicContinuousCallback([x ~ 1.0] => [a ~ -Pre(a)], discrete_parameters = [a])
    function save_affect!(mod, obs, ctx, integ)
        return (; b = 5.0)
    end
    cb2 = [x ~ 0.5] => (f = save_affect!, modified = (; b))
    cb3 = SymbolicDiscreteCallback(1.0 => [c ~ t], discrete_parameters = [c])

    @mtkcompile sys = System(D(x) ~ cos(t), t, [x], [a, b, c];
        continuous_events = [cb1, cb2], discrete_events = [cb3])
    prob = ODEProblem(sys, [x => 1.0, a => 1.0, b => 2.0, c => 0.0], (0.0, 2pi))
    @test sort(canonicalize(Discrete(), prob.p)[1]) == [0.0, 1.0, 2.0]
    sol = solve(prob, Tsit5())

    @test sol[a] == [1.0, -1.0]
    @test sol[b] == [2.0, 5.0, 5.0]
    @test sol[c] == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
end

@testset "Heater" begin
    @variables temp(t)
    params = @parameters furnace_on_threshold=0.5 furnace_off_threshold=0.7 furnace_power=1.0 leakage=0.1 furnace_on::Bool=false
    eqs = [
        D(temp) ~ furnace_on * furnace_power - temp^2 * leakage
    ]

    furnace_off = ModelingToolkit.SymbolicContinuousCallback(
        [temp ~ furnace_off_threshold],
        ModelingToolkit.ImperativeAffect(modified = (; furnace_on)) do x, o, i, c
            @set! x.furnace_on = false
        end)
    furnace_enable = ModelingToolkit.SymbolicContinuousCallback(
        [temp ~ furnace_on_threshold],
        ModelingToolkit.ImperativeAffect(modified = (; furnace_on)) do x, o, i, c
            @set! x.furnace_on = true
        end)
    @named sys = System(
        eqs, t, [temp], params; continuous_events = [furnace_off, furnace_enable])
    ss = mtkcompile(sys)
    prob = ODEProblem(ss, [temp => 0.0, furnace_on => true], (0.0, 100.0))
    sol = solve(prob, Tsit5(); dtmax = 0.01)
    @test all(sol[temp][sol.t .> 1.0] .<= 0.79) && all(sol[temp][sol.t .> 1.0] .>= 0.49)

    furnace_off = ModelingToolkit.SymbolicContinuousCallback(
        [temp ~ furnace_off_threshold],
        ModelingToolkit.ImperativeAffect(modified = (; furnace_on)) do x, o, c, i
            @set! x.furnace_on = false
        end; initialize = ModelingToolkit.ImperativeAffect(modified = (;
            temp)) do x, o, c, i
            @set! x.temp = 0.2
        end)
    furnace_enable = ModelingToolkit.SymbolicContinuousCallback(
        [temp ~ furnace_on_threshold],
        ModelingToolkit.ImperativeAffect(modified = (; furnace_on)) do x, o, c, i
            @set! x.furnace_on = true
        end)
    @named sys = System(
        eqs, t, [temp], params; continuous_events = [furnace_off, furnace_enable])
    ss = mtkcompile(sys)
    prob = ODEProblem(ss, [temp => 0.0, furnace_on => true], (0.0, 100.0))
    sol = solve(prob, Tsit5(); dtmax = 0.01)
    @test all(sol[temp][sol.t .> 1.0] .<= 0.79) && all(sol[temp][sol.t .> 1.0] .>= 0.49)
    @test all(sol[temp][sol.t .!= 0.0] .<= 0.79) && all(sol[temp][sol.t .!= 0.0] .>= 0.2)
end

@testset "ImperativeAffect errors and warnings" begin
    @variables temp(t)
    params = @parameters furnace_on_threshold=0.5 furnace_off_threshold=0.7 furnace_power=1.0 leakage=0.1 furnace_on::Bool=false
    eqs = [
        D(temp) ~ furnace_on * furnace_power - temp^2 * leakage
    ]

    furnace_off = ModelingToolkit.SymbolicContinuousCallback(
        [temp ~ furnace_off_threshold],
        ModelingToolkit.ImperativeAffect(
            modified = (; furnace_on), observed = (; furnace_on)) do x, o, c, i
            @set! x.furnace_on = false
        end)
    @named sys = System(eqs, t, [temp], params; continuous_events = [furnace_off])
    ss = mtkcompile(sys)
    @test_logs (:warn,
        "The symbols Any[:furnace_on] are declared as both observed and modified; this is a code smell because it becomes easy to confuse them and assign/not assign a value.") prob=ODEProblem(
        ss, [temp => 0.0, furnace_on => true], (0.0, 100.0))

    @variables tempsq(t) # trivially eliminated
    eqs = [tempsq ~ temp^2
           D(temp) ~ furnace_on * furnace_power - temp^2 * leakage]

    furnace_off = ModelingToolkit.SymbolicContinuousCallback(
        [temp ~ furnace_off_threshold],
        ModelingToolkit.ImperativeAffect(
            modified = (; furnace_on, tempsq), observed = (; furnace_on)) do x, o, c, i
            @set! x.furnace_on = false
        end)
    @named sys = System(
        eqs, t, [temp, tempsq], params; continuous_events = [furnace_off])
    ss = mtkcompile(sys)
    @test_throws "refers to missing variable(s)" prob=ODEProblem(
        ss, [temp => 0.0, furnace_on => true], (0.0, 100.0))

    @parameters not_actually_here
    furnace_off = ModelingToolkit.SymbolicContinuousCallback(
        [temp ~ furnace_off_threshold],
        ModelingToolkit.ImperativeAffect(modified = (; furnace_on),
            observed = (; furnace_on, not_actually_here)) do x, o, c, i
            @set! x.furnace_on = false
        end)
    @named sys = System(
        eqs, t, [temp, tempsq], params; continuous_events = [furnace_off])
    ss = mtkcompile(sys)
    @test_throws "refers to missing variable(s)" prob=ODEProblem(
        ss, [temp => 0.0, furnace_on => true], (0.0, 100.0))

    furnace_off = ModelingToolkit.SymbolicContinuousCallback(
        [temp ~ furnace_off_threshold],
        ModelingToolkit.ImperativeAffect(modified = (; furnace_on),
            observed = (; furnace_on)) do x, o, c, i
            return (; fictional2 = false)
        end)
    @named sys = System(
        eqs, t, [temp, tempsq], params; continuous_events = [furnace_off])
    ss = mtkcompile(sys)
    prob = ODEProblem(
        ss, [temp => 0.0, furnace_on => true], (0.0, 100.0))
    @test_throws "Tried to write back to" solve(prob, Tsit5())
end

@testset "Quadrature" begin
    @variables theta(t) omega(t)
    params = @parameters qA=0 qB=0 hA=0 hB=0 cnt::Int=0
    eqs = [D(theta) ~ omega
           omega ~ 1.0]
    function decoder(oldA, oldB, newA, newB)
        state = (oldA, oldB, newA, newB)
        if state == (0, 0, 1, 0) || state == (1, 0, 1, 1) || state == (1, 1, 0, 1) ||
           state == (0, 1, 0, 0)
            return 1
        elseif state == (0, 0, 0, 1) || state == (0, 1, 1, 1) || state == (1, 1, 1, 0) ||
               state == (1, 0, 0, 0)
            return -1
        elseif state == (0, 0, 0, 0) || state == (0, 1, 0, 1) || state == (1, 0, 1, 0) ||
               state == (1, 1, 1, 1)
            return 0
        else
            return 0 # err is interpreted as no movement
        end
    end
    qAevt = ModelingToolkit.SymbolicContinuousCallback([cos(100 * theta) ~ 0],
        ModelingToolkit.ImperativeAffect((; qA, hA, hB, cnt), (; qB)) do x, o, c, i
            @set! x.hA = x.qA
            @set! x.hB = o.qB
            @set! x.qA = 1
            @set! x.cnt += decoder(x.hA, x.hB, x.qA, o.qB)
            x
        end,
        affect_neg = ModelingToolkit.ImperativeAffect(
            (; qA, hA, hB, cnt), (; qB)) do x, o, c, i
            @set! x.hA = x.qA
            @set! x.hB = o.qB
            @set! x.qA = 0
            @set! x.cnt += decoder(x.hA, x.hB, x.qA, o.qB)
            x
        end; rootfind = SciMLBase.RightRootFind)
    qBevt = ModelingToolkit.SymbolicContinuousCallback([cos(100 * theta - π / 2) ~ 0],
        ModelingToolkit.ImperativeAffect((; qB, hA, hB, cnt), (; qA)) do x, o, c, i
            @set! x.hA = o.qA
            @set! x.hB = x.qB
            @set! x.qB = 1
            @set! x.cnt += decoder(x.hA, x.hB, o.qA, x.qB)
            x
        end,
        affect_neg = ModelingToolkit.ImperativeAffect(
            (; qB, hA, hB, cnt), (; qA)) do x, o, c, i
            @set! x.hA = o.qA
            @set! x.hB = x.qB
            @set! x.qB = 0
            @set! x.cnt += decoder(x.hA, x.hB, o.qA, x.qB)
            x
        end; rootfind = SciMLBase.RightRootFind)
    @named sys = System(
        eqs, t, [theta, omega], params; continuous_events = [qAevt, qBevt])
    ss = mtkcompile(sys)
    prob = ODEProblem(ss, [theta => 1e-5], (0.0, pi))
    sol = solve(prob, Tsit5(); dtmax = 0.01)
    @test getp(sol, cnt)(sol) == 198 # we get 2 pulses per phase cycle (cos 0 crossing) and we go to 100 cycles; we miss a few due to the initial state
end

@testset "Initialization" begin
    @variables x(t)
    seen = false
    f = ModelingToolkit.ImperativeAffect(f = (m, o, ctx, int) -> (seen = true; return (;)))
    cb1 = ModelingToolkit.SymbolicContinuousCallback(
        [x ~ 0], nothing, initialize = [x ~ 1.5], finalize = f)
    @mtkcompile sys = System(D(x) ~ -1, t, [x], []; continuous_events = [cb1])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 2))
    sol = solve(prob, Tsit5(); dtmax = 0.01)
    @test sol[x][1] ≈ 1.0
    @test sol[x][2] ≈ 1.5 # the initialize affect has been applied
    @test seen == true

    @variables x(t)
    seen = false
    f = ModelingToolkit.ImperativeAffect(f = (m, o, ctx, int) -> (seen = true; return (;)))
    cb1 = ModelingToolkit.SymbolicContinuousCallback(
        [x ~ 0], nothing, initialize = [x ~ 1.5], finalize = f)
    inited = false
    finaled = false
    a = ModelingToolkit.ImperativeAffect(f = (
        m, o, ctx, int) -> (inited = true; return (;)))
    b = ModelingToolkit.ImperativeAffect(f = (
        m, o, ctx, int) -> (finaled = true; return (;)))
    cb2 = ModelingToolkit.SymbolicContinuousCallback(
        [x ~ 0.1], nothing, initialize = a, finalize = b)
    @mtkcompile sys = System(D(x) ~ -1, t, [x], []; continuous_events = [cb1, cb2])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 2))
    sol = solve(prob, Tsit5())
    @test sol[x][1] ≈ 1.0
    @test sol[x][2] ≈ 1.5 # the initialize affect has been applied
    @test seen == true
    @test inited == true
    @test finaled == true

    #periodic
    inited = false
    finaled = false
    cb3 = ModelingToolkit.SymbolicDiscreteCallback(
        1.0, [x ~ 2], initialize = a, finalize = b)
    @mtkcompile sys = System(D(x) ~ -1, t, [x], []; discrete_events = [cb3])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 2))
    sol = solve(prob, Tsit5())
    @test inited == true
    @test finaled == true
    @test isapprox(sol[x][3], 0.0, atol = 1e-9)
    @test sol[x][4] ≈ 2.0
    @test sol[x][5] ≈ 1.0

    seen = false
    inited = false
    finaled = false
    cb3 = ModelingToolkit.SymbolicDiscreteCallback(1.0, f, initialize = a, finalize = b)
    @mtkcompile sys = System(D(x) ~ -1, t, [x], []; discrete_events = [cb3])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 2))
    sol = solve(prob, Tsit5())
    @test seen == true
    @test inited == true

    #preset
    seen = false
    inited = false
    finaled = false
    cb3 = ModelingToolkit.SymbolicDiscreteCallback([1.0], f, initialize = a, finalize = b)
    @mtkcompile sys = System(D(x) ~ -1, t, [x], []; discrete_events = [cb3])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 2))
    sol = solve(prob, Tsit5())
    @test seen == true
    @test inited == true
    @test finaled == true

    #equational
    seen = false
    inited = false
    finaled = false
    cb3 = ModelingToolkit.SymbolicDiscreteCallback(
        t == 1.0, f, initialize = a, finalize = b)
    @mtkcompile sys = System(D(x) ~ -1, t, [x], []; discrete_events = [cb3])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 2))
    sol = solve(prob, Tsit5(); tstops = 1.0)
    @test seen == true
    @test inited == true
    @test finaled == true
end

@testset "Bump" begin
    @variables x(t) [irreducible = true] y(t) [irreducible = true]
    eqs = [x ~ y, D(x) ~ -1]
    cb = [x ~ 0.0] => [x ~ 0, y ~ 1]
    @mtkcompile pend = System(eqs, t; continuous_events = [cb])
    prob = ODEProblem(pend, [x => 1], (0.0, 3.0), guesses = [y => x])
    @test_broken !SciMLBase.successful_retcode(solve(prob, Rodas5()))

    cb = [x ~ 0.0] => [y ~ 1]
    @mtkcompile pend = System(eqs, t; continuous_events = [cb])
    prob = ODEProblem(pend, [x => 1], (0.0, 3.0), guesses = [y => x])
    @test_broken !SciMLBase.successful_retcode(solve(prob, Rodas5()))

    cb = [x ~ 0.0] => [x ~ 1, y ~ 1]
    @mtkcompile pend = System(eqs, t; continuous_events = [cb])
    prob = ODEProblem(pend, [x => 1], (0.0, 3.0), guesses = [y => x])
    @test all(≈(0.0; atol = 1e-9), solve(prob, Rodas5())[[x, y]][end])
end

@testset "Issue#3154 Array variable in discrete condition" begin
    @mtkmodel DECAY begin
        @parameters begin
            unrelated[1:2] = zeros(2)
            k(t) = 0.0
        end
        @variables begin
            x(t) = 10.0
        end
        @equations begin
            D(x) ~ -k * x
        end
        @discrete_events begin
            (t == 1.0) => [k ~ 1.0], [discrete_parameters = k]
        end
    end
    @mtkcompile decay = DECAY()
    prob = ODEProblem(decay, [], (0.0, 10.0))
    @test_nowarn solve(prob, Tsit5(), tstops = [1.0])
end

@testset "Array parameter updates in ImperativeAffect" begin
    function weird1(max_time; name)
        params = @parameters begin
            θ(t) = 0.0
        end
        vars = @variables begin
            x(t) = 0.0
        end
        eqs = reduce(vcat, Symbolics.scalarize.([
            D(x) ~ 1.0
        ]))
        reset = ModelingToolkit.ImperativeAffect(
            modified = (; x, θ)) do m, o, _, i
            @set! m.θ = 0.0
            @set! m.x = 0.0
            return m
        end
        return System(eqs, t, vars, params; name = name,
            continuous_events = [[x ~ max_time] => reset])
    end

    function weird2(max_time; name)
        params = @parameters begin
            θ(t) = 0.0
        end
        vars = @variables begin
            x(t) = 0.0
        end
        eqs = reduce(vcat, Symbolics.scalarize.([
            D(x) ~ 1.0
        ]))
        return System(eqs, t, vars, params; name = name) # note no event
    end

    @named wd1 = weird1(0.021)
    @named wd2 = weird2(0.021)

    sys1 = mtkcompile(System(Equation[], t; name = :parent,
        discrete_events = [0.01 => ModelingToolkit.ImperativeAffect(
            modified = (; θs = reduce(vcat, [[wd1.θ]])), ctx = [1]) do m, o, c, i
            @set! m.θs[1] = c[] += 1
        end],
        systems = [wd1]))
    sys2 = mtkcompile(System(Equation[], t; name = :parent,
        discrete_events = [0.01 => ModelingToolkit.ImperativeAffect(
            modified = (; θs = reduce(vcat, [[wd2.θ]])), ctx = [1]) do m, o, c, i
            @set! m.θs[1] = c[] += 1
        end],
        systems = [wd2]))

    sol1 = solve(ODEProblem(sys1, [], (0.0, 1.0)), Tsit5())
    @test 100.0 ∈ sol1[sys1.wd1.θ]

    sol2 = solve(ODEProblem(sys2, [], (0.0, 1.0)), Tsit5())
    @test 100.0 ∈ sol2[sys2.wd2.θ]
end

@testset "Implicit affects with Pre" begin
    using ModelingToolkit: UnsolvableCallbackError
    @parameters g
    @variables x(t) y(t) λ(t)
    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ 1]
    c_evt = [t ~ 5.0] => [x ~ Pre(x) + 0.1]
    @mtkcompile pend = System(eqs, t, continuous_events = c_evt)
    prob = ODEProblem(pend, [x => -1, y => 0, g => 1], (0.0, 10.0), guesses = [λ => 1])
    sol = solve(prob, FBDF())
    @test ≈(sol(5.000001, idxs = x) - sol(4.999999, idxs = x), 0.1, rtol = 1e-4)
    @test ≈(sol(5.000001, idxs = x)^2 + sol(5.000001, idxs = y)^2, 1, rtol = 1e-4)

    # Implicit affect with Pre
    c_evt = [t ~ 5.0] => [x ~ Pre(x) + y^2]
    @mtkcompile pend = System(eqs, t, continuous_events = c_evt)
    prob = ODEProblem(pend, [x => 1, y => 0, g => 1], (0.0, 10.0), guesses = [λ => 1])
    sol = solve(prob, FBDF())
    @test ≈(sol(5.000001, idxs = y)^2 + sol(4.999999, idxs = x),
        sol(5.000001, idxs = x), rtol = 1e-4)
    @test ≈(sol(5.000001, idxs = x)^2 + sol(5.000001, idxs = y)^2, 1, rtol = 1e-4)

    # Impossible affect errors
    c_evt = [t ~ 5.0] => [x ~ Pre(x) + 2]
    @mtkcompile pend = System(eqs, t, continuous_events = c_evt)
    prob = ODEProblem(pend, [x => 1, y => 0, g => 1], (0.0, 10.0), guesses = [λ => 1])
    @test_throws UnsolvableCallbackError sol=solve(prob, FBDF())

    # Changing both variables and parameters in the same affect.
    @parameters g(t)
    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ 1]
    c_evt = SymbolicContinuousCallback(
        [t ~ 5.0], [x ~ Pre(x) + 0.1, g ~ Pre(g) + 1], discrete_parameters = [g], iv = t)
    @mtkcompile pend = System(eqs, t, continuous_events = c_evt)
    prob = ODEProblem(pend, [x => 1, y => 0, g => 1], (0.0, 10.0), guesses = [λ => 1])
    sol = solve(prob, FBDF())
    @test sol.ps[g] ≈ [1, 2]
    @test ≈(sol(5.0000001, idxs = x) - sol(4.999999, idxs = x), 0.1, rtol = 1e-4)

    # Proper re-initialization after parameter change
    eqs = [y ~ g^2, D(x) ~ x]
    c_evt = SymbolicContinuousCallback(
        [t ~ 5.0], [x ~ Pre(x) + 1, g ~ Pre(g) + 1], discrete_parameters = [g], iv = t)
    @mtkcompile sys = System(eqs, t, continuous_events = c_evt)
    prob = ODEProblem(sys, [x => 1.0, g => 2], (0.0, 10.0))
    sol = solve(prob, FBDF())
    @test sol.ps[g] ≈ [2.0, 3.0]
    @test ≈(sol(5.00000001, idxs = x) - sol(4.9999999, idxs = x), 1; rtol = 1e-4)
    @test ≈(sol(5.00000001, idxs = y), 9, rtol = 1e-4)

    # Parameters that don't appear in affects should not be mutated.
    c_evt = [t ~ 5.0] => [x ~ Pre(x) + 1]
    @mtkcompile sys = System(eqs, t, continuous_events = c_evt)
    prob = ODEProblem(sys, [x => 0.5, g => 2], (0.0, 10.0), guesses = [y => 0])
    sol = solve(prob, FBDF())
    @test prob.ps[g] == sol.ps[g]
end

@testset "Array parameter updates of parent components in ImperativeEffect" begin
    function child(vals; name, max_time = 0.1)
        vars = @variables begin
            x(t) = 0.0
        end
        eqs = reduce(vcat, Symbolics.scalarize.([
            D(x) ~ 1.0
        ]))
        reset = ModelingToolkit.ImperativeAffect(
            modified = (; vals = Symbolics.scalarize(ParentScope.(vals)), x)) do m, o, _, i
            @set! m.vals = m.vals .+ 1
            @set! m.x = 0.0
            return m
        end
        return System(eqs, t, vars, []; name = name,
            continuous_events = [[x ~ max_time] => reset])
    end
    shared_pars = @parameters begin
        vals(t)[1:2] = 0.0
    end

    @named sys = System(Equation[], t, [], Symbolics.scalarize(vals);
        systems = [child(vals; name = :child)])
    sys = mtkcompile(sys)
    sol = solve(ODEProblem(sys, [], (0.0, 1.0)), Tsit5())
end

@testset "non-floating-point discretes and namespaced affects" begin
    function Inner(; name)
        @parameters p(t)::Int
        @variables x(t)
        cevs = ModelingToolkit.SymbolicContinuousCallback(
            [x ~ 1.0], [p ~ Pre(p) + 1]; iv = t, discrete_parameters = [p])
        System([D(x) ~ 1], t, [x], [p]; continuous_events = [cevs], name)
    end
    @named inner = Inner()
    @mtkcompile sys = System(Equation[], t; systems = [inner])
    prob = ODEProblem(sys, [inner.x => 0.0, inner.p => 0], (0.0, 5.0))
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
    @test sol[inner.p][end] ≈ 1.0
end

mutable struct ParamTest
    y::Any
end

@testset "callable parameter and symbolic affect" begin
    (pt::ParamTest)(x) = pt.y - x

    p1 = ParamTest(1)
    tp1 = typeof(p1)
    @parameters (p_1::tp1)(..) = p1
    @parameters p2(t) = 1.0
    @variables x(t) = 0.0
    @variables x2(t)
    event = [0.5] => [p2 ~ Pre(t)]

    eq = [
        D(x) ~ p2,
        x2 ~ p_1(x)
    ]
    @mtkcompile sys = System(eq, t, [x, x2], [p_1, p2], discrete_events = [event])

    prob = ODEProblem(sys, [], (0.0, 1.0))
    sol = solve(prob)
    @test SciMLBase.successful_retcode(sol)
    @test sol[x, end]≈1.0 atol=1e-6
end

@testset "Symbolic affects are compiled in `complete`" begin
    @parameters g
    @variables x(t) [state_priority = 10.0] y(t) [guess = 1.0]
    @variables λ(t) [guess = 1.0]
    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ 1]
    cevts = [[x ~ 0.0] => [D(x) ~ Pre(D(x)) + 1sign(Pre(D(x)))]]
    @named pend = System(eqs, t; continuous_events = cevts)

    scc = only(continuous_events(pend))
    @test scc.affect isa ModelingToolkit.SymbolicAffect

    pend = mtkcompile(pend)

    scc = only(continuous_events(pend))
    @test scc.affect isa ModelingToolkit.AffectSystem
    @test length(ModelingToolkit.all_equations(scc.affect)) == 5 # 1 affect, 3 algebraic, 1 observed

    u0 = [x => -1/2, D(x) => 1/2, g => 1]
    prob = ODEProblem(pend, u0, (0.0, 5.0))
    sol = solve(prob, FBDF())
    @test SciMLBase.successful_retcode(sol)
end

@testset "Algebraic equation with input variable in symbolic affect" begin
    # Specifically happens when the variable marked as an input is an algebraic variable
    # in the affect system.
    @variables x(t) [input = true] y(t)
    dev = ModelingToolkit.SymbolicDiscreteCallback(1.0, [y ~ Pre(y) + 1])
    @named sys = System([D(y) ~ 2x + 1, x^2 ~ 2y^3], t; discrete_events = [dev])
    sys = @test_nowarn mtkcompile(sys)
end

@testset "Non-`Real` symtype parameters in callback with unknown" begin
    @mtkmodel MWE begin
        @variables begin
            x(t) = 1.0
        end
        @parameters begin
            p(t)::Int = 1
        end
        @equations begin
            D(x) ~ p * x
        end
        @discrete_events begin
            1.0 => [x ~ p * Pre(x) * sin(x)]
        end
    end
    @mtkcompile sys = MWE()
    @test_nowarn ODEProblem(sys, [], (0.0, 1.0))
end

@testset "Test erroneously created events yields errors" begin
    @parameters p(t) d
    @variables X(t)
    @test_throws "Vectors of symbolic conditions are not allowed" SymbolicDiscreteCallback([X <
                                                                                            5.0] => [X ~
                                                                                                     Pre(X) +
                                                                                                     10.0])
    @test_throws "Vectors of symbolic conditions are not allowed" SymbolicDiscreteCallback([
        X < 5.0, X > 10.0] => [X ~ Pre(X) + 10.0])
    @test_throws "MethodError: no method" SymbolicContinuousCallback((X <
                                                                      5.0) => [X ~
                                                                               Pre(X) +
                                                                               10.0])
end

@testset "Issue#3990: Scalarized array passed to `discrete_parameters` of symbolic affect" begin
    N = 2
    @parameters v(t)[1:N]
    @parameters M(t)[1:N, 1:N]

    @variables x(t)

    Mini = rand(N, N) ./ (N^2)
    vini = vec(sum(Mini, dims = 1))

    v_eq = [D(x) ~ x * Symbolics.scalarize(sum(v))]
    M_eq = [D(x) ~ x * Symbolics.scalarize(sum(M))]

    v_event = ModelingToolkit.SymbolicDiscreteCallback(
        1.0,
        [v ~ -Pre(v)],
        discrete_parameters = [v]
    )

    M_event = ModelingToolkit.SymbolicDiscreteCallback(
        1.0,
        [M ~ -Pre(M)],
        discrete_parameters = [M]
    )

    @mtkcompile v_sys = System(v_eq, t; discrete_events = v_event)
    @mtkcompile M_sys = System(M_eq, t; discrete_events = M_event)

    u0p0_map = Dict(x => 1.0, M => Mini, v => vini)

    v_prob = ODEProblem(v_sys, u0p0_map, (0.0, 2.5))
    M_prob = ODEProblem(M_sys, u0p0_map, (0.0, 2.5))

    v_sol = solve(v_prob, Tsit5())
    M_sol = solve(M_prob, Tsit5())

    @test v_sol[v] ≈ [vini, -vini, vini]
    @test M_sol[M] ≈ [Mini, -Mini, Mini]
end

@testset "Issue#3990: Scalarized array passed to `discrete_parameters` of symbolic affect" begin
    N = 2
    @parameters v(t)[1:N]
    @parameters M(t)[1:N, 1:N]

    @variables x(t)

    Mini = rand(N, N) ./ (N^2)
    vini = vec(sum(Mini, dims = 1))

    v_eq = [D(x) ~ x * Symbolics.scalarize(sum(v))]
    M_eq = [D(x) ~ x * Symbolics.scalarize(sum(M))]

    v_event = ModelingToolkit.SymbolicDiscreteCallback(
        1.0,
        [v ~ -Pre(v)],
        discrete_parameters = collect(v)
    )

    M_event = ModelingToolkit.SymbolicDiscreteCallback(
        1.0,
        [M ~ -Pre(M)],
        discrete_parameters = vec(collect(M))
    )

    @mtkcompile v_sys = System(v_eq, t; discrete_events = v_event)
    @mtkcompile M_sys = System(M_eq, t; discrete_events = M_event)

    u0p0_map = Dict(x => 1.0, M => Mini, v => vini)

    v_prob = ODEProblem(v_sys, u0p0_map, (0.0, 2.5))
    M_prob = ODEProblem(M_sys, u0p0_map, (0.0, 2.5))

    v_sol = solve(v_prob, Tsit5())
    M_sol = solve(M_prob, Tsit5())

    @test v_sol[v] ≈ [vini, -vini, vini]
    @test M_sol[M] ≈ [Mini, -Mini, Mini]
end
