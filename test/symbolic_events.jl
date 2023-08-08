using ModelingToolkit, OrdinaryDiffEq, StochasticDiffEq, JumpProcesses, Test
using ModelingToolkit: SymbolicContinuousCallback,
    SymbolicContinuousCallbacks, NULL_AFFECT,
    get_callback
using StableRNGs
rng = StableRNG(12345)

@parameters t
@variables x(t) = 0
D = Differential(t)

eqs = [D(x) ~ 1]
affect = [x ~ 0]

## Test SymbolicContinuousCallback
@testset "SymbolicContinuousCallback constructors" begin
    e = SymbolicContinuousCallback(eqs[])
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT

    e = SymbolicContinuousCallback(eqs)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT

    e = SymbolicContinuousCallback(eqs, NULL_AFFECT)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT

    e = SymbolicContinuousCallback(eqs[], NULL_AFFECT)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT

    e = SymbolicContinuousCallback(eqs => NULL_AFFECT)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT

    e = SymbolicContinuousCallback(eqs[] => NULL_AFFECT)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT

    ## With affect

    e = SymbolicContinuousCallback(eqs[], affect)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect

    e = SymbolicContinuousCallback(eqs, affect)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect

    e = SymbolicContinuousCallback(eqs, affect)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect

    e = SymbolicContinuousCallback(eqs[], affect)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect

    e = SymbolicContinuousCallback(eqs => affect)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect

    e = SymbolicContinuousCallback(eqs[] => affect)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect

    # test plural constructor

    e = SymbolicContinuousCallbacks(eqs[])
    @test e isa Vector{SymbolicContinuousCallback}
    @test isequal(e[].eqs, eqs)
    @test e[].affect == NULL_AFFECT

    e = SymbolicContinuousCallbacks(eqs)
    @test e isa Vector{SymbolicContinuousCallback}
    @test isequal(e[].eqs, eqs)
    @test e[].affect == NULL_AFFECT

    e = SymbolicContinuousCallbacks(eqs[] => affect)
    @test e isa Vector{SymbolicContinuousCallback}
    @test isequal(e[].eqs, eqs)
    @test e[].affect == affect

    e = SymbolicContinuousCallbacks(eqs => affect)
    @test e isa Vector{SymbolicContinuousCallback}
    @test isequal(e[].eqs, eqs)
    @test e[].affect == affect

    e = SymbolicContinuousCallbacks([eqs[] => affect])
    @test e isa Vector{SymbolicContinuousCallback}
    @test isequal(e[].eqs, eqs)
    @test e[].affect == affect

    e = SymbolicContinuousCallbacks([eqs => affect])
    @test e isa Vector{SymbolicContinuousCallback}
    @test isequal(e[].eqs, eqs)
    @test e[].affect == affect

    e = SymbolicContinuousCallbacks(SymbolicContinuousCallbacks([eqs => affect]))
    @test e isa Vector{SymbolicContinuousCallback}
    @test isequal(e[].eqs, eqs)
    @test e[].affect == affect
end

##

@named sys = ODESystem(eqs, continuous_events = [x ~ 1])
@test getfield(sys, :continuous_events)[] ==
      SymbolicContinuousCallback(Equation[x ~ 1], NULL_AFFECT)
@test isequal(equations(getfield(sys, :continuous_events))[], x ~ 1)
fsys = flatten(sys)
@test isequal(equations(getfield(fsys, :continuous_events))[], x ~ 1)

@named sys2 = ODESystem([D(x) ~ 1], continuous_events = [x ~ 2], systems = [sys])
@test getfield(sys2, :continuous_events)[] ==
      SymbolicContinuousCallback(Equation[x ~ 2], NULL_AFFECT)
@test all(ModelingToolkit.continuous_events(sys2) .== [
    SymbolicContinuousCallback(Equation[x ~ 2], NULL_AFFECT),
    SymbolicContinuousCallback(Equation[sys.x ~ 1], NULL_AFFECT),
])

@test isequal(equations(getfield(sys2, :continuous_events))[1], x ~ 2)
@test length(ModelingToolkit.continuous_events(sys2)) == 2
@test isequal(ModelingToolkit.continuous_events(sys2)[1].eqs[], x ~ 2)
@test isequal(ModelingToolkit.continuous_events(sys2)[2].eqs[], sys.x ~ 1)

# Functions should be generated for root-finding equations
prob = ODEProblem(sys, Pair[], (0.0, 2.0))
p0 = 0
t0 = 0
@test get_callback(prob) isa ModelingToolkit.DiffEqCallbacks.ContinuousCallback
cb = ModelingToolkit.generate_rootfinding_callback(sys)
cond = cb.condition
out = [0.0]
cond.rf_ip(out, [0], p0, t0)
@test out[] ≈ -1 # signature is u,p,t
cond.rf_ip(out, [1], p0, t0)
@test out[] ≈ 0  # signature is u,p,t
cond.rf_ip(out, [2], p0, t0)
@test out[] ≈ 1  # signature is u,p,t

prob = ODEProblem(sys, Pair[], (0.0, 2.0))
sol = solve(prob, Tsit5())
@test minimum(t -> abs(t - 1), sol.t) < 1e-10 # test that the solver stepped at the root

# Test that a user provided callback is respected
test_callback = DiscreteCallback(x -> x, x -> x)
prob = ODEProblem(sys, Pair[], (0.0, 2.0), callback = test_callback)
cbs = get_callback(prob)
@test cbs isa CallbackSet
@test cbs.discrete_callbacks[1] == test_callback

prob = ODEProblem(sys2, Pair[], (0.0, 3.0))
cb = get_callback(prob)
@test cb isa ModelingToolkit.DiffEqCallbacks.VectorContinuousCallback

cond = cb.condition
out = [0.0, 0.0]
# the root to find is 2
cond.rf_ip(out, [0, 0], p0, t0)
@test out[1] ≈ -2 # signature is u,p,t
cond.rf_ip(out, [1, 0], p0, t0)
@test out[1] ≈ -1  # signature is u,p,t
cond.rf_ip(out, [2, 0], p0, t0) # this should return 0
@test out[1] ≈ 0  # signature is u,p,t

# the root to find is 1
out = [0.0, 0.0]
cond.rf_ip(out, [0, 0], p0, t0)
@test out[2] ≈ -1 # signature is u,p,t
cond.rf_ip(out, [0, 1], p0, t0) # this should return 0
@test out[2] ≈ 0  # signature is u,p,t
cond.rf_ip(out, [0, 2], p0, t0)
@test out[2] ≈ 1  # signature is u,p,t

sol = solve(prob, Tsit5())
@test minimum(t -> abs(t - 1), sol.t) < 1e-10 # test that the solver stepped at the first root
@test minimum(t -> abs(t - 2), sol.t) < 1e-10 # test that the solver stepped at the second root

@named sys = ODESystem(eqs, continuous_events = [x ~ 1, x ~ 2]) # two root eqs using the same state
prob = ODEProblem(sys, Pair[], (0.0, 3.0))
@test get_callback(prob) isa ModelingToolkit.DiffEqCallbacks.VectorContinuousCallback
sol = solve(prob, Tsit5())
@test minimum(t -> abs(t - 1), sol.t) < 1e-10 # test that the solver stepped at the first root
@test minimum(t -> abs(t - 2), sol.t) < 1e-10 # test that the solver stepped at the second root

## Test bouncing ball with equation affect
@variables t x(t)=1 v(t)=0
D = Differential(t)

root_eqs = [x ~ 0]
affect = [v ~ -v]

@named ball = ODESystem([D(x) ~ v
        D(v) ~ -9.8], t, continuous_events = root_eqs => affect)

@test getfield(ball, :continuous_events)[] ==
      SymbolicContinuousCallback(Equation[x ~ 0], Equation[v ~ -v])
ball = structural_simplify(ball)

@test length(ModelingToolkit.continuous_events(ball)) == 1

tspan = (0.0, 5.0)
prob = ODEProblem(ball, Pair[], tspan)
sol = solve(prob, Tsit5())
@test 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close
# plot(sol)

## Test bouncing ball in 2D with walls
@variables t x(t)=1 y(t)=0 vx(t)=0 vy(t)=1
D = Differential(t)

continuous_events = [[x ~ 0] => [vx ~ -vx]
    [y ~ -1.5, y ~ 1.5] => [vy ~ -vy]]

@named ball = ODESystem([D(x) ~ vx
        D(y) ~ vy
        D(vx) ~ -9.8
        D(vy) ~ -0.01vy], t; continuous_events)

ball = structural_simplify(ball)

tspan = (0.0, 5.0)
prob = ODEProblem(ball, Pair[], tspan)

cb = get_callback(prob)
@test cb isa ModelingToolkit.DiffEqCallbacks.VectorContinuousCallback
@test getfield(ball, :continuous_events)[1] ==
      SymbolicContinuousCallback(Equation[x ~ 0], Equation[vx ~ -vx])
@test getfield(ball, :continuous_events)[2] ==
      SymbolicContinuousCallback(Equation[y ~ -1.5, y ~ 1.5], Equation[vy ~ -vy])
cond = cb.condition
out = [0.0, 0.0, 0.0]
cond.rf_ip(out, [0, 0, 0, 0], p0, t0)
@test out ≈ [0, 1.5, -1.5]

sol = solve(prob, Tsit5())
@test 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close
@test minimum(sol[y]) ≈ -1.5 # check wall conditions
@test maximum(sol[y]) ≈ 1.5  # check wall conditions

# tv = sort([LinRange(0, 5, 200); sol.t])
# plot(sol(tv)[y], sol(tv)[x], line_z=tv)
# vline!([-1.5, 1.5], l=(:black, 5), primary=false)
# hline!([0], l=(:black, 5), primary=false)

## Test multi-variable affect
# in this test, there are two variables affected by a single event.
continuous_events = [
    [x ~ 0] => [vx ~ -vx, vy ~ -vy],
]

@named ball = ODESystem([D(x) ~ vx
        D(y) ~ vy
        D(vx) ~ -1
        D(vy) ~ 0], t; continuous_events)

ball = structural_simplify(ball)

tspan = (0.0, 5.0)
prob = ODEProblem(ball, Pair[], tspan)
sol = solve(prob, Tsit5())
@test 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close
@test -minimum(sol[y]) ≈ maximum(sol[y]) ≈ sqrt(2)  # the ball will never go further than √2 in either direction (gravity was changed to 1 to get this particular number)

# tv = sort([LinRange(0, 5, 200); sol.t])
# plot(sol(tv)[y], sol(tv)[x], line_z=tv)
# vline!([-1.5, 1.5], l=(:black, 5), primary=false)
# hline!([0], l=(:black, 5), primary=false)

# issue https://github.com/SciML/ModelingToolkit.jl/issues/1386
# tests that it works for ODAESystem
@variables vs(t) v(t) vmeasured(t)
eq = [vs ~ sin(2pi * t)
    D(v) ~ vs - v
    D(vmeasured) ~ 0.0]
ev = [sin(20pi * t) ~ 0.0] => [vmeasured ~ v]
@named sys = ODESystem(eq, continuous_events = ev)
sys = structural_simplify(sys)
prob = ODAEProblem(sys, zeros(2), (0.0, 5.1))
sol = solve(prob, Tsit5())
@test all(minimum((0:0.1:5) .- sol.t', dims = 2) .< 0.0001) # test that the solver stepped every 0.1s as dictated by event
@test sol([0.25])[vmeasured][] == sol([0.23])[vmeasured][] # test the hold property

##  https://github.com/SciML/ModelingToolkit.jl/issues/1528
Dₜ = Differential(t)

@parameters u(t) [input = true]  # Indicate that this is a controlled input
@parameters y(t) [output = true] # Indicate that this is a measured output

function Mass(; name, m = 1.0, p = 0, v = 0)
    ps = @parameters m = m
    sts = @variables pos(t)=p vel(t)=v
    eqs = Dₜ(pos) ~ vel
    ODESystem(eqs, t, [pos, vel], ps; name)
end
function Spring(; name, k = 1e4)
    ps = @parameters k = k
    @variables x(t) = 0 # Spring deflection
    ODESystem(Equation[], t, [x], ps; name)
end
function Damper(; name, c = 10)
    ps = @parameters c = c
    @variables vel(t) = 0
    ODESystem(Equation[], t, [vel], ps; name)
end
function SpringDamper(; name, k = false, c = false)
    spring = Spring(; name = :spring, k)
    damper = Damper(; name = :damper, c)
    compose(ODESystem(Equation[], t; name),
        spring, damper)
end
connect_sd(sd, m1, m2) = [sd.spring.x ~ m1.pos - m2.pos, sd.damper.vel ~ m1.vel - m2.vel]
sd_force(sd) = -sd.spring.k * sd.spring.x - sd.damper.c * sd.damper.vel
@named mass1 = Mass(; m = 1)
@named mass2 = Mass(; m = 1)
@named sd = SpringDamper(; k = 1000, c = 10)
function Model(u, d = 0)
    eqs = [connect_sd(sd, mass1, mass2)
        Dₜ(mass1.vel) ~ (sd_force(sd) + u) / mass1.m
        Dₜ(mass2.vel) ~ (-sd_force(sd) + d) / mass2.m]
    @named _model = ODESystem(eqs, t; observed = [y ~ mass2.pos])
    @named model = compose(_model, mass1, mass2, sd)
end
model = Model(sin(30t))
sys = structural_simplify(model)
@test isempty(ModelingToolkit.continuous_events(sys))

let
    function testsol(osys, u0, p, tspan; tstops = Float64[], skipparamtest = false,
        kwargs...)
        oprob = ODEProblem(osys, u0, tspan, p; kwargs...)
        sol = solve(oprob, Tsit5(); tstops = tstops, abstol = 1e-10, reltol = 1e-10)
        @test isapprox(sol(1.0000000001)[1] - sol(0.999999999)[1], 1.0; rtol = 1e-6)
        !skipparamtest && (@test oprob.p[1] == 1.0)
        @test isapprox(sol(4.0)[1], 2 * exp(-2.0))
        sol
    end

    @parameters k t1 t2
    @variables t A(t) B(t)

    cond1 = (t == t1)
    affect1 = [A ~ A + 1]
    cb1 = cond1 => affect1
    cond2 = (t == t2)
    affect2 = [k ~ 1.0]
    cb2 = cond2 => affect2

    ∂ₜ = Differential(t)
    eqs = [∂ₜ(A) ~ -k * A]
    @named osys = ODESystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2])
    u0 = [A => 1.0]
    p = [k => 0.0, t1 => 1.0, t2 => 2.0]
    tspan = (0.0, 4.0)
    testsol(osys, u0, p, tspan; tstops = [1.0, 2.0])

    cond1a = (t == t1)
    affect1a = [A ~ A + 1, B ~ A]
    cb1a = cond1a => affect1a
    @named osys1 = ODESystem(eqs, t, [A, B], [k, t1, t2], discrete_events = [cb1a, cb2])
    u0′ = [A => 1.0, B => 0.0]
    sol = testsol(osys1, u0′, p, tspan; tstops = [1.0, 2.0], check_length = false)
    @test sol(1.0000001, idxs = B) == 2.0

    # same as above - but with set-time event syntax
    cb1‵ = [1.0] => affect1 # needs to be a Vector for the event to happen only once
    cb2‵ = [2.0] => affect2
    @named osys‵ = ODESystem(eqs, t, [A], [k], discrete_events = [cb1‵, cb2‵])
    testsol(osys‵, u0, p, tspan)

    # mixing discrete affects
    @named osys3 = ODESystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵])
    testsol(osys3, u0, p, tspan; tstops = [1.0])

    # mixing with a func affect
    function affect!(integrator, u, p, ctx)
        integrator.p[p.k] = 1.0
        nothing
    end
    cb2‵‵ = [2.0] => (affect!, [], [k], nothing)
    @named osys4 = ODESystem(eqs, t, [A], [k, t1], discrete_events = [cb1, cb2‵‵])
    oprob4 = ODEProblem(osys4, u0, tspan, p)
    testsol(osys4, u0, p, tspan; tstops = [1.0])

    # mixing with symbolic condition in the func affect
    cb2‵‵‵ = (t == t2) => (affect!, [], [k], nothing)
    @named osys5 = ODESystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵‵‵])
    testsol(osys5, u0, p, tspan; tstops = [1.0, 2.0])
    @named osys6 = ODESystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb2‵‵‵, cb1])
    testsol(osys6, u0, p, tspan; tstops = [1.0, 2.0])

    # mix a continuous event too
    cond3 = A ~ 0.1
    affect3 = [k ~ 0.0]
    cb3 = cond3 => affect3
    @named osys7 = ODESystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵‵‵],
        continuous_events = [cb3])
    sol = testsol(osys7, u0, p, (0.0, 10.0); tstops = [1.0, 2.0], skipparamtest = true)
    @test isapprox(sol(10.0)[1], 0.1; atol = 1e-10, rtol = 1e-10)
end

let
    function testsol(ssys, u0, p, tspan; tstops = Float64[], skipparamtest = false,
        kwargs...)
        sprob = SDEProblem(ssys, u0, tspan, p; kwargs...)
        sol = solve(sprob, RI5(); tstops = tstops, abstol = 1e-10, reltol = 1e-10)
        @test isapprox(sol(1.0000000001)[1] - sol(0.999999999)[1], 1.0; rtol = 1e-4)
        !skipparamtest && (@test sprob.p[1] == 1.0)
        @test isapprox(sol(4.0)[1], 2 * exp(-2.0), atol = 1e-4)
        sol
    end

    @parameters k t1 t2
    @variables t A(t) B(t)

    cond1 = (t == t1)
    affect1 = [A ~ A + 1]
    cb1 = cond1 => affect1
    cond2 = (t == t2)
    affect2 = [k ~ 1.0]
    cb2 = cond2 => affect2

    ∂ₜ = Differential(t)
    eqs = [∂ₜ(A) ~ -k * A]
    @named ssys = SDESystem(eqs, Equation[], t, [A], [k, t1, t2],
        discrete_events = [cb1, cb2])
    u0 = [A => 1.0]
    p = [k => 0.0, t1 => 1.0, t2 => 2.0]
    tspan = (0.0, 4.0)
    testsol(ssys, u0, p, tspan; tstops = [1.0, 2.0])

    cond1a = (t == t1)
    affect1a = [A ~ A + 1, B ~ A]
    cb1a = cond1a => affect1a
    @named ssys1 = SDESystem(eqs, Equation[], t, [A, B], [k, t1, t2],
        discrete_events = [cb1a, cb2])
    u0′ = [A => 1.0, B => 0.0]
    sol = testsol(ssys1, u0′, p, tspan; tstops = [1.0, 2.0], check_length = false)
    @test sol(1.0000001, idxs = 2) == 2.0

    # same as above - but with set-time event syntax
    cb1‵ = [1.0] => affect1 # needs to be a Vector for the event to happen only once
    cb2‵ = [2.0] => affect2
    @named ssys‵ = SDESystem(eqs, Equation[], t, [A], [k], discrete_events = [cb1‵, cb2‵])
    testsol(ssys‵, u0, p, tspan)

    # mixing discrete affects
    @named ssys3 = SDESystem(eqs, Equation[], t, [A], [k, t1, t2],
        discrete_events = [cb1, cb2‵])
    testsol(ssys3, u0, p, tspan; tstops = [1.0])

    # mixing with a func affect
    function affect!(integrator, u, p, ctx)
        integrator.p[p.k] = 1.0
        nothing
    end
    cb2‵‵ = [2.0] => (affect!, [], [k], nothing)
    @named ssys4 = SDESystem(eqs, Equation[], t, [A], [k, t1],
        discrete_events = [cb1, cb2‵‵])
    testsol(ssys4, u0, p, tspan; tstops = [1.0])

    # mixing with symbolic condition in the func affect
    cb2‵‵‵ = (t == t2) => (affect!, [], [k], nothing)
    @named ssys5 = SDESystem(eqs, Equation[], t, [A], [k, t1, t2],
        discrete_events = [cb1, cb2‵‵‵])
    testsol(ssys5, u0, p, tspan; tstops = [1.0, 2.0])
    @named ssys6 = SDESystem(eqs, Equation[], t, [A], [k, t1, t2],
        discrete_events = [cb2‵‵‵, cb1])
    testsol(ssys6, u0, p, tspan; tstops = [1.0, 2.0])

    # mix a continuous event too
    cond3 = A ~ 0.1
    affect3 = [k ~ 0.0]
    cb3 = cond3 => affect3
    @named ssys7 = SDESystem(eqs, Equation[], t, [A], [k, t1, t2],
        discrete_events = [cb1, cb2‵‵‵],
        continuous_events = [cb3])
    sol = testsol(ssys7, u0, p, (0.0, 10.0); tstops = [1.0, 2.0], skipparamtest = true)
    @test isapprox(sol(10.0)[1], 0.1; atol = 1e-10, rtol = 1e-10)
end

let rng = rng
    function testsol(jsys, u0, p, tspan; tstops = Float64[], skipparamtest = false,
        N = 40000, kwargs...)
        dprob = DiscreteProblem(jsys, u0, tspan, p)
        jprob = JumpProblem(jsys, dprob, Direct(); kwargs...)
        sol = solve(jprob, SSAStepper(); tstops = tstops)
        @test (sol(1.000000000001)[1] - sol(0.99999999999)[1]) == 1
        !skipparamtest && (@test dprob.p[1] == 1.0)
        @test sol(40.0)[1] == 0
        sol
    end

    @parameters k t1 t2
    @variables t A(t) B(t)

    cond1 = (t == t1)
    affect1 = [A ~ A + 1]
    cb1 = cond1 => affect1
    cond2 = (t == t2)
    affect2 = [k ~ 1.0]
    cb2 = cond2 => affect2

    eqs = [MassActionJump(k, [A => 1], [A => -1])]
    @named jsys = JumpSystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2])
    u0 = [A => 1]
    p = [k => 0.0, t1 => 1.0, t2 => 2.0]
    tspan = (0.0, 40.0)
    testsol(jsys, u0, p, tspan; tstops = [1.0, 2.0], rng)

    cond1a = (t == t1)
    affect1a = [A ~ A + 1, B ~ A]
    cb1a = cond1a => affect1a
    @named jsys1 = JumpSystem(eqs, t, [A, B], [k, t1, t2], discrete_events = [cb1a, cb2])
    u0′ = [A => 1, B => 0]
    sol = testsol(jsys1, u0′, p, tspan; tstops = [1.0, 2.0], check_length = false, rng)
    @test sol(1.000000001, idxs = B) == 2

    # same as above - but with set-time event syntax
    cb1‵ = [1.0] => affect1 # needs to be a Vector for the event to happen only once
    cb2‵ = [2.0] => affect2
    @named jsys‵ = JumpSystem(eqs, t, [A], [k], discrete_events = [cb1‵, cb2‵])
    testsol(jsys‵, u0, [p[1]], tspan; rng)

    # mixing discrete affects
    @named jsys3 = JumpSystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵])
    testsol(jsys3, u0, p, tspan; tstops = [1.0], rng)

    # mixing with a func affect
    function affect!(integrator, u, p, ctx)
        integrator.p[p.k] = 1.0
        reset_aggregated_jumps!(integrator)
        nothing
    end
    cb2‵‵ = [2.0] => (affect!, [], [k], nothing)
    @named jsys4 = JumpSystem(eqs, t, [A], [k, t1], discrete_events = [cb1, cb2‵‵])
    testsol(jsys4, u0, p, tspan; tstops = [1.0], rng)

    # mixing with symbolic condition in the func affect
    cb2‵‵‵ = (t == t2) => (affect!, [], [k], nothing)
    @named jsys5 = JumpSystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵‵‵])
    testsol(jsys5, u0, p, tspan; tstops = [1.0, 2.0], rng)
    @named jsys6 = JumpSystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb2‵‵‵, cb1])
    testsol(jsys6, u0, p, tspan; tstops = [1.0, 2.0], rng)
end

let
    @variables t
    D = Differential(t)

    function oscillator_ce(k = 1.0; name)
        sts = @variables x(t)=1.0 v(t)=0.0 F(t)
        ps = @parameters k=k Θ=0.5
        eqs = [D(x) ~ v, D(v) ~ -k * x + F]
        ev = [x ~ Θ] => [x ~ 1.0, v ~ 0.0]
        ODESystem(eqs, t, sts, ps, continuous_events = [ev]; name)
    end

    @named oscce = oscillator_ce()
    eqs = [oscce.F ~ 0]
    @named eqs_sys = ODESystem(eqs, t)
    @named oneosc_ce = compose(eqs_sys, oscce)
    oneosc_ce_simpl = structural_simplify(oneosc_ce)

    prob = ODEProblem(oneosc_ce_simpl, [], (0.0, 2.0), [])
    sol = solve(prob, Tsit5(), saveat = 0.1)

    @test typeof(oneosc_ce_simpl) == ODESystem
    @test sol[1, 6] < 1.0 # test whether x(t) decreases over time
    @test sol[1, 18] > 0.5 # test whether event happened
end
