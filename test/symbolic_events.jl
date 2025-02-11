using ModelingToolkit, OrdinaryDiffEq, StochasticDiffEq, JumpProcesses, Test
using SciMLStructures: canonicalize, Discrete
using ModelingToolkit: SymbolicContinuousCallback,
                       SymbolicContinuousCallbacks, NULL_AFFECT,
                       get_callback,
                       t_nounits as t,
                       D_nounits as D
using StableRNGs
import SciMLBase
using SymbolicIndexingInterface
using Setfield
rng = StableRNG(12345)

@variables x(t) = 0

eqs = [D(x) ~ 1]
affect = [x ~ 0]
affect_neg = [x ~ 1]

## Test SymbolicContinuousCallback
@testset "SymbolicContinuousCallback constructors" begin
    e = SymbolicContinuousCallback(eqs[])
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT
    @test e.affect_neg == NULL_AFFECT
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT
    @test e.affect_neg == NULL_AFFECT
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs, NULL_AFFECT)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT
    @test e.affect_neg == NULL_AFFECT
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs[], NULL_AFFECT)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT
    @test e.affect_neg == NULL_AFFECT
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs => NULL_AFFECT)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT
    @test e.affect_neg == NULL_AFFECT
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs[] => NULL_AFFECT)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT
    @test e.affect_neg == NULL_AFFECT
    @test e.rootfind == SciMLBase.LeftRootFind

    ## With affect

    e = SymbolicContinuousCallback(eqs[], affect)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test e.affect_neg == affect
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs, affect)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test e.affect_neg == affect
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs, affect)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test e.affect_neg == affect
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs[], affect)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test e.affect_neg == affect
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs => affect)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test e.affect_neg == affect
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs[] => affect)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test e.affect_neg == affect
    @test e.rootfind == SciMLBase.LeftRootFind

    # with only positive edge affect

    e = SymbolicContinuousCallback(eqs[], affect, affect_neg = nothing)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test isnothing(e.affect_neg)
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs, affect, affect_neg = nothing)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test isnothing(e.affect_neg)
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs, affect, affect_neg = nothing)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test isnothing(e.affect_neg)
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs[], affect, affect_neg = nothing)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test isnothing(e.affect_neg)
    @test e.rootfind == SciMLBase.LeftRootFind

    # with explicit edge affects

    e = SymbolicContinuousCallback(eqs[], affect, affect_neg = affect_neg)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test e.affect_neg == affect_neg
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs, affect, affect_neg = affect_neg)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test e.affect_neg == affect_neg
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs, affect, affect_neg = affect_neg)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test e.affect_neg == affect_neg
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(eqs[], affect, affect_neg = affect_neg)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test e.affect_neg == affect_neg
    @test e.rootfind == SciMLBase.LeftRootFind

    # with different root finding ops

    e = SymbolicContinuousCallback(
        eqs[], affect, affect_neg = affect_neg, rootfind = SciMLBase.LeftRootFind)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test e.affect_neg == affect_neg
    @test e.rootfind == SciMLBase.LeftRootFind

    e = SymbolicContinuousCallback(
        eqs[], affect, affect_neg = affect_neg, rootfind = SciMLBase.RightRootFind)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test e.affect_neg == affect_neg
    @test e.rootfind == SciMLBase.RightRootFind

    e = SymbolicContinuousCallback(
        eqs[], affect, affect_neg = affect_neg, rootfind = SciMLBase.NoRootFind)
    @test e isa SymbolicContinuousCallback
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
    @test e.affect_neg == affect_neg
    @test e.rootfind == SciMLBase.NoRootFind
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

##

@named sys = ODESystem(eqs, t, continuous_events = [x ~ 1])
@test getfield(sys, :continuous_events)[] ==
      SymbolicContinuousCallback(Equation[x ~ 1], NULL_AFFECT)
@test isequal(equations(getfield(sys, :continuous_events))[], x ~ 1)
fsys = flatten(sys)
@test isequal(equations(getfield(fsys, :continuous_events))[], x ~ 1)

@named sys2 = ODESystem([D(x) ~ 1], t, continuous_events = [x ~ 2], systems = [sys])
@test getfield(sys2, :continuous_events)[] ==
      SymbolicContinuousCallback(Equation[x ~ 2], NULL_AFFECT)
@test all(ModelingToolkit.continuous_events(sys2) .== [
    SymbolicContinuousCallback(Equation[x ~ 2], NULL_AFFECT),
    SymbolicContinuousCallback(Equation[sys.x ~ 1], NULL_AFFECT)
])

@test isequal(equations(getfield(sys2, :continuous_events))[1], x ~ 2)
@test length(ModelingToolkit.continuous_events(sys2)) == 2
@test isequal(ModelingToolkit.continuous_events(sys2)[1].eqs[], x ~ 2)
@test isequal(ModelingToolkit.continuous_events(sys2)[2].eqs[], sys.x ~ 1)

sys = complete(sys)
sys_nosplit = complete(sys; split = false)
sys2 = complete(sys2)
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
prob_nosplit = ODEProblem(sys_nosplit, Pair[], (0.0, 2.0))
sol = solve(prob, Tsit5())
sol_nosplit = solve(prob_nosplit, Tsit5())
@test minimum(t -> abs(t - 1), sol.t) < 1e-10 # test that the solver stepped at the root
@test minimum(t -> abs(t - 1), sol_nosplit.t) < 1e-10 # test that the solver stepped at the root

# Test that a user provided callback is respected
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

@named sys = ODESystem(eqs, t, continuous_events = [x ~ 1, x ~ 2]) # two root eqs using the same unknown
sys = complete(sys)
prob = ODEProblem(sys, Pair[], (0.0, 3.0))
@test get_callback(prob) isa ModelingToolkit.DiffEqCallbacks.VectorContinuousCallback
sol = solve(prob, Tsit5())
@test minimum(t -> abs(t - 1), sol.t) < 1e-10 # test that the solver stepped at the first root
@test minimum(t -> abs(t - 2), sol.t) < 1e-10 # test that the solver stepped at the second root

## Test bouncing ball with equation affect
@variables x(t)=1 v(t)=0

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
@variables x(t)=1 y(t)=0 vx(t)=0 vy(t)=1

continuous_events = [[x ~ 0] => [vx ~ -vx]
                     [y ~ -1.5, y ~ 1.5] => [vy ~ -vy]]

@named ball = ODESystem(
    [D(x) ~ vx
     D(y) ~ vy
     D(vx) ~ -9.8
     D(vy) ~ -0.01vy], t; continuous_events)

_ball = ball
ball = structural_simplify(_ball)
ball_nosplit = structural_simplify(_ball; split = false)

tspan = (0.0, 5.0)
prob = ODEProblem(ball, Pair[], tspan)
prob_nosplit = ODEProblem(ball_nosplit, Pair[], tspan)

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
sol_nosplit = solve(prob_nosplit, Tsit5())
@test 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close
@test minimum(sol[y]) ≈ -1.5 # check wall conditions
@test maximum(sol[y]) ≈ 1.5  # check wall conditions
@test 0 <= minimum(sol_nosplit[x]) <= 1e-10 # the ball never went through the floor but got very close
@test minimum(sol_nosplit[y]) ≈ -1.5 # check wall conditions
@test maximum(sol_nosplit[y]) ≈ 1.5  # check wall conditions

# tv = sort([LinRange(0, 5, 200); sol.t])
# plot(sol(tv)[y], sol(tv)[x], line_z=tv)
# vline!([-1.5, 1.5], l=(:black, 5), primary=false)
# hline!([0], l=(:black, 5), primary=false)

## Test multi-variable affect
# in this test, there are two variables affected by a single event.
continuous_events = [
    [x ~ 0] => [vx ~ -vx, vy ~ -vy]
]

@named ball = ODESystem([D(x) ~ vx
                         D(y) ~ vy
                         D(vx) ~ -1
                         D(vy) ~ 0], t; continuous_events)

ball_nosplit = structural_simplify(ball)
ball = structural_simplify(ball)

tspan = (0.0, 5.0)
prob = ODEProblem(ball, Pair[], tspan)
prob_nosplit = ODEProblem(ball_nosplit, Pair[], tspan)
sol = solve(prob, Tsit5())
sol_nosplit = solve(prob_nosplit, Tsit5())
@test 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close
@test -minimum(sol[y]) ≈ maximum(sol[y]) ≈ sqrt(2)  # the ball will never go further than √2 in either direction (gravity was changed to 1 to get this particular number)
@test 0 <= minimum(sol_nosplit[x]) <= 1e-10 # the ball never went through the floor but got very close
@test -minimum(sol_nosplit[y]) ≈ maximum(sol_nosplit[y]) ≈ sqrt(2)  # the ball will never go further than √2 in either direction (gravity was changed to 1 to get this particular number)

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
@named sys = ODESystem(eq, t, continuous_events = ev)
sys = structural_simplify(sys)
prob = ODEProblem(sys, zeros(2), (0.0, 5.1))
sol = solve(prob, Tsit5())
@test all(minimum((0:0.1:5) .- sol.t', dims = 2) .< 0.0001) # test that the solver stepped every 0.1s as dictated by event
@test sol([0.25])[vmeasured][] == sol([0.23])[vmeasured][] # test the hold property

##  https://github.com/SciML/ModelingToolkit.jl/issues/1528
Dₜ = D

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
    function testsol(osys, u0, p, tspan; tstops = Float64[], paramtotest = nothing,
            kwargs...)
        oprob = ODEProblem(complete(osys), u0, tspan, p; kwargs...)
        sol = solve(oprob, Tsit5(); tstops = tstops, abstol = 1e-10, reltol = 1e-10)
        @test isapprox(sol(1.0000000001)[1] - sol(0.999999999)[1], 1.0; rtol = 1e-6)
        paramtotest === nothing || (@test sol.ps[paramtotest] == 1.0)
        @test isapprox(sol(4.0)[1], 2 * exp(-2.0))
        sol
    end

    @parameters k t1 t2
    @variables A(t) B(t)

    cond1 = (t == t1)
    affect1 = [A ~ A + 1]
    cb1 = cond1 => affect1
    cond2 = (t == t2)
    affect2 = [k ~ 1.0]
    cb2 = cond2 => affect2

    ∂ₜ = D
    eqs = [∂ₜ(A) ~ -k * A]
    @named osys = ODESystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2])
    u0 = [A => 1.0]
    p = [k => 0.0, t1 => 1.0, t2 => 2.0]
    tspan = (0.0, 4.0)
    testsol(osys, u0, p, tspan; tstops = [1.0, 2.0], paramtotest = k)

    cond1a = (t == t1)
    affect1a = [A ~ A + 1, B ~ A]
    cb1a = cond1a => affect1a
    @named osys1 = ODESystem(eqs, t, [A, B], [k, t1, t2], discrete_events = [cb1a, cb2])
    u0′ = [A => 1.0, B => 0.0]
    sol = testsol(
        osys1, u0′, p, tspan; tstops = [1.0, 2.0], check_length = false, paramtotest = k)
    @test sol(1.0000001, idxs = B) == 2.0

    # same as above - but with set-time event syntax
    cb1‵ = [1.0] => affect1 # needs to be a Vector for the event to happen only once
    cb2‵ = [2.0] => affect2
    @named osys‵ = ODESystem(eqs, t, [A], [k], discrete_events = [cb1‵, cb2‵])
    testsol(osys‵, u0, p, tspan; paramtotest = k)

    # mixing discrete affects
    @named osys3 = ODESystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵])
    testsol(osys3, u0, p, tspan; tstops = [1.0], paramtotest = k)

    # mixing with a func affect
    function affect!(integrator, u, p, ctx)
        integrator.ps[p.k] = 1.0
        nothing
    end
    cb2‵‵ = [2.0] => (affect!, [], [k], [k], nothing)
    @named osys4 = ODESystem(eqs, t, [A], [k, t1], discrete_events = [cb1, cb2‵‵])
    oprob4 = ODEProblem(complete(osys4), u0, tspan, p)
    testsol(osys4, u0, p, tspan; tstops = [1.0], paramtotest = k)

    # mixing with symbolic condition in the func affect
    cb2‵‵‵ = (t == t2) => (affect!, [], [k], [k], nothing)
    @named osys5 = ODESystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵‵‵])
    testsol(osys5, u0, p, tspan; tstops = [1.0, 2.0])
    @named osys6 = ODESystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb2‵‵‵, cb1])
    testsol(osys6, u0, p, tspan; tstops = [1.0, 2.0], paramtotest = k)

    # mix a continuous event too
    cond3 = A ~ 0.1
    affect3 = [k ~ 0.0]
    cb3 = cond3 => affect3
    @named osys7 = ODESystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵‵‵],
        continuous_events = [cb3])
    sol = testsol(osys7, u0, p, (0.0, 10.0); tstops = [1.0, 2.0])
    @test isapprox(sol(10.0)[1], 0.1; atol = 1e-10, rtol = 1e-10)
end

let
    function testsol(ssys, u0, p, tspan; tstops = Float64[], paramtotest = nothing,
            kwargs...)
        sprob = SDEProblem(complete(ssys), u0, tspan, p; kwargs...)
        sol = solve(sprob, RI5(); tstops = tstops, abstol = 1e-10, reltol = 1e-10)
        @test isapprox(sol(1.0000000001)[1] - sol(0.999999999)[1], 1.0; rtol = 1e-4)
        paramtotest === nothing || (@test sol.ps[paramtotest] == 1.0)
        @test isapprox(sol(4.0)[1], 2 * exp(-2.0), atol = 1e-4)
        sol
    end

    @parameters k t1 t2
    @variables A(t) B(t)

    cond1 = (t == t1)
    affect1 = [A ~ A + 1]
    cb1 = cond1 => affect1
    cond2 = (t == t2)
    affect2 = [k ~ 1.0]
    cb2 = cond2 => affect2

    ∂ₜ = D
    eqs = [∂ₜ(A) ~ -k * A]
    @named ssys = SDESystem(eqs, [0.0], t, [A], [k, t1, t2],
        discrete_events = [cb1, cb2])
    u0 = [A => 1.0]
    p = [k => 0.0, t1 => 1.0, t2 => 2.0]
    tspan = (0.0, 4.0)
    testsol(ssys, u0, p, tspan; tstops = [1.0, 2.0], paramtotest = k)

    cond1a = (t == t1)
    affect1a = [A ~ A + 1, B ~ A]
    cb1a = cond1a => affect1a
    @named ssys1 = SDESystem(eqs, [0.0], t, [A, B], [k, t1, t2],
        discrete_events = [cb1a, cb2])
    u0′ = [A => 1.0, B => 0.0]
    sol = testsol(
        ssys1, u0′, p, tspan; tstops = [1.0, 2.0], check_length = false, paramtotest = k)
    @test sol(1.0000001, idxs = 2) == 2.0

    # same as above - but with set-time event syntax
    cb1‵ = [1.0] => affect1 # needs to be a Vector for the event to happen only once
    cb2‵ = [2.0] => affect2
    @named ssys‵ = SDESystem(eqs, [0.0], t, [A], [k], discrete_events = [cb1‵, cb2‵])
    testsol(ssys‵, u0, p, tspan; paramtotest = k)

    # mixing discrete affects
    @named ssys3 = SDESystem(eqs, [0.0], t, [A], [k, t1, t2],
        discrete_events = [cb1, cb2‵])
    testsol(ssys3, u0, p, tspan; tstops = [1.0], paramtotest = k)

    # mixing with a func affect
    function affect!(integrator, u, p, ctx)
        setp(integrator, p.k)(integrator, 1.0)
        nothing
    end
    cb2‵‵ = [2.0] => (affect!, [], [k], [k], nothing)
    @named ssys4 = SDESystem(eqs, [0.0], t, [A], [k, t1],
        discrete_events = [cb1, cb2‵‵])
    testsol(ssys4, u0, p, tspan; tstops = [1.0], paramtotest = k)

    # mixing with symbolic condition in the func affect
    cb2‵‵‵ = (t == t2) => (affect!, [], [k], [k], nothing)
    @named ssys5 = SDESystem(eqs, [0.0], t, [A], [k, t1, t2],
        discrete_events = [cb1, cb2‵‵‵])
    testsol(ssys5, u0, p, tspan; tstops = [1.0, 2.0], paramtotest = k)
    @named ssys6 = SDESystem(eqs, [0.0], t, [A], [k, t1, t2],
        discrete_events = [cb2‵‵‵, cb1])
    testsol(ssys6, u0, p, tspan; tstops = [1.0, 2.0], paramtotest = k)

    # mix a continuous event too
    cond3 = A ~ 0.1
    affect3 = [k ~ 0.0]
    cb3 = cond3 => affect3
    @named ssys7 = SDESystem(eqs, [0.0], t, [A], [k, t1, t2],
        discrete_events = [cb1, cb2‵‵‵],
        continuous_events = [cb3])
    sol = testsol(ssys7, u0, p, (0.0, 10.0); tstops = [1.0, 2.0])
    @test isapprox(sol(10.0)[1], 0.1; atol = 1e-10, rtol = 1e-10)
end

let rng = rng
    function testsol(jsys, u0, p, tspan; tstops = Float64[], paramtotest = nothing,
            N = 40000, kwargs...)
        jsys = complete(jsys)
        dprob = DiscreteProblem(jsys, u0, tspan, p)
        jprob = JumpProblem(jsys, dprob, Direct(); kwargs...)
        sol = solve(jprob, SSAStepper(); tstops = tstops)
        @test (sol(1.000000000001)[1] - sol(0.99999999999)[1]) == 1
        paramtotest === nothing || (@test sol.ps[paramtotest] == 1.0)
        @test sol(40.0)[1] == 0
        sol
    end

    @parameters k t1 t2
    @variables A(t) B(t)

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
    testsol(jsys, u0, p, tspan; tstops = [1.0, 2.0], rng, paramtotest = k)

    cond1a = (t == t1)
    affect1a = [A ~ A + 1, B ~ A]
    cb1a = cond1a => affect1a
    @named jsys1 = JumpSystem(eqs, t, [A, B], [k, t1, t2], discrete_events = [cb1a, cb2])
    u0′ = [A => 1, B => 0]
    sol = testsol(jsys1, u0′, p, tspan; tstops = [1.0, 2.0],
        check_length = false, rng, paramtotest = k)
    @test sol(1.000000001, idxs = B) == 2

    # same as above - but with set-time event syntax
    cb1‵ = [1.0] => affect1 # needs to be a Vector for the event to happen only once
    cb2‵ = [2.0] => affect2
    @named jsys‵ = JumpSystem(eqs, t, [A], [k], discrete_events = [cb1‵, cb2‵])
    testsol(jsys‵, u0, [p[1]], tspan; rng, paramtotest = k)

    # mixing discrete affects
    @named jsys3 = JumpSystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵])
    testsol(jsys3, u0, p, tspan; tstops = [1.0], rng, paramtotest = k)

    # mixing with a func affect
    function affect!(integrator, u, p, ctx)
        integrator.ps[p.k] = 1.0
        reset_aggregated_jumps!(integrator)
        nothing
    end
    cb2‵‵ = [2.0] => (affect!, [], [k], [k], nothing)
    @named jsys4 = JumpSystem(eqs, t, [A], [k, t1], discrete_events = [cb1, cb2‵‵])
    testsol(jsys4, u0, p, tspan; tstops = [1.0], rng, paramtotest = k)

    # mixing with symbolic condition in the func affect
    cb2‵‵‵ = (t == t2) => (affect!, [], [k], [k], nothing)
    @named jsys5 = JumpSystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb1, cb2‵‵‵])
    testsol(jsys5, u0, p, tspan; tstops = [1.0, 2.0], rng, paramtotest = k)
    @named jsys6 = JumpSystem(eqs, t, [A], [k, t1, t2], discrete_events = [cb2‵‵‵, cb1])
    testsol(jsys6, u0, p, tspan; tstops = [1.0, 2.0], rng, paramtotest = k)
end

let
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

@testset "Additional SymbolicContinuousCallback options" begin
    # baseline affect (pos + neg + left root find)
    @variables c1(t)=1.0 c2(t)=1.0 # c1 = cos(t), c2 = cos(3t)
    eqs = [D(c1) ~ -sin(t); D(c2) ~ -3 * sin(3 * t)]
    record_crossings(i, u, _, c) = push!(c, i.t => i.u[u.v])
    cr1 = []
    cr2 = []
    evt1 = ModelingToolkit.SymbolicContinuousCallback(
        [c1 ~ 0], (record_crossings, [c1 => :v], [], [], cr1))
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (record_crossings, [c2 => :v], [], [], cr2))
    @named trigsys = ODESystem(eqs, t; continuous_events = [evt1, evt2])
    trigsys_ss = structural_simplify(trigsys)
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
        [c1 ~ 0], (record_crossings, [c1 => :v], [], [], cr1p);
        affect_neg = (record_crossings, [c1 => :v], [], [], cr1n))
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (record_crossings, [c2 => :v], [], [], cr2p);
        affect_neg = (record_crossings, [c2 => :v], [], [], cr2n))
    @named trigsys = ODESystem(eqs, t; continuous_events = [evt1, evt2])
    trigsys_ss = structural_simplify(trigsys)
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
        [c1 ~ 0], (record_crossings, [c1 => :v], [], [], cr1p); affect_neg = nothing)
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (record_crossings, [c2 => :v], [], [], cr2p); affect_neg = nothing)
    @named trigsys = ODESystem(eqs, t; continuous_events = [evt1, evt2])
    trigsys_ss = structural_simplify(trigsys)
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
        [c1 ~ 0], (record_crossings, [c1 => :v], [], [], cr1p); affect_neg = nothing)
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (record_crossings, [c2 => :v], [], [], cr2p);
        affect_neg = (record_crossings, [c2 => :v], [], [], cr2n))
    @named trigsys = ODESystem(eqs, t; continuous_events = [evt1, evt2])
    trigsys_ss = structural_simplify(trigsys)
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
        [c1 ~ 0], (record_crossings, [c1 => :v], [], [], cr1);
        rootfind = SciMLBase.RightRootFind)
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (record_crossings, [c2 => :v], [], [], cr2);
        rootfind = SciMLBase.RightRootFind)
    @named trigsys = ODESystem(eqs, t; continuous_events = [evt1, evt2])
    trigsys_ss = structural_simplify(trigsys)
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
        [c1 ~ 0], (record_crossings, [c1 => :v], [], [], cr1);
        rootfind = SciMLBase.LeftRootFind)
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (record_crossings, [c2 => :v], [], [], cr2);
        rootfind = SciMLBase.RightRootFind)
    @named trigsys = ODESystem(eqs, t; continuous_events = [evt1, evt2])
    trigsys_ss = structural_simplify(trigsys)
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
        [c1 ~ 0], (record_crossings, [c1 => :v], [], [], cr1);
        rootfind = SciMLBase.LeftRootFind)
    evt2 = ModelingToolkit.SymbolicContinuousCallback(
        [c2 ~ 0], (record_crossings, [c2 => :v], [], [], cr2);
        rootfind = SciMLBase.RightRootFind)
    @named trigsys = ODESystem(eqs, t; continuous_events = [evt2, evt1])
    trigsys_ss = structural_simplify(trigsys)
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
            [60] => [
                binary_valve_1.S ~ 1.0, binary_valve_2.S ~ 0.0, binary_valve_2.Δp ~ 1.0]
            [120] => [binary_valve_1.S ~ 0.0, binary_valve_2.Δp ~ 0.0]
        end
    end

    # Test Simulation
    @mtkbuild sys = TestSystem()

    # Test Simulation
    prob = ODEProblem(sys, [], (0.0, 150.0))
    sol = solve(prob)
    @test sol[end] == [0.0, 0.0, 0.0]
end

@testset "Discrete variable timeseries" begin
    @variables x(t)
    @parameters a(t) b(t) c(t)
    cb1 = [x ~ 1.0] => [a ~ -a]
    function save_affect!(integ, u, p, ctx)
        integ.ps[p.b] = 5.0
    end
    cb2 = [x ~ 0.5] => (save_affect!, [], [b], [b], nothing)
    cb3 = 1.0 => [c ~ t]

    @mtkbuild sys = ODESystem(D(x) ~ cos(t), t, [x], [a, b, c];
        continuous_events = [cb1, cb2], discrete_events = [cb3])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 2pi), [a => 1.0, b => 2.0, c => 0.0])
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
    @named sys = ODESystem(
        eqs, t, [temp], params; continuous_events = [furnace_off, furnace_enable])
    ss = structural_simplify(sys)
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
    @named sys = ODESystem(
        eqs, t, [temp], params; continuous_events = [furnace_off, furnace_enable])
    ss = structural_simplify(sys)
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
    @named sys = ODESystem(eqs, t, [temp], params; continuous_events = [furnace_off])
    ss = structural_simplify(sys)
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
    @named sys = ODESystem(
        eqs, t, [temp, tempsq], params; continuous_events = [furnace_off])
    ss = structural_simplify(sys)
    @test_throws "refers to missing variable(s)" prob=ODEProblem(
        ss, [temp => 0.0, furnace_on => true], (0.0, 100.0))

    @parameters not_actually_here
    furnace_off = ModelingToolkit.SymbolicContinuousCallback(
        [temp ~ furnace_off_threshold],
        ModelingToolkit.ImperativeAffect(modified = (; furnace_on),
            observed = (; furnace_on, not_actually_here)) do x, o, c, i
            @set! x.furnace_on = false
        end)
    @named sys = ODESystem(
        eqs, t, [temp, tempsq], params; continuous_events = [furnace_off])
    ss = structural_simplify(sys)
    @test_throws "refers to missing variable(s)" prob=ODEProblem(
        ss, [temp => 0.0, furnace_on => true], (0.0, 100.0))

    furnace_off = ModelingToolkit.SymbolicContinuousCallback(
        [temp ~ furnace_off_threshold],
        ModelingToolkit.ImperativeAffect(modified = (; furnace_on),
            observed = (; furnace_on)) do x, o, c, i
            return (; fictional2 = false)
        end)
    @named sys = ODESystem(
        eqs, t, [temp, tempsq], params; continuous_events = [furnace_off])
    ss = structural_simplify(sys)
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
    @named sys = ODESystem(
        eqs, t, [theta, omega], params; continuous_events = [qAevt, qBevt])
    ss = structural_simplify(sys)
    prob = ODEProblem(ss, [theta => 1e-5], (0.0, pi))
    sol = solve(prob, Tsit5(); dtmax = 0.01)
    @test getp(sol, cnt)(sol) == 198 # we get 2 pulses per phase cycle (cos 0 crossing) and we go to 100 cycles; we miss a few due to the initial state
end

@testset "Initialization" begin
    @variables x(t)
    seen = false
    f = ModelingToolkit.FunctionalAffect(
        f = (i, u, p, c) -> seen = true, sts = [], pars = [], discretes = [])
    cb1 = ModelingToolkit.SymbolicContinuousCallback(
        [x ~ 0], Equation[], initialize = [x ~ 1.5], finalize = f)
    @mtkbuild sys = ODESystem(D(x) ~ -1, t, [x], []; continuous_events = [cb1])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 2), [])
    sol = solve(prob, Tsit5(); dtmax = 0.01)
    @test sol[x][1] ≈ 1.0
    @test sol[x][2] ≈ 1.5 # the initialize affect has been applied
    @test seen == true

    @variables x(t)
    seen = false
    f = ModelingToolkit.FunctionalAffect(
        f = (i, u, p, c) -> seen = true, sts = [], pars = [], discretes = [])
    cb1 = ModelingToolkit.SymbolicContinuousCallback(
        [x ~ 0], Equation[], initialize = [x ~ 1.5], finalize = f)
    inited = false
    finaled = false
    a = ModelingToolkit.FunctionalAffect(
        f = (i, u, p, c) -> inited = true, sts = [], pars = [], discretes = [])
    b = ModelingToolkit.FunctionalAffect(
        f = (i, u, p, c) -> finaled = true, sts = [], pars = [], discretes = [])
    cb2 = ModelingToolkit.SymbolicContinuousCallback(
        [x ~ 0.1], Equation[], initialize = a, finalize = b)
    @mtkbuild sys = ODESystem(D(x) ~ -1, t, [x], []; continuous_events = [cb1, cb2])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 2), [])
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
    @mtkbuild sys = ODESystem(D(x) ~ -1, t, [x], []; discrete_events = [cb3])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 2), [])
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
    @mtkbuild sys = ODESystem(D(x) ~ -1, t, [x], []; discrete_events = [cb3])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 2), [])
    sol = solve(prob, Tsit5())
    @test seen == true
    @test inited == true

    #preset
    seen = false
    inited = false
    finaled = false
    cb3 = ModelingToolkit.SymbolicDiscreteCallback([1.0], f, initialize = a, finalize = b)
    @mtkbuild sys = ODESystem(D(x) ~ -1, t, [x], []; discrete_events = [cb3])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 2), [])
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
    @mtkbuild sys = ODESystem(D(x) ~ -1, t, [x], []; discrete_events = [cb3])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 2), [])
    sol = solve(prob, Tsit5(); tstops = 1.0)
    @test seen == true
    @test inited == true
    @test finaled == true
end

@testset "Bump" begin
    @variables x(t) [irreducible = true] y(t) [irreducible = true]
    eqs = [x ~ y, D(x) ~ -1]
    cb = [x ~ 0.0] => [x ~ 0, y ~ 1]
    @mtkbuild pend = ODESystem(eqs, t; continuous_events = [cb])
    prob = ODEProblem(pend, [x => 1], (0.0, 3.0), guesses = [y => x])
    @test_throws "DAE initialization failed" solve(prob, Rodas5())

    cb = [x ~ 0.0] => [y ~ 1]
    @mtkbuild pend = ODESystem(eqs, t; continuous_events = [cb])
    prob = ODEProblem(pend, [x => 1], (0.0, 3.0), guesses = [y => x])
    @test_broken !SciMLBase.successful_retcode(solve(prob, Rodas5()))

    cb = [x ~ 0.0] => [x ~ 1, y ~ 1]
    @mtkbuild pend = ODESystem(eqs, t; continuous_events = [cb])
    prob = ODEProblem(pend, [x => 1], (0.0, 3.0), guesses = [y => x])
    @test all(≈(0.0; atol = 1e-9), solve(prob, Rodas5())[[x, y]][end])
end

@testset "Issue#3154 Array variable in discrete condition" begin
    @mtkmodel DECAY begin
        @parameters begin
            unrelated[1:2] = zeros(2)
            k = 0.0
        end
        @variables begin
            x(t) = 10.0
        end
        @equations begin
            D(x) ~ -k * x
        end
        @discrete_events begin
            (t == 1.0) => [k ~ 1.0]
        end
    end
    @mtkbuild decay = DECAY()
    prob = ODEProblem(decay, [], (0.0, 10.0), [])
    @test_nowarn solve(prob, Tsit5(), tstops = [1.0])
end

@testset "Array parameter updates in ImperativeEffect" begin
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
        return ODESystem(eqs, t, vars, params; name = name,
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
        return ODESystem(eqs, t, vars, params; name = name) # note no event
    end

    @named wd1 = weird1(0.021)
    @named wd2 = weird2(0.021)

    sys1 = structural_simplify(ODESystem([], t; name = :parent,
        discrete_events = [0.01 => ModelingToolkit.ImperativeAffect(
            modified = (; θs = reduce(vcat, [[wd1.θ]])), ctx = [1]) do m, o, c, i
            @set! m.θs[1] = c[] += 1
        end],
        systems = [wd1]))
    sys2 = structural_simplify(ODESystem([], t; name = :parent,
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
