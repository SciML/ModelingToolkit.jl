using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: SymbolicContinuousCallback, SymbolicContinuousCallbacks, NULL_AFFECT, get_callback

@parameters t
@variables x(t)=0 
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
@test getfield(sys, :continuous_events)[] == SymbolicContinuousCallback(Equation[x ~ 1], NULL_AFFECT)
@test isequal(equations(getfield(sys, :continuous_events))[], x ~ 1)
fsys = flatten(sys)
@test isequal(equations(getfield(fsys, :continuous_events))[], x ~ 1)

@named sys2 = ODESystem([D(x) ~ 1], continuous_events = [x ~ 2], systems=[sys])
@test getfield(sys2, :continuous_events)[] == SymbolicContinuousCallback(Equation[x ~ 2], NULL_AFFECT)
@test all(ModelingToolkit.continuous_events(sys2) .== [SymbolicContinuousCallback(Equation[x ~ 2], NULL_AFFECT), SymbolicContinuousCallback(Equation[sys.x ~ 1], NULL_AFFECT)])

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
@test minimum(t->abs(t-1), sol.t) < 1e-10 # test that the solver stepped at the root


# Test that a user provided callback is respected
test_callback = DiscreteCallback(x->x, x->x)
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
cond.rf_ip(out, [0,0], p0, t0)
@test out[1] ≈ -2 # signature is u,p,t
cond.rf_ip(out, [1,0], p0, t0)
@test out[1] ≈ -1  # signature is u,p,t
cond.rf_ip(out, [2,0], p0, t0) # this should return 0
@test out[1] ≈ 0  # signature is u,p,t

# the root to find is 1
out = [0.0, 0.0]
cond.rf_ip(out, [0,0], p0, t0)
@test out[2] ≈ -1 # signature is u,p,t
cond.rf_ip(out, [0,1], p0, t0) # this should return 0
@test out[2] ≈ 0  # signature is u,p,t
cond.rf_ip(out, [0,2], p0, t0)
@test out[2] ≈ 1  # signature is u,p,t

sol = solve(prob, Tsit5())
@test minimum(t->abs(t-1), sol.t) < 1e-10 # test that the solver stepped at the first root
@test minimum(t->abs(t-2), sol.t) < 1e-10 # test that the solver stepped at the second root

@named sys = ODESystem(eqs, continuous_events = [x ~ 1, x ~ 2]) # two root eqs using the same state
prob = ODEProblem(sys, Pair[], (0.0, 3.0))
@test get_callback(prob) isa ModelingToolkit.DiffEqCallbacks.VectorContinuousCallback
sol = solve(prob, Tsit5())
@test minimum(t->abs(t-1), sol.t) < 1e-10 # test that the solver stepped at the first root
@test minimum(t->abs(t-2), sol.t) < 1e-10 # test that the solver stepped at the second root


## Test bouncing ball with equation affect
@variables t x(t)=1 v(t)=0
D = Differential(t)

root_eqs = [x ~ 0]
affect   = [v ~ -v]

@named ball = ODESystem([
    D(x) ~ v
    D(v) ~ -9.8
], t, continuous_events = root_eqs => affect)

@test getfield(ball, :continuous_events)[] == SymbolicContinuousCallback(Equation[x ~ 0], Equation[v ~ -v])
ball = structural_simplify(ball)

tspan = (0.0,5.0)
prob = ODEProblem(ball, Pair[], tspan)
sol = solve(prob,Tsit5())
@test 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close
# plot(sol)


## Test bouncing ball in 2D with walls
@variables t x(t)=1 y(t)=0 vx(t)=0 vy(t)=1
D = Differential(t)

continuous_events = [
    [x ~ 0] => [vx ~ -vx]
    [y ~ -1.5, y ~ 1.5] => [vy ~ -vy]
]

@named ball = ODESystem([
    D(x) ~ vx
    D(y) ~ vy
    D(vx) ~ -9.8
    D(vy) ~ -0.01vy # there is some small air resistance
], t, continuous_events = continuous_events)



ball = structural_simplify(ball)

tspan = (0.0,5.0)
prob = ODEProblem(ball, Pair[], tspan)


cb = get_callback(prob)
@test cb isa ModelingToolkit.DiffEqCallbacks.VectorContinuousCallback
@test getfield(ball, :continuous_events)[1] == SymbolicContinuousCallback(Equation[x ~ 0], Equation[vx ~ -vx])
@test getfield(ball, :continuous_events)[2] == SymbolicContinuousCallback(Equation[y ~ -1.5, y ~ 1.5], Equation[vy ~ -vy])
cond = cb.condition
out = [0.0, 0.0, 0.0]
cond.rf_ip(out, [0,0,0,0], p0, t0)
@test out ≈ [0, 1.5, -1.5] 


sol = solve(prob,Tsit5())
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
    [x ~ 0] => [vx ~ -vx, vy ~ -vy]
]

@named ball = ODESystem([
    D(x) ~ vx
    D(y) ~ vy
    D(vx) ~ -1
    D(vy) ~ 0
], t, continuous_events = continuous_events)

ball = structural_simplify(ball)

tspan = (0.0,5.0)
prob = ODEProblem(ball, Pair[], tspan)
sol = solve(prob,Tsit5())
@test 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close
@test -minimum(sol[y]) ≈ maximum(sol[y]) ≈ sqrt(2)  # the ball will never go further than √2 in either direction (gravity was changed to 1 to get this particular number)

# tv = sort([LinRange(0, 5, 200); sol.t])
# plot(sol(tv)[y], sol(tv)[x], line_z=tv)
# vline!([-1.5, 1.5], l=(:black, 5), primary=false)
# hline!([0], l=(:black, 5), primary=false)