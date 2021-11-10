using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: EqAffectPair, NULL_AFFECT

@parameters t
@variables x(t)=0 
D = Differential(t)

eqs = [D(x) ~ 1]
affect = [x ~ 0]

## Test EqAffectPair
@testset "EqAffectPair constructors" begin
    e = EqAffectPair(eqs[]) 
    @test e isa EqAffectPair
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT

    e = EqAffectPair(eqs) 
    @test e isa EqAffectPair
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT

    e = EqAffectPair(eqs, NULL_AFFECT) 
    @test e isa EqAffectPair
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT

    e = EqAffectPair(eqs[], NULL_AFFECT) 
    @test e isa EqAffectPair
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT

    e = EqAffectPair(eqs => NULL_AFFECT) 
    @test e isa EqAffectPair
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT

    e = EqAffectPair(eqs[] => NULL_AFFECT) 
    @test e isa EqAffectPair
    @test isequal(e.eqs, eqs)
    @test e.affect == NULL_AFFECT

    ## With affect

    e = EqAffectPair(eqs[], affect) 
    @test e isa EqAffectPair
    @test isequal(e.eqs, eqs)
    @test e.affect == affect

    e = EqAffectPair(eqs, affect) 
    @test e isa EqAffectPair
    @test isequal(e.eqs, eqs)
    @test e.affect == affect

    e = EqAffectPair(eqs, affect) 
    @test e isa EqAffectPair
    @test isequal(e.eqs, eqs)
    @test e.affect == affect

    e = EqAffectPair(eqs[], affect) 
    @test e isa EqAffectPair
    @test isequal(e.eqs, eqs)
    @test e.affect == affect

    e = EqAffectPair(eqs => affect) 
    @test e isa EqAffectPair
    @test isequal(e.eqs, eqs)
    @test e.affect == affect

    e = EqAffectPair(eqs[] => affect) 
    @test e isa EqAffectPair
    @test isequal(e.eqs, eqs)
    @test e.affect == affect
end


##

@named sys = ODESystem(eqs, root_eqs = [x ~ 1])
@test getfield(sys, :root_eqs)[] == EqAffectPair(Equation[x ~ 1], NULL_AFFECT)
@test isequal(equations(getfield(sys, :root_eqs))[], x ~ 1)
fsys = flatten(sys)
@test isequal(equations(getfield(fsys, :root_eqs))[], x ~ 1)

@named sys2 = ODESystem([D(x) ~ 1], root_eqs = [x ~ 2], systems=[sys])
@test getfield(sys2, :root_eqs)[] == EqAffectPair(Equation[x ~ 2], NULL_AFFECT)
@test all(ModelingToolkit.root_eqs(sys2) .== [EqAffectPair(Equation[x ~ 2], NULL_AFFECT), EqAffectPair(Equation[sys.x ~ 1], NULL_AFFECT)])

@test isequal(equations(getfield(sys2, :root_eqs))[1], x ~ 2)
@test length(ModelingToolkit.root_eqs(sys2)) == 2
@test isequal(ModelingToolkit.root_eqs(sys2)[1].eqs[], x ~ 2)
@test isequal(ModelingToolkit.root_eqs(sys2)[2].eqs[], sys.x ~ 1)

# Functions should be generated for root-finding equations
prob = ODEProblem(sys, Pair[], (0.0, 2.0))
p0 = 0
t0 = 0
@test prob.kwargs[:callback] isa ModelingToolkit.DiffEqCallbacks.ContinuousCallback
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


prob = ODEProblem(sys2, Pair[], (0.0, 3.0))
cbs = prob.kwargs[:callback]
@test cbs isa ModelingToolkit.DiffEqCallbacks.CallbackSet

cb = cbs.continuous_callbacks[1]
cond = cb.condition
out = [0.0]
# the root to find is 2
cond.rf_ip(out, [0], p0, t0)
@test out[] ≈ -2 # signature is u,p,t
cond.rf_ip(out, [1], p0, t0)
@test out[] ≈ -1  # signature is u,p,t
cond.rf_ip(out, [2], p0, t0) # this should return 0
@test out[] ≈ 0  # signature is u,p,t

# the root to find is 1
cb = cbs.continuous_callbacks[2]
cond = cb.condition
out = [0.0]
cond.rf_ip(out, [0], p0, t0)
@test out[] ≈ -1 # signature is u,p,t
cond.rf_ip(out, [1], p0, t0) # this should return 0
@test out[] ≈ 0  # signature is u,p,t
cond.rf_ip(out, [2], p0, t0)
@test out[] ≈ 1  # signature is u,p,t

sol = solve(prob, Tsit5())
@test minimum(t->abs(t-1), sol.t) < 1e-10 # test that the solver stepped at the first root
@test minimum(t->abs(t-2), sol.t) < 1e-10 # test that the solver stepped at the second root

@named sys = ODESystem(eqs, root_eqs = [x ~ 1, x ~ 2]) # two root eqs using the same state
prob = ODEProblem(sys, Pair[], (0.0, 3.0))
@test prob.kwargs[:callback] isa ModelingToolkit.DiffEqCallbacks.VectorContinuousCallback
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
], t, root_eqs = root_eqs => affect)

@test getfield(ball, :root_eqs)[] == EqAffectPair(Equation[x ~ 0], Equation[v ~ -v])
ball = structural_simplify(ball)

tspan = (0.0,5.0)
prob = ODEProblem(ball, Pair[], tspan)
sol = solve(prob,Tsit5())
@test 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close
# plot(sol)
