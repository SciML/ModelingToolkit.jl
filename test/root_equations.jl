using ModelingToolkit, OrdinaryDiffEq, Test

@parameters t
@variables x(t)=0 
D = Differential(t)

eqs = [D(x) ~ 1]

@named sys = ODESystem(eqs, root_eqs = [x ~ 1])
@test isequal(sys.root_eqs[], x ~ 1)
fsys = flatten(sys)
@test isequal(fsys.root_eqs[], x ~ 1)
prob = ODEProblem(sys, Pair[], (0.0, 2.0))

@named sys2 = ODESystem([D(x) ~ sys.x], root_eqs = [x ~ 2], systems=[sys])
@test isequal(sys2.root_eqs[], x ~ 2)
@test length(ModelingToolkit.root_eqs(sys2)) == 2
@test isequal(ModelingToolkit.root_eqs(sys2)[1], x ~ 2)
@test isequal(ModelingToolkit.root_eqs(sys2)[2], sys.x ~ 1)

# Functions should be generated for root-finding equations
@test prob.kwargs[:callback] isa ModelingToolkit.DiffEqCallbacks.ContinuousCallback
cb = ModelingToolkit.generate_rootfinding_callback(sys)
cond = cb.condition
out = [0.0]
cond.rf_ip(out, [0], 1.2, 1)
@test out[] ≈ -1 # signature is u,p,t
cond.rf_ip(out, [1], 1.2, 1)
@test out[] ≈ 0  # signature is u,p,t
cond.rf_ip(out, [2], 1.2, 1)
@test out[] ≈ 1  # signature is u,p,t


prob = ODEProblem(sys, Pair[], (0.0, 2.0))
sol = solve(prob, Tsit5())
@test minimum(t->abs(t-1), sol.t) < 1e-10 # test that the solver stepped at the root


prob = ODEProblem(sys2, Pair[], (0.0, 3.0))
@test prob.kwargs[:callback] isa ModelingToolkit.DiffEqCallbacks.VectorContinuousCallback
sol = solve(prob, Tsit5())
@test minimum(t->abs(t-1), sol.t) < 1e-10 # test that the solver stepped at the first root
@test minimum(t->abs(t-2), sol.t) < 1e-10 # test that the solver stepped at the second root

@named sys = ODESystem(eqs, root_eqs = [x ~ 1, x ~ 2]) # two root eqs using the same state
prob = ODEProblem(sys, Pair[], (0.0, 3.0))
@test prob.kwargs[:callback] isa ModelingToolkit.DiffEqCallbacks.VectorContinuousCallback
sol = solve(prob, Tsit5())
@test minimum(t->abs(t-1), sol.t) < 1e-10 # test that the solver stepped at the first root
@test minimum(t->abs(t-2), sol.t) < 1e-10 # test that the solver stepped at the second root
