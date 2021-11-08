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
cb = ModelingToolkit.generate_rootfinding_callback(sys)
cond = cb.condition
@test cond.rf([0], 1.2, 1.3) ≈ -1
@test cond.rf([1], 1.2, 1.3) ≈ 0
@test cond.rf([2], 1.2, 1.3) ≈ 1


prob = ODEProblem(sys, Pair[], (0.0, 2.0))
sol = solve(prob, Tsit5())
@test minimum(t->abs(t-1), sol.t) < 1e-10 # test that the solver stepped at the root


# TODO: the problem is that the ContinuousCallback only handles scalar conditions. NEed to adjust the condition function to always be multivariate and use VectorContinuousCallback
prob = ODEProblem(sys2, Pair[], (0.0, 2.0))
sol = solve(prob, Tsit5())
@test minimum(t->abs(t-1), sol.t) < 1e-10 # test that the solver stepped at the first root
@test minimum(t->abs(t-2), sol.t) < 1e-10 # test that the solver stepped at the second root