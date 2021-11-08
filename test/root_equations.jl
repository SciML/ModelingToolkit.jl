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