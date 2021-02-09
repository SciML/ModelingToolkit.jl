using ModelingToolkit
using NonlinearSolve
using Test

@variables x y z
@parameters σ ρ β

eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]

par = [
    σ => 1,
    ρ => 0.1+σ,
    β => ρ*1.1
]
ns = NonlinearSystem(eqs, [x,y,z],[σ,ρ,β], name=:ns, default_p=par)
ModelingToolkit.default_p(ns)
resolved = ModelingToolkit.varmap_to_vars(Dict(), parameters(ns), defaults=ModelingToolkit.default_p(ns))
@test resolved == [1, 0.1+1, (0.1+1)*1.1]

prob = NonlinearProblem(ns, ones(3), Pair[])
@show sol = solve(prob,NewtonRaphson())

@variables a
@parameters b
top = NonlinearSystem([0 ~ -a + ns.x+1], [a], [b], systems=[ns], name=:top)
flat = flatten(top)
flat.b = ns.σ*0.5

flatres = ModelingToolkit.varmap_to_vars(Dict(), parameters(flat), defaults=ModelingToolkit.default_p(flat))
@test flatres == [0.5, 1, 0.1+1, (0.1+1)*1.1]

prob = NonlinearProblem(flat, ones(4), Pair[])
@show sol = solve(prob,NewtonRaphson())