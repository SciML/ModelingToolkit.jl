using ModelingToolkit
using NonlinearSolve
using Test
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables x y z u
@parameters σ ρ β

eqs = [0 ~ σ * (y - x),
    0 ~ x * (ρ - z) - y,
    0 ~ x * y - β * z]

par = [
    σ => 1,
    ρ => 0.1 + σ,
    β => ρ * 1.1
]
u0 = [
    x => u,
    y => σ, # default u0 from default p
    z => u - 0.1
]
ns = System(eqs, [x, y, z], [σ, ρ, β], name = :ns, defaults = [par; u0])
ModelingToolkit.get_defaults(ns)[y] = u * 1.1
resolved = ModelingToolkit.varmap_to_vars(defaults(ns), parameters(ns))
@test resolved == [1, 0.1 + 1, (0.1 + 1) * 1.1]

prob = NonlinearProblem(complete(ns), [u => 1.0])
@test prob.u0 == [1.0, 1.1, 0.9]
sol = solve(prob, NewtonRaphson())

@variables a
@parameters b
top = System([0 ~ -a + ns.x + b], [a], [b], systems = [ns], name = :top)
ModelingToolkit.get_defaults(top)[b] = ns.σ * 0.5
ModelingToolkit.get_defaults(top)[ns.x] = unknowns(ns, u) * 0.5

res = ModelingToolkit.varmap_to_vars(defaults(top), parameters(top))
@test res == [0.5, 1, 0.1 + 1, (0.1 + 1) * 1.1]

top = complete(top)
prob = NonlinearProblem(top, [unknowns(ns, u) => 1.0, a => 1.0])
@test prob.u0 == [1.0, 0.5, 1.1, 0.9]
sol = solve(prob, NewtonRaphson())

# test NullParameters+defaults
prob = NonlinearProblem(top, [unknowns(ns, u) => 1.0, a => 1.0])
@test prob.u0 == [1.0, 0.5, 1.1, 0.9]
sol = solve(prob, NewtonRaphson())

# test initial conditions and parameters at the problem level
pars = @parameters(begin
    x0
end)
vars = @variables(begin
    x(ModelingToolkit.t_nounits)
end)
der = Differential(t)
eqs = [der(x) ~ x]
@named sys = System(eqs, t, vars, [x0])
sys = complete(sys)
initialValues = [x => x0
                 x0 => 10.0]
tspan = (0.0, 1.0)
problem = ODEProblem(sys, initialValues, tspan)
@test problem.u0 isa Vector{Float64}
