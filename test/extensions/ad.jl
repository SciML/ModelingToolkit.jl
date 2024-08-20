using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Zygote
using SymbolicIndexingInterface
using SciMLStructures
using OrdinaryDiffEq
using SciMLSensitivity

@variables x(t)[1:3] y(t)
@parameters p[1:3, 1:3] q
eqs = [
         D(x) ~ p * x
         D(y) ~ sum(p) + q * y
]
u0 = [x => zeros(3),
      y => 1.]
ps = [p => zeros(3, 3),
    q => 1.]
tspan = (0., 10.)
@mtkbuild sys = ODESystem(eqs, t)
prob = ODEProblem(sys, u0, tspan, ps)
sol = solve(prob, Tsit5())

mtkparams = parameter_values(prob)
new_p = rand(10)
gs = gradient(new_p) do new_p
    new_params = SciMLStructures.replace(SciMLStructures.Tunable(), mtkparams, new_p)
    new_prob = remake(prob, p = new_params)
    new_sol = solve(new_prob, Tsit5())
    sum(new_sol)
end
