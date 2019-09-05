using Distributed
# add processes to workspace
addprocs(2)

using ModelingToolkit
using OrdinaryDiffEq

# create the Lorenz system
@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

de = ODESystem(eqs)
ode_func = ODEFunction(de, [x,y,z], [σ, ρ, β])

u0 = [19.,20.,50.]
params = [16.,45.92,4]

ode_prob = ODEProblem(ode_func, u0, (0., 10.),params)

@everywhere begin

    using DifferentialEquations
    using ModelingToolkit

    function solve_lorenz(ode_problem)
        print(solve(ode_problem,Tsit5()))
    end
end

solve_lorenz(ode_prob)

future = @spawn solve_lorenz(ode_prob)
fetch(future)
