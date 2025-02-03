using Distributed
# add processes to workspace
addprocs(2)

@everywhere using ModelingToolkit, OrdinaryDiffEq
@everywhere using ModelingToolkit: t_nounits as t, D_nounits as D

# create the Lorenz system
@everywhere @parameters σ ρ β
@everywhere @variables x(t) y(t) z(t)

@everywhere eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@everywhere @named de = ODESystem(eqs, t)
@everywhere de = complete(de)

@everywhere u0 = [19.0, 20.0, 50.0]
@everywhere params = [16.0, 45.92, 4]

@everywhere ode_prob = ODEProblem(de, u0, (0.0, 10.0), params)

@everywhere begin
    using OrdinaryDiffEq
    using ModelingToolkit

    function solve_lorenz(ode_problem)
        print(solve(ode_problem, Tsit5()))
    end
end

solve_lorenz(ode_prob)

future = @spawn solve_lorenz(ode_prob)
fetch(future)
