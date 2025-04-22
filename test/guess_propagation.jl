using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: D, t_nounits as t
using Test

# Standard case

@variables x(t) [guess = 2]
@variables y(t)
eqs = [D(x) ~ 1
       x ~ y]
initialization_eqs = [1 ~ exp(1 + x)]

@named sys = ODESystem(eqs, t; initialization_eqs)
sys = complete(structural_simplify(sys))
tspan = (0.0, 0.2)
prob = ODEProblem(sys, [], tspan, [])

@test prob.f.initializeprob[y] == 2.0
@test prob.f.initializeprob[x] == 2.0
sol = solve(prob.f.initializeprob; show_trace = Val(true))

# Guess via observed

@variables x(t)
@variables y(t) [guess = 2]
eqs = [D(x) ~ 1
       x ~ y]
initialization_eqs = [1 ~ exp(1 + x)]

@named sys = ODESystem(eqs, t; initialization_eqs)
sys = complete(structural_simplify(sys))
tspan = (0.0, 0.2)
prob = ODEProblem(sys, [], tspan, [])

@test prob.f.initializeprob[x] == 2.0
@test prob.f.initializeprob[y] == 2.0
sol = solve(prob.f.initializeprob; show_trace = Val(true))

# Guess via parameter

@parameters a = -1.0
@variables x(t) [guess = a]

eqs = [D(x) ~ a]

initialization_eqs = [1 ~ exp(1 + x)]

@named sys = ODESystem(eqs, t; initialization_eqs)
sys = complete(structural_simplify(sys))

tspan = (0.0, 0.2)
prob = ODEProblem(sys, [], tspan, [])

@test prob.f.initializeprob[x] == -1.0
sol = solve(prob.f.initializeprob; show_trace = Val(true))

# Guess via observed parameter

@parameters a = -1.0
@variables x(t)
@variables y(t) [guess = a]

eqs = [D(x) ~ a,
    y ~ x]

initialization_eqs = [1 ~ exp(1 + x)]

@named sys = ODESystem(eqs, t; initialization_eqs)
sys = complete(structural_simplify(sys))

tspan = (0.0, 0.2)
prob = ODEProblem(sys, [], tspan, [])

@test prob.f.initializeprob[x] == -1.0
sol = solve(prob.f.initializeprob; show_trace = Val(true))

# Test parameters + defaults
# https://github.com/SciML/ModelingToolkit.jl/issues/2774

@parameters x0
@variables x(t)
@variables y(t) = x
@mtkbuild sys = ODESystem([x ~ x0, D(y) ~ x], t)
prob = ODEProblem(sys, [], (0.0, 1.0), [x0 => 1.0])
@test prob[x] == 1.0
@test prob[y] == 1.0

@parameters x0
@variables x(t)
@variables y(t) = x0
@mtkbuild sys = ODESystem([x ~ x0, D(y) ~ x], t)
prob = ODEProblem(sys, [], (0.0, 1.0), [x0 => 1.0])
@test prob[x] == 1.0
@test prob[y] == 1.0

@parameters x0
@variables x(t)
@variables y(t) = x0
@mtkbuild sys = ODESystem([x ~ y, D(y) ~ x], t)
prob = ODEProblem(sys, [], (0.0, 1.0), [x0 => 1.0])
@test prob[x] == 1.0
@test prob[y] == 1.0

@parameters x0
@variables x(t) = x0
@variables y(t) = x
@mtkbuild sys = ODESystem([x ~ y, D(y) ~ x], t)
prob = ODEProblem(sys, [], (0.0, 1.0), [x0 => 1.0])
@test prob[x] == 1.0
@test prob[y] == 1.0
