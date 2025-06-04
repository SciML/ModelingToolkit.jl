using ModelingToolkit, OrdinaryDiffEq, StochasticDiffEq
using Test, LinearAlgebra


# Change of variables: z = log(x)
# (this implies that x = exp(z) is automatically non-negative)
@independent_variables t
@variables z(t)[1:2, 1:2]
D = Differential(t)
eqs = [D(D(z)) ~ ones(2, 2)]
@mtkcompile sys = System(eqs, t)
@test_nowarn ODEProblem(sys, [z => zeros(2, 2), D(z) => ones(2, 2)], (0.0, 10.0))

@parameters α
@variables x(t)
D = Differential(t)
eqs = [D(x) ~ α*x]

tspan = (0., 1.)
def = [x => 1.0, α => -0.5]

@mtkcompile sys = System(eqs, t;defaults=def)
prob = ODEProblem(sys, [], tspan)
sol = solve(prob, Tsit5())

@variables z(t)
forward_subs  = [log(x) => z]
backward_subs = [x => exp(z)]
new_sys = changeofvariables(sys, t, forward_subs, backward_subs)
@test equations(new_sys)[1] == (D(z) ~ α)

new_prob = ODEProblem(new_sys, [], tspan)
new_sol = solve(new_prob, Tsit5())

@test isapprox(new_sol[x][end], sol[x][end], atol=1e-4)



# Riccati equation
@parameters α
@variables x(t)
D = Differential(t)
eqs = [D(x) ~ t^2 + α - x^2]
def = [x=>1., α => 1.]
@mtkcompile sys = System(eqs, t; defaults=def)

@variables z(t)
forward_subs  = [t + α/(x+t) =>    z       ]
backward_subs = [       x    => α/(z-t) - t]

new_sys = changeofvariables(sys, t, forward_subs, backward_subs; simplify=true, t0=0.)
# output should be equivalent to
# t^2 + α - z^2 + 2   (but this simplification is not found automatically)

tspan = (0., 1.)
prob = ODEProblem(sys,[],tspan)
new_prob = ODEProblem(new_sys,[],tspan)

sol = solve(prob, Tsit5())
new_sol = solve(new_prob, Tsit5())

@test isapprox(sol[x][end], new_sol[x][end], rtol=1e-4)


# # Linear transformation to diagonal system
# @variables x(t)[1:3]
# D = Differential(t)
# A = [0. -1. 0.; -0.5 0.5 0.; 0. 0. -1.]
# right = A.*transpose(x)
# eqs = [D(x[1]) ~ sum(right[1, 1:3]), D(x[2]) ~ sum(right[2, 1:3]), D(x[3]) ~ sum(right[3, 1:3])]

# tspan = (0., 10.)
# u0 = [x[1] => 1.0, x[2] => 2.0, x[3] => -1.0]

# @mtkcompile sys = System(eqs, t; defaults=u0)
# prob = ODEProblem(sys,[],tspan)
# sol = solve(prob, Tsit5())

# T = eigen(A).vectors

# @variables z(t)[1:3]
# forward_subs  = T \ x .=> z
# backward_subs = x     .=> T*z

# new_sys = changeofvariables(sys, t, forward_subs, backward_subs; simplify=true)

# new_prob = ODEProblem(new_sys, [], tspan)
# new_sol = solve(new_prob, Tsit5())

# # test RHS
# new_rhs = [eq.rhs for eq in equations(new_sys)]
# new_A = Symbolics.value.(Symbolics.jacobian(new_rhs, z))
# @test isapprox(diagm(eigen(A).values), new_A, rtol = 1e-10)
# @test isapprox( new_sol[x[1],end], sol[x[1],end], rtol=1e-4)

# Change of variables for sde
@independent_variables t
@brownian B
@parameters μ σ
@variables x(t) y(t)
D = Differential(t)
eqs = [D(x) ~ μ*x + σ*x*B]

def = [x=>0., μ => 2., σ=>1.]
@mtkcompile sys = System(eqs, t; defaults=def)
forward_subs = [log(x) => y]
backward_subs = [x => exp(y)]
new_sys = change_of_variable_SDE(sys, t, forward_subs, backward_subs)
@test equations(new_sys)[1] == (D(y) ~ μ - 1/2*σ^2)
@test ModelingToolkit.get_noise_eqs(new_sys)[1] === ModelingToolkit.value(σ)


