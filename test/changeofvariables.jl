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

tspan = (0.0, 1.0)
def = [x => 1.0, α => -0.5]

@mtkcompile sys = System(eqs, t; defaults = def)
prob = ODEProblem(sys, [], tspan)
sol = solve(prob, Tsit5())

@variables z(t)
forward_subs = [log(x) => z]
backward_subs = [x => exp(z)]
new_sys = change_of_variables(sys, t, forward_subs, backward_subs)
@test equations(new_sys)[1] == (D(z) ~ α)

new_prob = ODEProblem(new_sys, [], tspan)
new_sol = solve(new_prob, Tsit5())

@test isapprox(new_sol[x][end], sol[x][end], atol = 1e-4)

# Riccati equation
@parameters α
@variables x(t)
D = Differential(t)
eqs = [D(x) ~ t^2 + α - x^2]
def = [x=>1.0, α => 1.0]
@mtkcompile sys = System(eqs, t; defaults = def)

@variables z(t)
forward_subs = [t + α/(x+t) => z]
backward_subs = [x => α/(z-t) - t]

new_sys = change_of_variables(
    sys, t, forward_subs, backward_subs; simplify = true, t0 = 0.0)
# output should be equivalent to
# t^2 + α - z^2 + 2   (but this simplification is not found automatically)

tspan = (0.0, 1.0)
prob = ODEProblem(sys, [], tspan)
new_prob = ODEProblem(new_sys, [], tspan)

sol = solve(prob, Tsit5())
new_sol = solve(new_prob, Tsit5())

@test isapprox(sol[x][end], new_sol[x][end], rtol = 1e-4)

# Linear transformation to diagonal system
@independent_variables t
@variables x(t)[1:3]
x = reshape(x, 3, 1)
D = Differential(t)
A = [0.0 -1.0 0.0; -0.5 0.5 0.0; 0.0 0.0 -1.0]
right = A*x
eqs = vec(D.(x) .~ right)

tspan = (0.0, 10.0)
u0 = [x[1] => 1.0, x[2] => 2.0, x[3] => -1.0]

@mtkcompile sys = System(eqs, t; defaults = u0)
prob = ODEProblem(sys, [], tspan)
sol = solve(prob, Tsit5())

T = eigen(A).vectors
T_inv = inv(T)

@variables z(t)[1:3]
z = reshape(z, 3, 1)
forward_subs = vec(T_inv*x .=> z)
backward_subs = vec(x .=> T*z)

new_sys = change_of_variables(sys, t, forward_subs, backward_subs; simplify = true)

new_prob = ODEProblem(new_sys, [], tspan)
new_sol = solve(new_prob, Tsit5())

# test RHS
new_rhs = [eq.rhs for eq in equations(new_sys)]
new_A = Symbolics.value.(Symbolics.jacobian(new_rhs, z))
A = diagm(eigen(A).values)
A = sortslices(A, dims = 1)
new_A = sortslices(new_A, dims = 1)
@test isapprox(A, new_A, rtol = 1e-10)
@test isapprox(new_sol[x[1], end], sol[x[1], end], rtol = 1e-4)

# Change of variables for sde
noise_eqs = ModelingToolkit.get_noise_eqs
value = ModelingToolkit.value

@independent_variables t
@brownians B
@parameters μ σ
@variables x(t) y(t)
D = Differential(t)
eqs = [D(x) ~ μ*x + σ*x*B]

def = [x=>0.0, μ => 2.0, σ=>1.0]
@mtkcompile sys = System(eqs, t; defaults = def)
forward_subs = [log(x) => y]
backward_subs = [x => exp(y)]
new_sys = change_of_variables(sys, t, forward_subs, backward_subs)
@test equations(new_sys)[1] == (D(y) ~ μ - 1/2*σ^2)
@test noise_eqs(new_sys)[1] === value(σ)

#Multiple Brownian and equations
@independent_variables t
@brownians Bx By
@parameters μ σ α
@variables x(t) y(t) z(t) w(t) u(t) v(t)
D = Differential(t)
eqs = [D(x) ~ μ*x + σ*x*Bx, D(y) ~ α*By, D(u) ~ μ*u + σ*u*Bx + α*u*By]
def = [x=>0.0, y => 0.0, u=>0.0, μ => 2.0, σ=>1.0, α=>3.0]
forward_subs = [log(x) => z, y^2 => w, log(u) => v]
backward_subs = [x => exp(z), y => w^0.5, u => exp(v)]

@mtkcompile sys = System(eqs, t; defaults = def)
new_sys = change_of_variables(sys, t, forward_subs, backward_subs)
@test equations(new_sys)[1] == (D(z) ~ μ - 1/2*σ^2)
@test equations(new_sys)[2] == (D(w) ~ α^2)
@test equations(new_sys)[3] == (D(v) ~ μ - 1/2*(α^2 + σ^2))
@test noise_eqs(new_sys)[1, 1] === value(σ)
@test noise_eqs(new_sys)[1, 2] === value(0)
@test noise_eqs(new_sys)[2, 1] === value(0)
@test noise_eqs(new_sys)[2, 2] === value(simplify(substitute(2*α*y, backward_subs[2])))
@test noise_eqs(new_sys)[3, 1] === value(σ)
@test noise_eqs(new_sys)[3, 2] === value(α)

# Test for  Brownian instead of noise
@named sys = System(eqs, t; defaults = def)
new_sys = change_of_variables(sys, t, forward_subs, backward_subs; simplify = false)
@test simplify(equations(new_sys)[1]) == simplify((D(z) ~ μ - 1/2*σ^2 + σ*Bx))
@test simplify(equations(new_sys)[2]) == simplify((D(w) ~ α^2 + 2*α*w^0.5*By))
@test simplify(equations(new_sys)[3]) ==
      simplify((D(v) ~ μ - 1/2*(α^2 + σ^2) + σ*Bx + α*By))
