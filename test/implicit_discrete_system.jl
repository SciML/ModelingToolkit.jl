using ModelingToolkit, Test
using ModelingToolkit: t_nounits as t

k = ShiftIndex(t)
@variables x(t) = 1
@mtkbuild sys = ImplicitDiscreteSystem([x(k) ~ x(k)*x(k-1) - 3], t)
tspan = (0, 10)

# Shift(t, -1)(x(t)) - x_{t-1}(t)
# -3 - x(t) + x(t)*x_{t-1}
f = ImplicitDiscreteFunction(sys)
u_next = [3., 1.5]
@test f(u_next, [2.,3.], [], t) ≈ [0., 0.]
u_next = [0., 0.]
@test f(u_next, [2.,3.], [], t) ≈ [3., -3.]

resid = rand(2)
f(resid, u_next, [2.,3.], [], t)
@test resid ≈ [3., -3.]

# Initialization cases.
prob = ImplicitDiscreteProblem(sys, [x(k-1) => 3.], tspan)
@test prob.u0 == [3., 1.]
prob = ImplicitDiscreteProblem(sys, [], tspan)
@test prob.u0 == [1., 1.]
@variables x(t)
@mtkbuild sys = ImplicitDiscreteSystem([x(k) ~ x(k)*x(k-1) - 3], t)
@test_throws ErrorException prob = ImplicitDiscreteProblem(sys, [], tspan)

# Test solvers
