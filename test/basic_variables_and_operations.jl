using ModelingToolkit
using Test

@IVar t
@DVar x(t) y(t) z(t)
@test_throws AssertionError @DVar w

@Deriv D'~t # Default of first derivative, Derivative(t,1)
@Param σ ρ β
@Const c=0

# Default values
p = Parameter(:p)
u = DependentVariable(:u, dependents = [t])

s = JumpVariable(:s, dependents=[t])
n = NoiseVariable(:n, dependents=[t])

σ*(y-x)
D(x)
D(x) ~ -σ*(y-x)
D(y) ~ x*(ρ-z)-sin(y)

@test D(t) == Constant(1)
