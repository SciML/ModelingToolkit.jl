using ModelingToolkit
using Test

@IVar t
@DVar x(t) y(t) z(t)
@test_throws AssertionError @DVar w

@Deriv D'~t # Default of first derivative, Derivative(t,1)
@Param σ ρ β
@Const c=0

# Default values
p = Parameter(:p, 1)
u = DependentVariable(:u, [1], dependents = [t])

s = JumpVariable(:s,3,dependents=[t])
n = NoiseVariable(:n,dependents=[t])

@test D(t) == @term(1)
