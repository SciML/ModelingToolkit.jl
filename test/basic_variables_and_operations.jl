using SciCompDSL
using Base.Test

@DVar x y z
@IVar t
@Deriv D'~t # Default of first derivative, Derivative(t,1)
@Param σ ρ β
@Const c=0

# Default values
p = Parameter(:p, 1)
u = DependentVariable(:u, [1])

s = JumpVariable(:s,3)
n = NoiseVariable(:n)

σ*(y-x)
D*x
D*x == -σ*(y-x)
D*y == x*(ρ-z)-sin(y)

@test isequal(D*t,Constant(1))
null_op = 0*t
@test isequal(simplify_constants(null_op),Constant(0))
