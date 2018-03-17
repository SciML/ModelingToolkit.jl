using SciCompDSL
using Base.Test

@IVar t
@DVar x(t)
macroexpand(:(@DVar x(t)))
macroexpand(:(@IVar t))
@DVar x y z

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
D*x ~ -σ*(y-x)
D*y ~ x*(ρ-z)-sin(y)

@test D*t == Constant(1)
null_op = 0*t
@test simplify_constants(null_op) == Constant(0)

macro escape_all(x...)
    :($(esc(x))...)
end
x = 1
y = 2
z = 3
macroexpand(:(@escape_all x y z))
