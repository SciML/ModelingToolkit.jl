using ModelingToolkit
using Test

@Param t σ ρ β
@Unknown x(t) y(t) z(t)
@Deriv D'~t # Default of first derivative, Derivative(t,1)
@Const c=0

# Default values
p = Parameter(:p)
u = Unknown(:u, dependents = [t])

σ*(y-x)
D(x)
D(x) ~ -σ*(y-x)
D(y) ~ x*(ρ-z)-sin(y)

@test D(t) == Constant(1)
