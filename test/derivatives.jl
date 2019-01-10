using ModelingToolkit
using Test

# Derivatives
@IVar t
@Var x(t) y(t) z(t)
@Param σ ρ β
@Deriv D'~t
dsin = @term(D(sin(t)))
expand_derivatives(dsin)

@test expand_derivatives(dsin) == @term(cos(t))
dcsch = @term(D(csch(t)))
@test expand_derivatives(dcsch) == simplify_constants(@term(coth(t) * csch(t) * -1))

# Chain rule
dsinsin = @term(D(sin(sin(t))))
@test expand_derivatives(dsinsin) == @term(cos(sin(t))*cos(t))
# Binary
dpow1 = Derivative(^,[x, y],Val(1))
dpow2 = Derivative(^,[x, y],Val(2))
@test dpow1 == @term(y * x ^ (y - 1))
@test dpow2 == @term(x ^ y * log(x))

d1 = @term(D(sin(t)*t))
d2 = @term(D(sin(t)*cos(t)))
@test expand_derivatives(d1) == @term(t*cos(t)+sin(t))
@test expand_derivatives(d2) == simplify_constants(@term(cos(t)*cos(t)+sin(t)*(-1*sin(t))))

eqs = [
    @term 0 ~ σ*(y-x)
    @term 0 ~ x*(ρ-z)-y
    @term 0 ~ x*y - β*z
]
sys = NonlinearSystem(eqs,[x,y,z],[σ,ρ,β])
jac = ModelingToolkit.calculate_jacobian(sys)
@test jac[1,1] == @term σ*-1
@test jac[1,2] == @term σ
@test jac[1,3] == @term 0
@test jac[2,1] == @term ρ-z
@test jac[2,2] == @term -1
@test jac[2,3] == @term x*-1
@test jac[3,1] == @term y
@test jac[3,2] == @term x
@test jac[3,3] == @term -1*β

# Variable dependence checking in differentiation
@Var a(t) b(a)
@test D(b) ≠ @term(0)
