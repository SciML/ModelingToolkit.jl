using ModelingToolkit
using Test

# Derivatives
@Param t σ ρ β
@Unknown x(t) y(t) z(t)
@Deriv D'~t
dsin = D(sin(t))
expand_derivatives(dsin)

@test expand_derivatives(dsin) == cos(t)
dcsch = D(csch(t))
@test expand_derivatives(dcsch) == simplify_constants(coth(t) * csch(t) * -1)

# Chain rule
dsinsin = D(sin(sin(t)))
@test expand_derivatives(dsinsin) == cos(sin(t))*cos(t)

d1 = D(sin(t)*t)
d2 = D(sin(t)*cos(t))
@test expand_derivatives(d1) == t*cos(t)+sin(t)
@test expand_derivatives(d2) == simplify_constants(cos(t)*cos(t)+sin(t)*(-1*sin(t)))

eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
sys = NonlinearSystem(eqs,[x,y,z],[σ,ρ,β])
jac = ModelingToolkit.calculate_jacobian(sys)
@test jac[1,1] == σ*-1
@test jac[1,2] == σ
@test jac[1,3] == 0
@test jac[2,1] == ρ-z
@test jac[2,2] == -1
@test jac[2,3] == x*-1
@test jac[3,1] == y
@test jac[3,2] == x
@test jac[3,3] == -1*β

# Variable dependence checking in differentiation
@Unknown a(t) b(a)
@test D(b) ≠ Constant(0)

