using ModelingToolkit
using Test

# Derivatives
@Param t σ ρ β
@Unknown x(t) y(t) z(t)
@Deriv D'~t D2''~t

@test expand_derivatives(D(t)) == 1
@test expand_derivatives(D(D(t))) == 0

dsin = D(sin(t))
@test expand_derivatives(dsin) == cos(t)

dcsch = D(csch(t))
@test expand_derivatives(dcsch) == simplify_constants(coth(t) * csch(t) * -1)

@test expand_derivatives(D(-7)) == 0
@test expand_derivatives(D(sin(2t))) == simplify_constants(cos(2t) * 2)
@test expand_derivatives(D2(sin(t))) == simplify_constants(-sin(t))
@test expand_derivatives(D2(sin(2t))) == simplify_constants(sin(2t) * -4)
@test expand_derivatives(D2(t)) == 0
@test expand_derivatives(D2(5)) == 0

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
jac = calculate_jacobian(sys)
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
@test D(b) ≠ 0

@test expand_derivatives(D(x * y)) == simplify_constants(y*D(x) + x*D(y))
@test_broken expand_derivatives(D(x * y)) == simplify_constants(D(x)*y + x*D(y))

@test expand_derivatives(D(2t)) == 2
@test expand_derivatives(D(2x)) == 2D(x)
