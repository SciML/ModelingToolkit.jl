using SciCompDSL
using Base.Test

# Derivatives
@IVar t
@Var x y z
@Param σ ρ β
@Deriv D'~t
dsin = D*sin(t)
expand_derivatives(dsin)

@test expand_derivatives(dsin) == cos(t)
dcsch = D*csch(t)
@test expand_derivatives(dcsch) == Operation(-coth(t)*csch(t))
# Chain rule
dsinsin = D*sin(sin(t))
@test expand_derivatives(dsinsin) == cos(sin(t))*cos(t)
# Binary
dpow1 = Derivative(^,[x, y],Val{1})
dpow2 = Derivative(^,[x, y],Val{2})
@test dpow1 == y*x^(y-1)
@test dpow2 == x^y*log(x)


eqs = [D*x ~ σ*(y-x),
       D*y ~ x*(ρ-z)-y,
       D*z ~ x*y - β*z]
de = NonlinearSystem(eqs,[x,y,z],[σ,ρ,β])
jac = SciCompDSL.generate_nlsys_jacobian(de)
@test_broken jac[1,1] == -σ
@test_broken jac[1,2] == σ
@test_broken jac[1,3] == 0
@test jac[2,1] == ρ-z
@test_broken jac[2,2] == -1
@test_broken jac[2,3] == -x
@test jac[3,1] == y
@test_broken jac[3,2] == x
@test_broken jac[3,3] == -β
