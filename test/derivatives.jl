using SciCompDSL
using Base.Test

# Derivatives
@IVar t
@Var x y z
@Param σ ρ β
@Deriv D'~t
dsin = D*sin(t)
expand_derivatives(dsin)

@test isequal(expand_derivatives(dsin),cos(t))
dcsch = D*csch(t)
@test isequal(expand_derivatives(dcsch),Operation(-coth(t)*csch(t)))
# Chain rule
dsinsin = D*sin(sin(t))
@test isequal(expand_derivatives(dsinsin),cos(sin(t))*cos(t))
# Binary
dpow1 = Derivative(^,[x, y],Val{1})
dpow2 = Derivative(^,[x, y],Val{2})
@test isequal(dpow1, y*x^(y-1))
@test isequal(dpow2, x^y*log(x))


eqs = [D*x == σ*(y-x),
       D*y == x*(ρ-z)-y,
       D*z == x*y - β*z]
de = DiffEqSystem(eqs,[t],[x,y,z],Variable[],[σ,ρ,β])
jac = SciCompDSL.generate_ode_jacobian(de,false)
@test_broken isequal(jac[1,1],-σ)
@test_broken isequal(jac[1,2],σ)
@test_broken isequal(jac[1,3],0)
@test_broken isequal(jac[2,1],ρ)
@test_broken isequal(jac[2,2],-1)
@test_broken isequal(jac[2,3],-x)
@test_broken isequal(jac[3,1],y)
@test_broken isequal(jac[3,2],x)
@test_broken isequal(jac[3,3],-β)
