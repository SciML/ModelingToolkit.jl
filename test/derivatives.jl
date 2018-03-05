using SciCompDSL
using Base.Test

# Derivatives
@IVar t
@Deriv D'~t
dsin = D*sin(t)
@test isequal(expand_derivatives(dsin),cos(t))
dcsch = D*csch(t)
@test isequal(expand_derivatives(dcsch),Operation(-coth(t)*csch(t)))
# Chain rule
dsinsin = D*sin(sin(t))
@test isequal(expand_derivatives(dsinsin),cos(sin(t))*cos(t))
# Binary
dpow1 = Derivative(^,[:x, :y],Val{1})
dpow2 = Derivative(^,[:x, :y],Val{2})
@test isequal(dpow1, Operation( :(y*x^(y-1)) ) )
@test isequal(dpow2, Operation( :(x^y*log(x)) ) )
