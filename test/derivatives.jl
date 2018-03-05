using SciCompDSL
using Base.Test

# Derivatives
@IVar t
@Deriv D'~t
dsin = D*sin(t)
@test isequal(expand_derivatives(dsin),cos(t))
dsinsin = D*sin(sin(t))
@test isequal(expand_derivatives(dsinsin),cos(sin(t))*cos(t))
