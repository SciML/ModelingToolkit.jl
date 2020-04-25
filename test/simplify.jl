using ModelingToolkit
using Test

@parameters t
@variables x(t) y(t) z(t)

null_op = 0*t
@test isequal(simplify_constants(null_op), 0)

one_op = 1*t
@test isequal(simplify_constants(one_op), t)

identity_op = Operation(identity,[x])
@test isequal(simplify_constants(identity_op), x)

minus_op = -x
@test isequal(simplify_constants(minus_op), -x)
simplify_constants(minus_op)

@variables x

@test simplified_expr(expand_derivatives(Differential(x)((x-2)^2))) == :(2 * (-2 + x))
@test simplified_expr(expand_derivatives(Differential(x)((x-2)^3))) == :(3 * (-2 + x)^2)
@test simplified_expr(simplify_constants(x+2+3)) == :(5 + x)

d1 = Differential(x)((-2 + x)^2)
d2 = Differential(x)(d1)
d3 = Differential(x)(d2)

@test simplified_expr(expand_derivatives(d3)) == :(0)
#@test simplify_constants(x^0) == 1
