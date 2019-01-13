using ModelingToolkit
using Test

@Param t
@Unknown x(t) y(t) z(t)

null_op = 0*t
@test simplify_constants(null_op) == Constant(0)

one_op = 1*t
@test simplify_constants(one_op) == t

identity_op = Operation(identity,[x])
@test simplify_constants(identity_op) == x

minus_op = -x
@test simplify_constants(minus_op) == -1*x
simplify_constants(minus_op)
