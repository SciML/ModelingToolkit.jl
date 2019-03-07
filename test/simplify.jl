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
@test isequal(simplify_constants(minus_op), -1*x)
simplify_constants(minus_op)
