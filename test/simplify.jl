using ModelingToolkit
using Test

@IVar t
@DVar x(t) y(t) z(t)

null_op = @term(0*t)
@test simplify_constants(null_op) == @term(0)

one_op = @term(1*t)
@test simplify_constants(one_op) == @term(t)

identity_op = @term(identity(x))
@test simplify_constants(identity_op) == @term(x)

minus_op = @term(-x)
@test simplify_constants(minus_op) == @term(-1*x)

fold_const_op = @term(x * 2 * 3 * y)
@test simplify_constants(fold_const_op) == simplify_constants(@term(x * y * 6))

@test_broken simplify_constants(@term(x * 2)) == simplify_constants(@term(2 * x))
