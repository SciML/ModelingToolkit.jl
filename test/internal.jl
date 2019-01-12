using ModelingToolkit
using Test

# `Expr`, `Number` -> `Operation`
@IVar a
@Param b
@DVar x(t) y()
@test convert(Expression, 2) == 2
expr = :(-inv(2sqrt(+($a, $b))))
op   = Operation(-, [Operation(inv,
                [Operation(*, [2, Operation(sqrt,
                              [Operation(+, [a, b])])])])])
@test convert(Expression, expr) == op
expr1 = :($x^($y-1))
op1   = Operation(^, [x, Operation(-, [y, 1])])
@test convert(Expression, expr1) == op1

@test_throws ArgumentError convert(Expression, :([a, b]))
@test_throws ArgumentError convert(Expression, :(a ? b : c))
