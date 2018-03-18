using ModelingToolkit
using Base.Test

# `Expr`, `Number` -> `Operation`
@IVar a
@Param b
@DVar x(t)
@Var y
@test parse(Operation, 2) == 2
expr = :(-inv(2sqrt(+($a, $b))))
op   = Operation(-, [Operation(inv,
                [Operation(*, [2, Operation(sqrt,
                              [Operation(+, [a, b])])])])])
@test parse(Operation, expr) == op
expr1 = :($x^($y-1))
op1   = Operation(^, [x, Operation(-, [y, 1])])
@test parse(Operation, expr1) == op1
