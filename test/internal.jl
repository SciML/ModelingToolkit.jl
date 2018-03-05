using SciCompDSL
using Base.Test

# `Expr`, `Number` -> `Operation`
@IVar a b x y
@test isequal(parse(Operation, 2), 2)
expr = :(-inv(2sqrt(+($a, $b))))
op   = Operation(-, [Operation(inv,
                [Operation(*, [2, Operation(sqrt,
                              [Operation(+, [a, b])])])])])
@test isequal(parse(Operation, expr), op)
expr1 = :($x^($y-1))
op1   = Operation(^, [x, Operation(-, [y, 1])])
@test isequal(parse(Operation, expr1), op1)
