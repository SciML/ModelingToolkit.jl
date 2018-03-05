using SciCompDSL
using Base.Test

# `Expr`, `Number` -> `Operation`
@test isequal(Operation(2), 2)
expr = :(-inv(2sqrt(+(a, b))))
op   = Operation(-, [Operation(inv,
                [Operation(*, [2, Operation(sqrt,
                              [Operation(+, [:a, :b])])])])])
expr_post = Expr(:call, :-, op.args...)
@test isequal(Operation(expr), op)
@test isequal(Expr(Operation(expr)), expr_post)
expr1 = :(x^(y-1))
op1   = Operation(:^, [:x, Operation(-, [:y, 1])])
@test isequal(Operation(expr1), op1)
