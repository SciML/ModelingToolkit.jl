using ModelingToolkit
using Test

@Param t()
@Unknown w
@Unknown x(t)
@Unknown y(t)
@Unknown z(t)
@Unknown g(3-t, t^2)
w1 = Unknown(:w)
x1 = Unknown(:x)(t)
y1 = Unknown(:y)(t)
z1 = Unknown(:z)(t)
g1 = Unknown(:g)(3-t, t^2)
@test w1 == w
@test w1(t) == w(t)
@test x1 == x
@test y1 == y
@test z1 == z
@test g1 == g
@test convert(Expr, w) == :w
@test convert(Expr, w(t)) == :(w(t()))
@test convert(Expr, x) == :(x(t()))
@test convert(Expr, y) == :(y(t()))
@test convert(Expr, z) == :(z(t()))
@test convert(Expr, g) == :(g(3-t(), t()^2))

@Param begin
    t()
    s()
end
t1 = Parameter(:t)()
s1 = Parameter(:s)()
@test t1 == t
@test s1 == s
@test convert(Expr, t) == :(t())
@test convert(Expr, s) == :(s())
@test convert(Expr, cos(t + sin(s))) == :(cos(t() + sin(s())))

@Deriv D''~t
D1 = Differential(t, 2)
@test D1 == D
@test convert(Expr, D) == D

@Const c=0 v=2
c1 = Constant(0)
v1 = Constant(2)
@test c1 == c
@test v1 == v
@test convert(Expr, c) == 0
@test convert(Expr, v) == 2
