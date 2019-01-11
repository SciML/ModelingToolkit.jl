using ModelingToolkit
using Test

@Var a b
a1 = Variable(:a)
@test a1 == a
@test convert(Expr, a) == :a

@Var begin
    a
    b
end

@IVar t
@DVar x(t)
@DVar y(t)
@DVar z(t)
x1 = DependentVariable(:x ,dependents = [t])
y1 = DependentVariable(:y ,dependents = [t])
z1 = DependentVariable(:z ,dependents = [t])
@test x1 == x
@test y1 == y
@test z1 == z
@test convert(Expr, x) == :x
@test convert(Expr, y) == :y
@test convert(Expr, z) == :z

@IVar begin
    t
    s
end
t1 = IndependentVariable(:t)
s1 = IndependentVariable(:s)
@test t1 == t
@test s1 == s
@test convert(Expr, t) == :t
@test convert(Expr, s) == :s
@test convert(Expr, cos(t + sin(s))) == :(cos(t + sin(s)))

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
