using ModelingToolkit
using Test

@param t
@variable x(t)
@variable y(t)
@variable z(t)
x1 = Variable(:x, [t])
y1 = Variable(:y, [t])
z1 = Variable(:z, [t])
@test isequal(x1, x)
@test isequal(y1, y)
@test isequal(z1, z)
@test convert(Expr, x) == :x
@test convert(Expr, y) == :y
@test convert(Expr, z) == :z

@param begin
    t
    s
end
t1 = Variable(:t; known = true)
s1 = Variable(:s; known = true)
@test isequal(t1, t)
@test isequal(s1, s)
@test convert(Expr, t) == :t
@test convert(Expr, s) == :s
@test convert(Expr, cos(t + sin(s))) == :(cos(t + sin(s)))

@deriv D'~t
D1 = Differential(t)
@test D1 == D
@test convert(Expr, D) == D

@test isequal(x ≤ y + 1, (x < y + 1) | (x == y + 1))
