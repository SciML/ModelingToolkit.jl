using SciCompDSL
using Base.Test

@DVar x y=sin(1)+exp(1) z
x1 = DependentVariable(:x)
y1 = DependentVariable(:y, sin(1) + exp(1))
z1 = DependentVariable(:z)
@test isequal(x1, x)
@test isequal(y1, y)
@test isequal(z1, z)
@IVar begin
    t
    s = cos(2.5)
end
t1 = IndependentVariable(:t)
s1 = IndependentVariable(:s, cos(2.5))
@test isequal(t1, t)
@test isequal(s1, s)
@Deriv D''~t
D1 = Differential(t, 2)
@test isequal(D1, D)
@Const c=0 v=2
c1 = Constant(0)
v1 = Constant(2)
@test isequal(c1, c)
@test isequal(v1, v)
