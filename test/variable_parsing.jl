using ModelingToolkit
using Base.Test

@Var a=1.0 b
a1 = Variable(:a,1.0)
@test a1 == a

@Var begin
    a = 1.0
    b
end

@IVar t
@DVar x(t)
@DVar y(t)=sin(1)+exp(1)
@DVar z(t)
x1 = DependentVariable(:x,dependents = [t])
y1 = DependentVariable(:y, sin(1) + exp(1),dependents = [t])
z1 = DependentVariable(:z,dependents = [t])
@test x1 == x
@test y1 == y
@test z1 == z
@IVar begin
    t
    s = cos(2.5)
end
t1 = IndependentVariable(:t)
s1 = IndependentVariable(:s, cos(2.5))
@test t1 == t
@test s1 == s
@Deriv D''~t
D1 = Differential(t, 2)
@test D1 == D
@Const c=0 v=2
c1 = Constant(0)
v1 = Constant(2)
@test c1 == c
@test v1 == v
