using ModelingToolkit
using LinearAlgebra
using SparseArrays: sparse
using Test

@variables a,b,c,d,e,f,g,h,i

# test hashing
aa = a; # old a

@variables a

@test isequal(a, aa)
@test hash(a) == hash(aa)

@test isequal(get_variables(a+aa+1), [a])

@test hash(a+b ~ c+d) == hash(a+b ~ c+d)

# test some matrix operations don't throw errors
X = [0 b c; d e f; g h i]
@test iszero(simplify(det(X) - ((d * ((b * i) - (c * h))) + (g * ((b * f) - (c * e))))))
F = lu(X)
@test F.p == [2, 1, 3]
R = simplify.(F.L * F.U - X[F.p, :], polynorm=true)
@test iszero(R)
@test simplify.(F \ X) == I
@test ModelingToolkit._solve(X, X, true) == I
inv(X)
qr(X)

X2 = [0 b c; 0 0 0; 0 h 0]
@test_throws SingularException lu(X2)
F2 = lu(X2, check=false)
@test F2.info == 1

# test operations with sparse arrays and Operations
# note `isequal` instead of `==` because `==` would give another Operation

# test that we can create a sparse array of Operation
Oarray = zeros(Num, 2,2)
Oarray[2,2] = a
@test isequal(sparse(Oarray), sparse([2], [2], [a]))

# test Operation * sparse
@test isequal(a * sparse([2], [2], [1]), sparse([2], [2], [a * 1]))

# test sparse{Operation} + sparse{Operation}
A = sparse([2], [2], [a])
B = sparse([2], [2], [b])
@test isequal(A + B, sparse([2], [2], [a+b]))

# test sparse{Operation} * sparse{Operation}
C = sparse([1, 2], [2, 1], [c, c])
D = sparse([1, 2], [2, 1], [d, d])

@test isequal(C * D, sparse([1,2], [1,2], [c * d, c * d]))

@parameters t σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)
Dx = Differential(x)
Dy = Differential(y)
Dz = Differential(z)

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
J = Num[Dx(eqs[1].rhs) Dy(eqs[1].rhs) Dz(eqs[1].rhs)
 Dx(eqs[2].rhs) Dy(eqs[2].rhs) Dz(eqs[2].rhs)
 Dx(eqs[3].rhs) Dy(eqs[3].rhs) Dz(eqs[3].rhs)]

J = expand_derivatives.(J)
using LinearAlgebra
luJ = lu(J,Val(false))

using ModelingToolkit
@variables M[1:2,1:2]
inv(M)

@variables b[1:2]
M = [1 0; 0 2]
M \ b
M \ reshape(b,2,1)
M = [1 1; 0 2]
M \ reshape(b,2,1)


M = [1 a; 0 2]
M \ b
M \ [1, 2]

# test det
@variables X[1:4,1:4]
d1 = det(X, laplace=true)
d2 = det(X, laplace=false)
_det1 = eval(build_function(d1, X))
_det2 = eval(build_function(d2, X))
A = [1 1 1 1
     1 0 1 1
     1 1 0 1
     1 1 1 0]
@test _det1(A) == -1
@test _det2(A) == -1

@variables X[1:3,1:3]
d1 = det(X, laplace=true)
d2 = det(X, laplace=false)
_det1 = eval(build_function(d1, X))
_det2 = eval(build_function(d2, X))
A = [1 1 1
     1 0 1
     1 1 1]
@test _det1(A) == 0
@test _det2(A) == 0

@variables a b c d
z1 = a + b * im
z2 = c + d * im
@test z1 * 2 - Complex(2a, 2b) == 0
@test isequal(2z1, Complex(2a, 2b))
@test isequal(z1 / z1, 1)
@test isequal(z1 / z2, Complex((a*c + b*d)/(c^2 + d^2), (b*c - a*d)/(c^2 + d^2)))
@test isequal(1 / z2, Complex(c/(c^2 + d^2), -d/(c^2 + d^2)))
@test isequal(z1 * z2, Complex(a*c - b*d, a*d + b*c))
@test isequal(z1 - z2, Complex(a - c, b - d))
@test isequal(z1 + z2, Complex(a + c, b + d))
@test isequal(z1 + 2, Complex(a + 2, b))
@test isequal(2 + z1, Complex(2 + a, b))
@test isequal(z1 - 2, Complex(a - 2, b))
@test isequal(2 - z1, Complex(2 - a, -b))
