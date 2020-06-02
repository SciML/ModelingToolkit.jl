using ModelingToolkit
using LinearAlgebra
using SparseArrays: sparse
using Test

@variables a,b,c,d

# test some matrix operations don't throw errors
X = [a b;c d]
det(X)
lu(X)
inv(X)

# test operations with sparse arrays and Operations
# note `isequal` instead of `==` because `==` would give another Operation

# test that we can create a sparse array of Operation
Oarray = zeros(Operation, 2,2)
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
