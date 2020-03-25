using ModelingToolkit
using LinearAlgebra
@variables a,b,c,d
X = [a b;c d]
det(X)
lu(X)
inv(X)
