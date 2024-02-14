using ModelingToolkit
using LinearAlgebra
using Test

A = [0 1 1 2 2 1 1 2 1 2
     0 1 -1 -3 -2 2 1 -5 0 -5
     0 1 2 2 1 1 2 1 1 2
     0 1 1 1 2 1 1 2 2 1
     0 2 1 2 2 2 2 1 1 1
     0 1 1 1 2 2 1 1 2 1
     0 2 1 2 2 1 2 1 1 2
     0 1 7 17 14 2 1 19 4 23
     0 1 -1 -3 -2 1 1 -4 0 -5
     0 1 1 2 2 1 1 2 2 2]
N = ModelingToolkit.nullspace(A)
@test size(N, 2) == 3
@test rank(N) == 3
@test iszero(A * N)

A = [0 1 2 0 1 0;
     0 0 0 0 0 1;
     0 0 0 0 0 1;
     1 0 1 2 0 1;
     0 0 0 2 1 0]
col_order = Int[]
N = ModelingToolkit.nullspace(A; col_order)
colspan = A[:, col_order[1:4]]   # rank is 4
@test iszero(ModelingToolkit.nullspace(colspan))
@test !iszero(ModelingToolkit.nullspace(A[:, col_order[1:5]]))
@test !iszero(ModelingToolkit.nullspace(A[:, [col_order[1:4]..., col_order[6]]]))
