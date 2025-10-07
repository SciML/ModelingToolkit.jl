using SparseArrays
using ModelingToolkit
import ModelingToolkit: bareiss!, find_pivot_col, bareiss_update!, swaprows!
import Base: swapcols!

function det_bareiss!(M)
    parity = 1
    _swaprows!(M, i, j) = (i != j && (parity = -parity); swaprows!(M, i, j))
    _swapcols!(M, i, j) = (i != j && (parity = -parity); swapcols!(M, i, j))
    # We only look at the last entry, so we don't care that the sub-diagonals are
    # garbage.
    zero!(M, i, j) = nothing
    rank = bareiss!(M, (_swapcols!, _swaprows!, bareiss_update!, zero!);
        find_pivot = find_pivot_col)
    return parity * M[end, end]
end

@testset "bareiss tests" begin
    # copy gives a dense matrix
    @testset "bareiss tests: $T" for T in (copy, sparse)
        # matrix determinant pairs
        for (M, d) in ((BigInt[9 1 8 0; 0 0 8 7; 7 6 8 3; 2 9 7 7], -1),
            (BigInt[1 big(2)^65+1; 3 4], 4 - 3 * (big(2)^65 + 1)))
            # test that the determinant was correctly computed
            @test det_bareiss!(T(M)) == d
        end
    end
end

# Bareiss regression test for relatively prime matrix
let M = reshape([2, 3, 5, 7, 11, 13, 17, 19, 23], 3, 3)
    morig = ModelingToolkit.SparseMatrixCLIL(M)
    m = copy(morig)
    # Primarily we care about that this does not throw, but it does also have full rank,
    # so let's test that
    @test ModelingToolkit.do_bareiss!(
        m, morig, [true for i in 1:3], [true for i in 1:3])[1:3] ==
          (3, 3, 3)
end
