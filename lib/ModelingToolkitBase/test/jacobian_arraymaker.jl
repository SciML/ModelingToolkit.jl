using ModelingToolkitBase, Test, SparseArrays, LinearAlgebra
using ModelingToolkitBase: t_nounits as t, D_nounits as D
import ModelingToolkitBase as MTKB

# A large in-place Jacobian is generated from an `ArrayMaker` instead of a matrix of
# expressions, if Symbolics has the option that makes that worth it.
function ring_system(n; name = :ring)
    @variables x(t)[1:n]
    @parameters k
    eqs = [
        D(x[i]) ~ k * (x[mod1(i + 1, n)] - 2x[i] + x[mod1(i - 1, n)])^2 / (1 + x[i]^2) +
            2.5 * x[mod1(i + 3, n)]
            for i in 1:n
    ]
    return mtkcompile(System(eqs, t; name))
end

run_expr(ex, args...) = Base.invokelatest(eval(ex), args...)

n = cld(MTKB.JACOBIAN_ARRAYMAKER_MIN, 4) + 10
sys = ring_system(n)
p = MTKB.MTKParameters(sys, [only(parameters(sys)) => 0.7])
u = rand(n)
jac(; kwargs...) = MTKB.generate_jacobian(sys; expression = Val{true}, chunk_size = 0, kwargs...)
# The entries of an `ArrayMaker` that are not constants are written in batches.
uses_arraymaker(ex) = occursin("batched_setindex!", string(ex))

if hasfield(MTKB.Symbolics.CodegenFunctionOptions, :eager_scalar_writes)
    @testset "sparse" begin
        oop, iip = jac(; sparse = true)
        @test uses_arraymaker(iip)
        Jref = run_expr(oop, u, p, 0.0)
        @test nnz(Jref) >= MTKB.JACOBIAN_ARRAYMAKER_MIN
        J = copy(Jref)
        fill!(nonzeros(J), NaN)
        @test run_expr(iip, J, u, p, 0.0) === J
        @test J == Jref

        # The pattern of the output is checked before anything is written.
        checked = jac(; sparse = true, checkbounds = true)[2]
        @test uses_arraymaker(checked)
        @test run_expr(checked, copy(Jref), u, p, 0.0) == Jref
        @test_throws AssertionError run_expr(checked, sparse(ones(n, n)), u, p, 0.0)

        # The writes are not checked, so the number of stored values of the output is.
        short = sparse(Diagonal(ones(n)))
        @test_throws DimensionMismatch run_expr(iip, short, u, p, 0.0)
        @test nonzeros(short) == ones(n)

        # Options that only the matrix of expressions is generated for.
        @test !uses_arraymaker(jac(; sparse = true, optimize = Any[])[2])
        @test !uses_arraymaker(jac(; sparse = true, skipzeros = true)[2])
    end

    @testset "dense" begin
        # A Jacobian with too few entries that are not zero.
        m = cld(MTKB.JACOBIAN_ARRAYMAKER_MIN, 8)
        small = ring_system(m; name = :small)
        @test !uses_arraymaker(MTKB.generate_jacobian(small; expression = Val{true})[2])

        oop, iip = jac()
        @test uses_arraymaker(iip)
        Jref = run_expr(oop, u, p, 0.0)
        J = fill(NaN, n, n)
        @test run_expr(iip, J, u, p, 0.0) === J
        @test J == Jref
        # Any other output is written by the function generated from the matrix of
        # expressions: this one does not accept the zeros written first on its diagonal,
        # and a view may share memory between its entries.
        matrix_iip = jac(; skipzeros = false, optimize = Any[])[2]
        # An element type that zeros of another type are not converted to is another case.
        makers = (
            () -> view(fill(NaN, n + 1, n), 1:n, :), () -> fill(NaN, n + 1, n),
            () -> Matrix{Union{Missing, Float64}}(missing, n, n),
        )
        for make in makers
            Jother = make()
            Jexpected = make()
            run_expr(matrix_iip, Jexpected, u, p, 0.0)
            @test run_expr(iip, Jother, u, p, 0.0) === Jother
            @test isequal(Jother, Jexpected)
        end
        lazy = MTKB.LazyGeneratedFunction(matrix_iip, MTKB.GeneratedFunctionOptions())
        @test lazy == MTKB.LazyGeneratedFunction(copy(matrix_iip), MTKB.GeneratedFunctionOptions())
        @test hash(lazy) == hash(MTKB.LazyGeneratedFunction(copy(matrix_iip), MTKB.GeneratedFunctionOptions()))
        Jlazy = fill(NaN, n, n)
        lazy(Jlazy, u, p, 0.0)
        @test Jlazy == Jref
        # Bounds checks of the writes are only generated for the matrix of expressions.
        @test !uses_arraymaker(jac(; checkbounds = true)[2])
    end

    @testset "through ODEProblem" begin
        for sparse in (false, true)
            prob = ODEProblem(sys, [unknowns(sys) .=> u; only(parameters(sys)) => 0.7], (0.0, 1.0); jac = true, sparse)
            Jref = run_expr(jac(; sparse)[1], prob.u0, prob.p, 0.0)
            J = sparse ? copy(prob.f.jac_prototype) : zeros(n, n)
            prob.f.jac(J, prob.u0, prob.p, 0.0)
            @test J == Jref
        end
    end
else
    @testset "not used without the Symbolics option" begin
        @test !uses_arraymaker(jac(; sparse = true)[2])
        @test !uses_arraymaker(jac()[2])
    end
end
