using ModelingToolkitBase, Test
using ModelingToolkitBase: unwrap, complete, shaped_parameter
using Symbolics
using SymbolicUtils
using SciMLBase
using SparseArrays, LinearAlgebra
using SymbolicIndexingInterface: getp
using OrdinaryDiffEqBDF: DFBDF
using DiffEqBase: BrownFullBasicInit

struct BlockOp
    A::SparseMatrixCSC{Float64, Int}
end
Base.size(op::BlockOp) = size(op.A)
apply_op(op::BlockOp, x::AbstractVector) = op.A * x
@register_array_symbolic apply_op(op::BlockOp, x::AbstractVector) begin
    size = (size(op, 1),)
    eltype = Float64
    ndims = 1
end

@independent_variables t
D = Differential(t)

@testset "declaration" begin
    K = shaped_parameter(:K, SparseMatrixCSC{Float64, Int}, 4, 6)
    @test size(K) == (4, 6)
    @test SymbolicUtils.symtype(unwrap(K)) == SparseMatrixCSC{Float64, Int}
    @test ModelingToolkitBase.isparameter(K)
    @test !ModelingToolkitBase.istunable(K)
    @test ModelingToolkitBase.istunable(shaped_parameter(:K2, SparseMatrixCSC{Float64, Int}, 4, 6; tunable = true))
    @variables u(t)[1:6] w(t)[1:4]
    @test size(K * u) == (4,)
    @test (w ~ K * u) isa Equation
    @test_throws ArgumentError shaped_parameter(:K, SparseMatrixCSC, 4, 6)
    @test_throws ArgumentError shaped_parameter(:K, BlockOp)
end

@testset "rectangular registered operator" begin
    Kop = shaped_parameter(:Kop, BlockOp, 4, 6)
    @variables u(t)[1:6] w(t)[1:4]
    @test size(apply_op(Kop, u)) == (4,)
    @test (w ~ apply_op(Kop, u)) isa Equation
end

@testset "two fields with unequal block sizes build a DAEProblem" begin
    nu, nw = 7, 3
    Kuu = sparse(2.0I, nu, nu) + 0.1 * sprand(nu, nu, 0.3)
    Kuw = sprand(nu, nw, 0.6)
    Kwu = sprand(nw, nu, 0.6)
    Kww = sparse(3.0I, nw, nw)
    Muu = sparse(1.0I, nu, nu)
    Mww = sparse(1.0I, nw, nw)
    blocks = (Muu = Muu, Mww = Mww, Kuu = Kuu, Kuw = Kuw, Kwu = Kwu, Kww = Kww)
    ps = Dict(
        k => shaped_parameter(k, BlockOp, size(A)...) for (k, A) in pairs(blocks)
    )
    @variables u(t)[1:nu] w(t)[1:nw]
    eqs = [
        zeros(nu) ~ apply_op(ps[:Muu], D(u)) + apply_op(ps[:Kuu], u) + apply_op(ps[:Kuw], w),
        zeros(nw) ~ apply_op(ps[:Mww], D(w)) + apply_op(ps[:Kwu], u) + apply_op(ps[:Kww], w),
    ]
    @named sys = System(eqs, t, [u, w], collect(values(ps)))
    csys = complete(sys)

    u0, w0 = rand(nu), rand(nw)
    du0 = -(Kuu * u0 + Kuw * w0)
    dw0 = -(Kwu * u0 + Kww * w0)
    op = Dict{Any, Any}(u => u0, w => w0, D(u) => du0, D(w) => dw0)
    for (k, A) in pairs(blocks)
        op[ps[k]] = BlockOp(A)
    end
    prob = DAEProblem(
        csys, op, (0.0, 1.0); initializealg = BrownFullBasicInit(),
        build_initializeprob = false
    )

    Kfull = [Kuu Kuw; Kwu Kww]
    hand!(r, dx, x) = (r .= dx; mul!(r, Kfull, x, 1, 1); r)
    res = zeros(nu + nw)
    prob.f(res, [du0; dw0], [u0; w0], prob.p, 0.0)
    @test norm(res) < 1.0e-12
    dx, x = rand(nu + nw), rand(nu + nw)
    prob.f(res, dx, x, prob.p, 0.0)
    @test res ≈ hand!(zeros(nu + nw), dx, x)

    sol = solve(prob, DFBDF(); reltol = 1.0e-10, abstol = 1.0e-10)
    @test SciMLBase.successful_retcode(sol)
    @test length(sol[w][end]) == nw
    @test sol.u[end] ≈ exp(-Matrix(Kfull)) * [u0; w0] rtol = 1.0e-5
end

@testset "declared shape is validated against supplied values" begin
    struct NoSizeOp
        c::Float64
    end
    @variables u(t)[1:6]
    Kop = shaped_parameter(:Kop, BlockOp, 4, 6)
    Nop = shaped_parameter(:Nop, NoSizeOp, 4, 6)
    @named sys = System([zeros(6) ~ D(u) + u], t, [u], [Kop, Nop])
    csys = complete(sys)
    good = BlockOp(sprand(4, 6, 0.5))
    badrows = BlockOp(sprand(5, 6, 0.5))
    badcols = BlockOp(sprand(4, 5, 0.5))
    u0 = Dict(u => zeros(6), D(u) => zeros(6))
    for bad in (badrows, badcols)
        @test_throws DimensionMismatch DAEProblem(
            csys, merge(u0, Dict(Kop => bad, Nop => NoSizeOp(1.0))), (0.0, 1.0);
            build_initializeprob = false
        )
    end
    prob = DAEProblem(
        csys, merge(u0, Dict(Kop => good, Nop => NoSizeOp(1.0))), (0.0, 1.0);
        build_initializeprob = false
    )
    @test getp(prob, Nop)(prob) == NoSizeOp(1.0)
    @test_throws DimensionMismatch remake(prob; p = [Kop => badrows])
    @test_throws DimensionMismatch remake(prob; p = [Kop => badcols])
    prob2 = remake(prob; p = [Kop => BlockOp(sprand(4, 6, 0.2))])
    @test getp(prob2, Kop)(prob2) isa BlockOp
end
