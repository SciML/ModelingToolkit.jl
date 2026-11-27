using ModelingToolkitBase, Test
using ModelingToolkitBase: unwrap, complete, shaped_parameter, MTKParameters
using Symbolics
using SymbolicUtils
using SciMLBase
using SparseArrays, LinearAlgebra
using SymbolicIndexingInterface: getp, setp, parameter_index, remake_buffer
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
    @test_throws "scalar type" shaped_parameter(:a, Float64, 4, 6)
    @test size(shaped_parameter(:A, Matrix{Float64}, 4, 6)) == (4, 6)
    for T in (Matrix{Float64}, SparseMatrixCSC{Float64, Int}, BlockOp)
        @test_throws "not tunable" shaped_parameter(:K2, T, 4, 6; tunable = true)
    end
    @variables u(t)[1:6] w(t)[1:4]
    @test size(K * u) == (4,)
    @test (w ~ K * u) isa Equation
    @test_throws ArgumentError shaped_parameter(:K, SparseMatrixCSC, 4, 6)
    @test_throws ArgumentError shaped_parameter(:K, BlockOp)
    @test_throws ArgumentError shaped_parameter(:K, BlockOp, 0, 3)
    @test_throws ArgumentError shaped_parameter(:K, BlockOp, 3, -1)
    @test_throws ArgumentError shaped_parameter(:K, BlockOp, 3, 1:0)
    @test size(shaped_parameter(:K, BlockOp, 1:4, 3)) == (4, 3)
    @test_throws "offset" shaped_parameter(:K, BlockOp, 2:5, 3)
    @test_throws "offset" shaped_parameter(:A, Matrix{Float64}, 2:5, 3)
    @test_throws "needs 2 dimensions" shaped_parameter(:A, Matrix{Float64}, 4)
    @test_throws "needs 2 dimensions" shaped_parameter(:A, SparseMatrixCSC{Float64, Int}, 4, 3, 2)
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
    replacement = BlockOp(sparse(1:4, 1:4, fill(2.0, 4), 4, 6))
    prob2 = remake(prob; p = [Kop => replacement])
    @test getp(prob2, Kop)(prob2).A == replacement.A
    @test getp(prob, Kop)(prob).A == good.A
    @test good.A != replacement.A

    ps = prob.p
    for bad in (badrows, badcols)
        @test_throws DimensionMismatch setp(csys, Kop)(ps, bad)
        @test_throws DimensionMismatch (ps[parameter_index(csys, Kop)] = bad)
        @test_throws DimensionMismatch remake_buffer(csys, ps, [parameter_index(csys, Kop)], [bad])
        # a non-concrete `vals` eltype takes the non-generated `remake_buffer` route
        @test_throws DimensionMismatch remake_buffer(csys, ps, [parameter_index(csys, Kop)], Any[bad])
    end
    setp(csys, Kop)(ps, good)
    @test getp(csys, Kop)(ps) === good
    setp(csys, Nop)(ps, NoSizeOp(2.0))
    @test getp(csys, Nop)(ps) == NoSizeOp(2.0)
end

@testset "bindings of shaped opaque parameters are rejected" begin
    @variables u(t)[1:6]
    Kop = shaped_parameter(:Kop, BlockOp, 4, 6)
    Nop = shaped_parameter(:Nop, NoSizeOp, 4, 6)
    Sop = shaped_parameter(:Sop, SparseMatrixCSC{Float64, Int}, 4, 6)
    Uop = shaped_parameter(:Uop, UniformScaling{Float64}, 4, 6)
    eqs = [zeros(6) ~ D(u) + u]
    cases = (
        (Kop, BlockOp(sprand(4, 6, 0.5))), (Nop, NoSizeOp(1.0)),
        (Sop, spzeros(4, 6)), (Uop, 1.0I), (Sop, missing),
    )
    for (p, val) in cases
        @test_throws r"bindings of shaped opaque parameters are not supported" System(
            eqs, t, [u], [p]; bindings = [p => val], name = :bad
        )
    end
end

@testset "shape-preserving bindings of parameters declared without shaped_parameter still work" begin
    A = sparse([1.0 2 3; 4 5 6])
    sym = SymbolicUtils.Sym{SymbolicUtils.SymReal}(
        :K; type = SparseMatrixCSC{Float64, Int},
        shape = SymbolicUtils.ShapeVecT([1:2, 1:3])
    )
    K = Symbolics.wrap(ModelingToolkitBase.toparam(sym))
    @variables u(t)[1:3] w(t)[1:2]
    sys = System([w ~ K * u], t, [u, w], [K]; bindings = [K => A], name = :legacy)
    @test sys isa System
    @test_throws "has size (3, 3), but K has size (2, 3)" System(
        [w ~ K * u], t, [u, w], [K]; bindings = [K => sparse(ones(3, 3))], name = :wrongsize
    )
end

@testset "shaped_parameter systems cannot be completed unsplit" begin
    @variables u(t)[1:6]
    Kop = shaped_parameter(:Kop, BlockOp, 4, 6)
    @named sys = System([zeros(6) ~ D(u) + u], t, [u], [Kop])
    @test_throws "split = false" complete(sys; split = false)
    @test_throws "split = false" complete(sys; split = false, flatten = false)
    @test complete(sys) isa System
    # the internal split = false snapshots of the parent and initialization system are not checked
    @variables x[1:6]
    @named init = System([zeros(6) ~ x], [x], [Kop]; is_initializesystem = true)
    @named isys = System([zeros(6) ~ D(u) + u], t, [u], [Kop]; initializesystem = init)
    @test complete(isys; split = true) isa System
    @test_throws "split = false" complete(isys; split = false)
end

struct VecFn end
(::VecFn)(x) = [x, 2x]
Base.size(::VecFn) = (5,)

@testset "array-returning callable parameters are not size-checked against the callable" begin
    @parameters (f::VecFn)(..)[1:2]
    @variables u(t)[1:2]
    @named sys = System([D(u) ~ f(t)], t, [u], [f])
    csys = complete(sys)
    ps = MTKParameters(csys, Dict(f => VecFn()))
    @test getp(csys, f)(ps) isa VecFn
    @test_nowarn System([D(u) ~ f(t)], t, [u], [f]; bindings = [f => VecFn()], name = :bound)
end

@testset "declared shape is validated for sparse matrix parameters" begin
    T = SparseMatrixCSC{Float64, Int}
    @variables u(t)[1:6]
    K = shaped_parameter(:K, T, 4, 6)
    @named sys = System([zeros(6) ~ D(u) + u], t, [u], [K])
    csys = complete(sys)
    u0 = Dict(u => zeros(6), D(u) => zeros(6))
    for bad in (spzeros(5, 6), spzeros(4, 5))
        @test_throws DimensionMismatch MTKParameters(csys, Dict(K => bad))
        @test_throws DimensionMismatch DAEProblem(
            csys, merge(u0, Dict(K => bad)), (0.0, 1.0); build_initializeprob = false
        )
    end
    prob = DAEProblem(
        csys, merge(u0, Dict(K => spzeros(4, 6))), (0.0, 1.0); build_initializeprob = false
    )
    @test_throws DimensionMismatch remake(prob; p = [K => spzeros(5, 6)])
    replacement = sparse(1:4, 1:4, fill(3.0, 4), 4, 6)
    prob2 = remake(prob; p = [K => replacement])
    @test getp(prob2, K)(prob2) == replacement
    @test getp(prob, K)(prob) == spzeros(4, 6)
    @test_throws DimensionMismatch remake_buffer(
        csys, prob.p, [parameter_index(csys, K)], Any[spzeros(5, 6)]
    )
end

@testset "parameters without a declared shape are not size-checked" begin
    @parameters legacy::BlockOp any_p::Any
    @variables u(t)[1:2]
    @named sys = System([zeros(2) ~ D(u) + u], t, [u], [legacy, any_p])
    csys = complete(sys)
    ps = MTKParameters(csys, Dict(legacy => BlockOp(sprand(2, 3, 0.5)), any_p => ones(2)))
    @test getp(csys, legacy)(ps) isa BlockOp
    @test getp(csys, any_p)(ps) == ones(2)
end

# OffsetArrays is a dependency of ModelingToolkitBase but not of the test environment.
const OffsetArrays = Base.loaded_modules[Base.PkgId(Base.UUID("6fe1bfb0-de20-5000-8ca7-80f57d26f881"), "OffsetArrays")]
const OffsetMat = OffsetArrays.OffsetMatrix{Float64, Matrix{Float64}}

# An operator whose size is that of a field that can be mutated after it is stored.
mutable struct MutableOp
    M::Matrix{Float64}
end
Base.size(op::MutableOp) = size(op.M)

@testset "shaped array parameters are validated by axes" begin
    A = shaped_parameter(:A, OffsetMat, 2, 2)
    @variables x(t)
    @named sys = System([D(x) ~ A[1, 1]], t, [x], [A])
    csys = complete(sys)
    raw = [1.0 2.0; 3.0 4.0]
    mat = OffsetArrays.OffsetArray(raw, 1:2, 1:2)
    offset = OffsetArrays.OffsetArray(raw, 0:1, 0:1)
    @test axes(offset) == (0:1, 0:1)
    du = zeros(1)

    ps = MTKParameters(csys, Dict(A => mat))
    prob = ODEProblem(csys, Dict(x => 0.0, A => mat), (0.0, 1.0))
    prob.f(du, [0.0], prob.p, 0.0)
    @test du[1] == 1.0

    idx = parameter_index(csys, A)
    @test_throws DimensionMismatch MTKParameters(csys, Dict(A => offset))
    @test_throws DimensionMismatch ODEProblem(csys, Dict(x => 0.0, A => offset), (0.0, 1.0))
    @test_throws DimensionMismatch setp(csys, A)(ps, offset)
    @test_throws DimensionMismatch (ps[idx] = offset)
    @test_throws DimensionMismatch remake(prob; p = [A => offset])
    @test_throws DimensionMismatch remake_buffer(csys, ps, [idx], [offset])
    @test_throws DimensionMismatch remake_buffer(csys, ps, [idx], Any[offset])
    @test_throws DimensionMismatch remake_buffer(csys, ps, [A], [offset])
    @test getp(csys, A)(ps) == mat
end

@testset "indexed updates are validated against the declared shape, not the stored value" begin
    Mop = shaped_parameter(:Mop, MutableOp, 2, 2)
    @variables u(t)[1:2]
    @named sys = System([zeros(2) ~ D(u) + u], t, [u], [Mop])
    csys = complete(sys)
    ps = MTKParameters(csys, Dict(Mop => MutableOp(ones(2, 2))))
    stored = getp(csys, Mop)(ps)
    stored.M = zeros(3, 2)
    idx = parameter_index(csys, Mop)
    bad = MutableOp(zeros(3, 2))
    good = MutableOp(fill(2.0, 2, 2))

    @test_throws DimensionMismatch setp(csys, Mop)(ps, bad)
    @test_throws DimensionMismatch (ps[idx] = bad)
    @test_throws DimensionMismatch remake_buffer(csys, ps, [idx], [bad])
    @test_throws DimensionMismatch remake_buffer(csys, ps, [idx], Any[bad])

    setp(csys, Mop)(ps, good)
    @test getp(csys, Mop)(ps) === good
    other = MutableOp(fill(3.0, 2, 2))
    ps[idx] = other
    @test getp(csys, Mop)(ps) === other
    for vals in ([good], Any[good])
        remade = remake_buffer(csys, ps, [idx], vals)
        @test getp(csys, Mop)(remade) === good
    end
end

@testset "bindings of dense shaped parameters are validated by axes" begin
    A = shaped_parameter(:A, Matrix{Float64}, 2, 2)
    @variables x(t)
    raw = [1.0 2.0; 3.0 4.0]
    offset = OffsetArrays.OffsetArray(raw, 0:1, 0:1)
    based = OffsetArrays.OffsetArray(raw, 1:2, 1:2)
    function rhs(binding)
        @named sys = System([D(x) ~ A[1, 1]], t, [x], [A]; bindings = [A => binding])
        prob = ODEProblem(complete(sys), Dict(x => 0.0), (0.0, 1.0))
        du = zeros(1)
        prob.f(du, [0.0], prob.p, 0.0)
        return du[1]
    end
    @test rhs(raw) == 1.0
    @test rhs(based) == 1.0
    @test_throws DimensionMismatch rhs(offset)
    @test_throws DimensionMismatch System(
        [D(x) ~ A[1, 1]], t, [x], [A]; bindings = [A => offset], name = :bad
    )
    # values supplied as initial conditions or defaults are validated when the buffer is built
    @named sys = System([D(x) ~ A[1, 1]], t, [x], [A]; initial_conditions = [A => offset])
    @test_throws DimensionMismatch ODEProblem(complete(sys), Dict(x => 0.0), (0.0, 1.0))
    @named plain = System([D(x) ~ A[1, 1]], t, [x], [A])
    @test_throws DimensionMismatch ODEProblem(
        complete(set_defaults(plain, [A => offset])), Dict(x => 0.0), (0.0, 1.0)
    )
end

@testset "completion enforces the binding restrictions on set_defaults bindings" begin
    @variables x(t)
    Sbig = shaped_parameter(:Sbig, SparseMatrixCSC{Float64, Int}, 2, 2)
    Ssmall = shaped_parameter(:Ssmall, SparseMatrixCSC{Float64, Int}, 1, 1)
    @named sys = System([D(x) ~ Sbig[2, 2]], t, [x], [Sbig, Ssmall])
    @test_throws r"bindings of shaped opaque parameters are not supported" complete(
        set_defaults(sys, [Sbig => Ssmall])
    )
    Abig = shaped_parameter(:Abig, Matrix{Float64}, 2, 2)
    Asmall = shaped_parameter(:Asmall, Matrix{Float64}, 1, 1)
    @named dense = System([D(x) ~ Abig[2, 2]], t, [x], [Abig, Asmall])
    @test_throws "has size (1, 1), but Abig has size (2, 2)" complete(
        set_defaults(dense, [Abig => Asmall])
    )
    @test complete(dense) isa System
end

@testset "dense shaped parameters cannot be bound to values of unknown shape" begin
    A = shaped_parameter(:A, Matrix{Float64}, 2, 2)
    @parameters B::Matrix{Float64}
    @variables x(t)
    for rhs in (B, 2 * B)
        @test_throws "unknown shape" System(
            [D(x) ~ A[2, 2]], t, [x], [A, B]; bindings = [A => rhs], name = :bad
        )
    end
    @named sys = System([D(x) ~ A[2, 2]], t, [x], [A, B])
    @test_throws "unknown shape" complete(set_defaults(sys, [A => B]))
    @test_throws "unknown shape" complete(set_defaults(sys, [A => 2 * B]))
    @test complete(sys) isa System
end
