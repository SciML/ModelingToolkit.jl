using ModelingToolkitBase, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D
import ModelingToolkitBase as MTKB

# A ring of coupled equations: neighbouring outputs share a little work, distant ones none.
function ring_system(n; name = :ring)
    @variables x(t)[1:n]
    @parameters k
    eqs = [
        D(x[i]) ~ k * (x[mod1(i + 1, n)] - 2x[i] + x[mod1(i - 1, n)])^2 / (1 + x[i]^2)
            for i in 1:n
    ]
    return mtkcompile(System(eqs, t; name))
end

function build_rhs(sys, rhss; kwargs...)
    p = MTKB.reorder_parameters(sys)
    args = collect(Any, (MTKB.flat_unknowns(sys), p..., t))
    opts = MTKB.BuildFunctionWrapperOptions(;
        p_start = 2, u_arg = 1, n_param_buffers = length(p), kwargs...
    )
    return MTKB.build_function_wrapper(sys, rhss, args, opts)
end

run_expr(ex, args...) = Base.invokelatest(eval(ex), args...)

# The helper calls of a combined function: the statements of an in-place one, the arguments of
# the `vcat` of an out-of-place one.
is_helper_call(a) = Meta.isexpr(a, :call) && Meta.isexpr(a.args[1], :function)
function count_helpers(ex)
    calls = filter(a -> a isa Expr, ex.args[2].args)
    vcats = filter(a -> Meta.isexpr(a, :call) && a.args[1] === :vcat, calls)
    return count(is_helper_call, calls) + sum(c -> count(is_helper_call, c.args), vcats; init = 0)
end
expect_chunked(ex) = count_helpers(ex) > 1

n = 200
sys = ring_system(n)
rhss = [eq.rhs for eq in equations(sys)]
p = MTKB.MTKParameters(sys, [only(parameters(sys)) => 0.7])
u = rand(n)
reference = zeros(n)
run_expr(build_rhs(sys, rhss; chunk_size = 0)[2], reference, u, p, 0.0)

@testset "explicit chunk_size splits the in-place function and matches" begin
    oop, iip = build_rhs(sys, rhss; chunk_size = 300)
    @test expect_chunked(iip)
    @test count_helpers(oop) == 0
    du = zeros(n)
    run_expr(iip, du, u, p, 0.0)
    @test du == reference
    @test run_expr(oop, u, p, 0.0) == reference
end

@testset "automatic mode leaves small functions alone" begin
    oop, iip = build_rhs(sys, rhss; chunk_size = -1)
    @test count_helpers(iip) == 0
    @test oop == build_rhs(sys, rhss; chunk_size = 0)[1]
    @test count_helpers(oop) == 0
end

@testset "chunking is off by default and with chunk_size = 0" begin
    @test MTKB.BuildFunctionWrapperOptions().chunk_size == 0
    @test count_helpers(build_rhs(sys, rhss)[2]) == 0
    @test count_helpers(build_rhs(sys, rhss; chunk_size = 0)[2]) == 0
    @test count_helpers(build_rhs(sys, rhss; chunk_size = 1)[2]) > 0
end

@testset "functions that cannot be split independently are left whole" begin
    @test count_helpers(build_rhs(sys, rhss; chunk_size = 300, create_bindings = true)[2]) == 0
    @test count_helpers(build_rhs(sys, rhss; chunk_size = 300, extra_assignments = [MTKB.Assignment(:unused, 1)])[2]) == 0
    @test count_helpers(build_rhs(sys, rhss; chunk_size = 300, output_type = Tuple)[2]) == 0
    # scalar (non-vector) results
    @test count_helpers(build_rhs(sys, first(rhss); chunk_size = 1)[2]) == 0
end

@testset "automatic mode declines when helpers would recompute shared work" begin
    # every output reads the same large subexpression, so the helpers would each redo it
    @variables y(t)[1:20]
    @parameters kk
    shared = sum(y[i]^2 * sin(y[i] + kk) for i in 1:20)
    eqs = [D(y[i]) ~ shared * y[i] for i in 1:20]
    sys2 = mtkcompile(System(eqs, t; name = :shared))
    rhss2 = [eq.rhs for eq in equations(sys2)]
    @test expect_chunked(build_rhs(sys2, rhss2; chunk_size = 20)[2])
    ir = MTKB.SU.subset_ir(MTKB.get_irstructure(sys2), rhss2)
    @test MTKB._chunk_ranges(ir, rhss2, -1) === nothing
end

@testset "automatic mode splits a large ODE right-hand side" begin
    nbig = 3000
    bigsys = ring_system(nbig; name = :bigring)
    ubig = rand(nbig)
    pbig = MTKB.MTKParameters(bigsys, [only(parameters(bigsys)) => 0.7])
    ir = MTKB.SU.subset_ir(MTKB.get_irstructure(bigsys), [eq.rhs for eq in equations(bigsys)])
    @test length(ir) > MTKB.DEFAULT_CHUNK_THRESHOLD
    big_rhss = [eq.rhs for eq in equations(bigsys)]
    # generic build_function_wrapper callers are not chunked by default
    @test count_helpers(build_rhs(bigsys, big_rhss)[2]) == 0
    expr_oop, expr_iip = MTKB.generate_rhs(bigsys, MTKB.GeneratedFunctionOptions(; expression = Val{true}))
    @test expect_chunked(expr_iip)
    expected = zeros(nbig)
    run_expr(build_rhs(bigsys, big_rhss)[2], expected, ubig, pbig, 0.0)
    du = zeros(nbig)
    run_expr(expr_iip, du, ubig, pbig, 0.0)
    @test du == expected
    prob = ODEProblem(bigsys, [bigsys.x => ubig, only(parameters(bigsys)) => 0.7], (0.0, 1.0))
    du_prob = zeros(nbig)
    prob.f(du_prob, ubig, prob.p, 0.0)
    @test du_prob == expected
end

@testset "aliased output still works without chunking; chunked with separate buffers is bitwise equal" begin
    @variables y(t)[1:4]
    ring = mtkcompile(System([D(y[i]) ~ y[mod1(i + 1, 4)]^2 for i in 1:4], t; name = :small))
    ring_rhss = [eq.rhs for eq in equations(ring)]
    ring_p = MTKB.MTKParameters(ring, [])
    for kwargs in ((;), (; chunk_size = 0))
        iip = build_rhs(ring, ring_rhss; kwargs...)[2]
        buffer = [1.0, 2.0, 3.0, 4.0]
        @test run_expr(iip, buffer, buffer, ring_p, 0.0) === buffer
        @test buffer == [4.0, 9.0, 16.0, 1.0]
    end
    chunked = build_rhs(ring, ring_rhss; chunk_size = 1)[2]
    @test expect_chunked(chunked)
    out = zeros(4)
    @test run_expr(chunked, out, [1.0, 2.0, 3.0, 4.0], ring_p, 0.0) === out
    @test out == [4.0, 9.0, 16.0, 1.0]
end

@testset "chunk sizes are the statements reachable from each slice" begin
    full_ir = MTKB.get_irstructure(sys)
    reachable(r) = length(MTKB.SU.subset_ir(full_ir, rhss[r]))
    total = reachable(1:n)
    for target in (150, 300, 1000)
        ranges = MTKB._chunk_ranges(full_ir, rhss, target)
        @test ranges !== nothing
        @test reduce(vcat, collect.(ranges)) == 1:n
        # a slice is closed as soon as it reaches the target, and not before
        @test all(r -> reachable(r) >= target, ranges[1:(end - 1)])
        @test all(r -> reachable(first(r):(last(r) - 1)) < target, ranges[1:(end - 1)])
    end
    # the unrelated statements of a shared IR do not affect the sizes
    subset_ir = MTKB.SU.subset_ir(full_ir, rhss)
    @test MTKB._chunk_ranges(subset_ir, rhss, 300) == MTKB._chunk_ranges(full_ir, rhss, 300)
    @test MTKB._chunk_ranges(full_ir, rhss, total) === nothing
    @test MTKB._chunk_ranges(full_ir, rhss, 0) === nothing
end

@testset "non-Array output types are not chunked" begin
    # a minimal labelled-vector type whose `vcat` would not reproduce the requested output
    struct LabelledVec <: AbstractVector{Float64}
        data::Vector{Float64}
        labels::Vector{Symbol}
    end
    Base.size(v::LabelledVec) = size(v.data)
    Base.getindex(v::LabelledVec, i::Int) = v.data[i]
    Base.setindex!(v::LabelledVec, x, i::Int) = (v.data[i] = x)
    Base.similar(v::LabelledVec, ::Type{Float64}, dims::Dims{1}) =
        LabelledVec(zeros(dims), v.labels[1:min(end, dims[1])])
    for output_type in (LabelledVec, Vector{Float64}, Array)
        kwargs = (; chunk_size = 1, output_type)
        expected = output_type === LabelledVec ? 0 : 1
        @test (count_helpers(build_rhs(sys, rhss; kwargs...)[2]) > 0) == (expected > 0)
    end
    @test build_rhs(sys, rhss; chunk_size = 1, output_type = LabelledVec) ==
        build_rhs(sys, rhss; chunk_size = 0, output_type = LabelledVec)
end

@testset "unspecified output type keeps the runtime container of the out-of-place function" begin
    # `similar` keeps the wrapper, `vcat` of wrappers does not
    struct LabelledU <: AbstractVector{Float64}
        data::Vector{Float64}
    end
    Base.size(v::LabelledU) = size(v.data)
    Base.getindex(v::LabelledU, i::Int) = v.data[i]
    Base.setindex!(v::LabelledU, x, i::Int) = (v.data[i] = x)
    Base.similar(::LabelledU, ::Type{Float64}, dims::Dims{1}) = LabelledU(zeros(dims))
    MTKB.SU.Code.create_array(::Type{LabelledU}, ::Nothing, ::Val{1}, ::Val, vals...) =
        LabelledU(collect(Float64, vals))
    Base.vcat(a::LabelledU, rest::LabelledU...) = vcat(a.data, (r.data for r in rest)...)

    oop, iip = build_rhs(sys, rhss; chunk_size = 1)
    oop0, iip0 = build_rhs(sys, rhss; chunk_size = 0)
    @test count_helpers(iip) > 1
    @test count_helpers(oop) == 0
    @test oop == oop0
    lu = LabelledU(u)
    result = run_expr(oop, lu, p, 0.0)
    @test result isa LabelledU
    @test typeof(result) === typeof(run_expr(oop0, lu, p, 0.0))
    @test result == reference
    du = zeros(n)
    run_expr(iip, du, lu, p, 0.0)
    @test du == reference
    # an explicit plain Array type is still `vcat`ed
    @test count_helpers(build_rhs(sys, rhss; chunk_size = 1, output_type = Array)[1]) > 1
end
