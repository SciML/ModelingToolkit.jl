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

# Chunking is intentionally disabled when the loaded RuntimeGeneratedFunctions cannot compile helpers.
const CHUNKING_AVAILABLE = MTKB._rgf_supports_local_functions()
expect_chunked(ex) = CHUNKING_AVAILABLE ? count_helpers(ex) > 1 : count_helpers(ex) == 0

count_helpers(ex) = count(a -> Meta.isexpr(a, :call) && Meta.isexpr(a.args[1], :function), ex.args[2].args)

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
    du = zeros(n)
    run_expr(iip, du, u, p, 0.0)
    @test du == reference
    @test run_expr(oop, u, p, 0.0) == reference
end

@testset "automatic mode leaves small functions alone" begin
    oop, iip = build_rhs(sys, rhss; chunk_size = -1)
    @test count_helpers(iip) == 0
    @test oop == build_rhs(sys, rhss; chunk_size = 0)[1]
end

@testset "chunking is off by default and with chunk_size = 0" begin
    @test MTKB.BuildFunctionWrapperOptions().chunk_size == 0
    @test count_helpers(build_rhs(sys, rhss)[2]) == 0
    @test count_helpers(build_rhs(sys, rhss; chunk_size = 0)[2]) == 0
    @test (count_helpers(build_rhs(sys, rhss; chunk_size = 1)[2]) > 0) == CHUNKING_AVAILABLE
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
