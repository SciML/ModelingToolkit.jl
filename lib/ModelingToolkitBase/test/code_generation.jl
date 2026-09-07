using ModelingToolkitBase, OrdinaryDiffEq, SymbolicIndexingInterface
using SciMLBase
using SciMLStructures
using SymbolicUtils: unwrap, iscall, operation
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using Test

@testset "`generate_custom_function`" begin
    @variables x(t) y(t)[1:3]
    @parameters p1 = 1.0 p2[1:3] = [1.0, 2.0, 3.0] p3::Int = 1 p4::Bool = false

    sys = complete(System(Equation[], t, [x; y], [p1, p2, p3, p4]; name = :sys))
    u0 = [1.0, 2.0, 3.0, 4.0]
    p = ModelingToolkitBase.MTKParameters(sys, [])

    fn1 = generate_custom_function(
        sys, x + y[1] + p1 + p2[1] + p3 * t; expression = Val(false)
    )
    @test fn1(u0, p, 0.0) == 5.0

    fn2 = generate_custom_function(
        sys, x + y[1] + p1 + p2[1] + p3 * t, [x], [p1, p2, p3]; expression = Val(false)
    )
    @test fn1(u0, p, 0.0) == 5.0

    fn3_oop,
        fn3_iip = generate_custom_function(
        sys, [x + y[2], y[3] + p2[2], p1 + p3, 3t]; expression = Val(false)
    )

    buffer = zeros(4)
    fn3_iip(buffer, u0, p, 1.0)
    @test buffer == [4.0, 6.0, 2.0, 3.0]
    @test fn3_oop(u0, p, 1.0) == [4.0, 6.0, 2.0, 3.0]

    fn4 = generate_custom_function(sys, ifelse(p4, p1, p2[2]); expression = Val(false))
    @test fn4(u0, p, 1.0) == 2.0
    fn5 = generate_custom_function(sys, ifelse(!p4, p1, p2[2]); expression = Val(false))
    @test fn5(u0, p, 1.0) == 1.0

    @variables x y[1:3]
    sys = complete(System(Equation[], [x; y], [p1, p2, p3, p4]; name = :sys))
    p = MTKParameters(sys, [])

    fn1 = generate_custom_function(sys, x + y[1] + p1 + p2[1] + p3; expression = Val(false))
    @test fn1(u0, p) == 6.0

    fn2 = generate_custom_function(
        sys, x + y[1] + p1 + p2[1] + p3, [x], [p1, p2, p3]; expression = Val(false)
    )
    @test fn1(u0, p) == 6.0

    fn3_oop,
        fn3_iip = generate_custom_function(
        sys, [x + y[2], y[3] + p2[2], p1 + p3]; expression = Val(false)
    )

    buffer = zeros(3)
    fn3_iip(buffer, u0, p)
    @test buffer == [4.0, 6.0, 2.0]
    @test fn3_oop(u0, p, 1.0) == [4.0, 6.0, 2.0]

    fn4 = generate_custom_function(sys, ifelse(p4, p1, p2[2]); expression = Val(false))
    @test fn4(u0, p, 1.0) == 2.0
    fn5 = generate_custom_function(sys, ifelse(!p4, p1, p2[2]); expression = Val(false))
    @test fn5(u0, p, 1.0) == 1.0
end

@testset "Non-standard array variables" begin
    @variables x(t)
    @parameters p[0:2] (f::Function)(..)
    @mtkcompile sys = System(D(x) ~ p[0] * x + p[1] * t + p[2] + f(p), t)
    prob = ODEProblem(sys, [x => 1.0, p => [1.0, 2.0, 3.0], f => sum], (0.0, 1.0))
    @test prob.ps[p] == [1.0, 2.0, 3.0]
    @test prob.ps[p[0]] == 1.0
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
end

@testset "scalarized array observed calling same function multiple times" begin
    @variables x(t) y(t)[1:2]
    @parameters foo(::Real)[1:2]
    val = Ref(0)
    function _tmp_fn2(x)
        val[] += 1
        return [x, 2x]
    end
    @mtkcompile sys = System([D(x) ~ y[1] + y[2], y ~ foo(x)], t)
    @test length(equations(sys)) == 1
    @test length(ModelingToolkitBase.observed(sys)) == 3
    prob = ODEProblem(sys, [x => 1.0, foo => _tmp_fn2], (0.0, 1.0))
    val[] = 0
    @test_nowarn prob.f(prob.u0, prob.p, 0.0)
    @test val[] == 1

    @testset "CSE in equations(sys)" begin
        val[] = 0
        @variables z(t)[1:2]
        @mtkcompile sys = System(
            [D(y) ~ foo(x), D(x) ~ sum(y), zeros(2) ~ foo(prod(z))], t
        )
        @test length(equations(sys)) == 5
        @test length(ModelingToolkitBase.observed(sys)) == 0
        prob = ODEProblem(
            sys, [y => ones(2), z => 2ones(2), x => 3.0, foo => _tmp_fn2], (0.0, 1.0)
        )
        val[] = 0
        @test_nowarn prob.f(prob.u0, prob.p, 0.0)
        @test val[] == 2
    end
end

@testset "Do not codegen redundant expressions" begin
    @variables v1(t) = 1
    @variables v2(t) [guess = 0]

    mutable struct Data
        count::Int
    end
    function update!(d::Data, t)
        d.count += 1 # Count the number of times the data gets updated.
    end
    function (d::Data)(t)
        update!(d, t)
        rand(1:10)
    end

    @parameters (d1::Data)(..) = Data(0)
    @parameters (d2::Data)(..) = Data(0)

    eqs = [
        D(v1) ~ d1(t),
    ]

    @named sys = System(eqs, t, [v1], [d1, d2]; observed = [v2 ~ d2(t)])
    sys = complete(sys)
    prob = ODEProblem(sys, [], (0.0, 1.0))
    # Manual solve because lack of tearing in MTKBase will cause `d2` to be called
    # when solving initialization.
    integ = init(prob, Tsit5())
    integ.ps[d2].count = 0
    solve!(integ)
    sol = integ.sol
    @test sol.ps[d2].count == 0
end

@testset "Derivatives dependent on observed" begin
    @variables x(t) y(t)
    @mtkcompile sys = System([D(x) ~ y, y ~ 2t + 1], t)
    @test length(equations(sys)) == 1
    v = ModelingToolkitBase.default_toterm(unwrap(D(x)))
    fn = ModelingToolkitBase.build_explicit_observed_function(sys, v)
    ps = MTKParameters(sys, nothing)
    @test fn([1.0], ps, 1.5) ≈ 4.0
end

@testset "Non-scalarized array observed without individual elements being unknowns/observables" begin
    @variables x(t)[1:3] y(t)
    @mtkcomplete sys = System([D(y) ~ 2y + sum(x)], t, [y], []; observed = [x ~ [y, y + 1, y + 2]])
    @test ModelingToolkitBase.observed_equations_used_by(sys, [x[1]]) == [1]
end

@testset "`observed_dependency_graph` result is cacheable" begin
    @variables x(t) y(t) z(t)
    @mtkcompile sys = System([D(x) ~ z, y ~ 2x + 1, z ~ 3y], t)
    obs = ModelingToolkitBase.observed(sys)
    graph = ModelingToolkitBase.observed_dependency_graph(sys, obs)
    @test graph isa fieldtype(ModelingToolkitBase.ObservedGraphCache, :graph)
    # this populates the observed graph cache, which requires the above type to match
    @test ModelingToolkitBase.observed_equations_used_by(sys, [equations(sys)[1].rhs]) ==
        [1, 2]
end

@testset "`BandedAMatrixWrapper` on an `Expr`" begin
    ex = ModelingToolkitBase.BandedAMatrixWrapper(:(f()), 3, (1, 1))
    @test ex == Expr(:call, ModelingToolkitBase.BandedAMatrixWrapper, :(f()), 3, (1, 1))
end

@testset "`calculate_paramjac`/`generate_paramjac`" begin
    @variables x(t) y(t)
    @parameters a b c[1:3] dd::Int e [tunable = false]
    @mtkcompile sys = System(
        [
            D(x) ~ a * x - b * x * y + c[1] * x + c[3] + dd + e,
            D(y) ~ -c[2] * y + x * y,
        ], t
    )
    opmap = [
        x => 1.0, y => 2.0, a => 10.0, b => 20.0, c => [31.0, 32.0, 33.0],
        dd => 2, e => 5.0,
    ]
    prob = ODEProblem(sys, opmap, (0.0, 1.0))
    p = prob.p
    tunable, repack, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), p)

    cols = ModelingToolkitBase.paramjac_parameters(sys)
    # the columns index the tunable buffer, so non-tunable, integer-valued and
    # `Initial(...)` parameters must not appear
    @test getp(sys, cols)(p) == tunable
    @test !any(isequal(unwrap(dd)), cols)
    @test !any(isequal(unwrap(e)), cols)
    @test !any(x -> iscall(x) && operation(x) isa Initial, cols)

    pjac = calculate_paramjac(sys)
    @test size(pjac) == (length(unknowns(sys)), length(tunable))

    fn = generate_paramjac(sys; expression = Val{false}, wrap_gfw = Val{true})
    u0 = prob.u0
    oop = fn(u0, p, 0.0)
    @test size(oop) == size(pjac)

    iip = zeros(size(pjac))
    fn(iip, u0, p, 0.0)
    @test iip == oop

    # central difference of the RHS over the tunable buffer, i.e. what a consumer
    # perturbing `canonicalize(Tunable(), p)[1]` would measure
    rhs = generate_rhs(sys; expression = Val{false}, wrap_gfw = Val{true})
    fd = zeros(size(pjac))
    for j in eachindex(tunable)
        h = sqrt(eps(Float64)) * max(abs(tunable[j]), 1.0)
        up = copy(tunable)
        dn = copy(tunable)
        up[j] += h
        dn[j] -= h
        fd[:, j] = (rhs(u0, repack(up), 0.0) - rhs(u0, repack(dn), 0.0)) ./ (2h)
    end
    @test oop ≈ fd atol = 1.0e-6

    sparse_fn = generate_paramjac(
        sys; expression = Val{false}, wrap_gfw = Val{true}, sparse = true
    )
    @test Array(sparse_fn(u0, p, 0.0)) == oop
    sparse_buf = similar(calculate_paramjac(sys; sparse = true), Float64)
    sparse_fn(sparse_buf, u0, p, 0.0)
    @test Array(sparse_buf) == oop

    # the wrapper must advertise both arities so `SciMLBase.isinplace` accepts it
    @test SciMLBase.numargs(fn) == (3, 4)
end

@testset "`paramjac` columns of an array parameter are column-major" begin
    # A problem cannot be built for a system with a matrix parameter yet, so pin the
    # column layout symbolically against the index cache instead. Every entry of `MM` has
    # to appear in the equations: an array parameter with an unused entry is missing from
    # the index cache entirely, which breaks `generate_rhs`/`generate_jacobian` too.
    @variables x(t) y(t)
    @parameters a MM[1:2, 1:2]
    @mtkcompile sys = System(
        [
            D(x) ~ a * x + MM[1, 1] * x + MM[1, 2] * y,
            D(y) ~ MM[2, 1] * x - MM[2, 2] * y,
        ], t
    )
    cols = ModelingToolkitBase.paramjac_parameters(sys)
    idx = parameter_index(sys, MM).idx
    # column `idx[i, j]` of `pJ` must be the derivative w.r.t. `MM[i, j]`
    @test all(cols[idx[i, j]] === unwrap(MM[i, j]) for i in 1:2, j in 1:2)
    # which for the index cache's column-major layout means consecutive columns
    @test all(
        isequal.(cols[vec(idx)], unwrap.([MM[1, 1], MM[2, 1], MM[1, 2], MM[2, 2]]))
    )
    @test size(calculate_paramjac(sys)) == (length(unknowns(sys)), length(cols))
end

@testset "`generate_paramjac` on a time-independent system" begin
    @variables xx yy
    @parameters aa bb
    @mtkcompile sys = System([0 ~ xx^2 - aa, 0 ~ yy - bb * xx])
    fn = generate_paramjac(sys; expression = Val{false}, wrap_gfw = Val{true})
    prob = NonlinearProblem(sys, [xx => 6.0, yy => 2.0, aa => 4.0, bb => 3.0])
    oop = fn(prob.u0, prob.p)
    iip = zeros(size(oop))
    fn(iip, prob.u0, prob.p)
    @test oop == iip
    # how far `mtkcompile` tears this system depends on whether the tearing
    # implementation is loaded, so derive the expected shape instead of fixing it
    @test size(oop) ==
        (length(unknowns(sys)), length(ModelingToolkitBase.paramjac_parameters(sys)))
end

@testset "`generate_paramjac` with no tunable parameters" begin
    @variables z(t)
    @parameters q [tunable = false]
    @mtkcompile sys = System([D(z) ~ -q * z], t)
    @test size(calculate_paramjac(sys)) == (length(unknowns(sys)), 0)
    @test_throws ArgumentError generate_paramjac(
        sys; expression = Val{false}, wrap_gfw = Val{true}
    )

    # a system with no parameters at all must reach the same error, not a `BoundsError`
    # from indexing the empty result of `reorder_parameters`
    @variables w(t)
    @mtkcompile nops = System([D(w) ~ -w], t)
    @test isempty(ModelingToolkitBase.paramjac_parameters(nops))
    @test_throws ArgumentError generate_paramjac(
        nops; expression = Val{false}, wrap_gfw = Val{true}
    )
end

@testset "`calculate_paramjac` on a non-split system" begin
    @variables x(t) y(t)
    @parameters a b
    sys = mtkcompile(
        System([D(x) ~ a * x - b * x * y, D(y) ~ -a * y + b * x * y], t; name = :sys);
        split = false
    )
    prob = ODEProblem(sys, [x => 1.0, y => 2.0, a => 1.5, b => 0.5], (0.0, 1.0))
    tunable, repack, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), prob.p)
    # without an index cache `p` is a flat vector of every parameter, so every one of
    # them gets a column
    @test size(calculate_paramjac(sys), 2) ==
        length(parameters(sys; initial_parameters = true)) == length(tunable)

    fn = generate_paramjac(sys; expression = Val{false}, wrap_gfw = Val{true})
    rhs = generate_rhs(sys; expression = Val{false}, wrap_gfw = Val{true})
    u0 = prob.u0
    oop = fn(u0, prob.p, 0.0)
    iip = zeros(size(oop))
    fn(iip, u0, prob.p, 0.0)
    @test iip == oop

    fd = zeros(size(oop))
    for j in eachindex(tunable)
        h = sqrt(eps(Float64)) * max(abs(tunable[j]), 1.0)
        up = copy(tunable)
        dn = copy(tunable)
        up[j] += h
        dn[j] -= h
        fd[:, j] = (rhs(u0, repack(up), 0.0) - rhs(u0, repack(dn), 0.0)) ./ (2h)
    end
    @test oop ≈ fd atol = 1.0e-6
end

@testset "`ODEFunction` `paramjac`" begin
    @variables x(t) y(t)
    @parameters a b
    @mtkcompile sys = System([D(x) ~ a * x - b * x * y, D(y) ~ -a * y + b * x * y], t)
    opmap = [x => 1.0, y => 2.0, a => 1.5, b => 0.5]

    @test ODEFunction(sys).paramjac === nothing
    for T in (true, false)
        fn = ODEFunction{T}(sys; paramjac = true)
        @test fn.paramjac !== nothing
    end

    ref = ODEProblem(sys, opmap, (0.0, 1.0))
    u0, p = ref.u0, ref.p
    expected = generate_paramjac(sys; expression = Val{false}, wrap_gfw = Val{true})(
        u0, p, 0.0
    )
    @test ODEFunction(sys; paramjac = true).paramjac(u0, p, 0.0) == expected

    prob = ODEProblem(sys, opmap, (0.0, 1.0); paramjac = true)
    @test prob.f.paramjac !== nothing
    @test !haskey(prob.kwargs, :paramjac)
end

@testset "`calculate_paramjac` respects `SymbolicADDisallowed`" begin
    @variables x(t)
    @parameters a
    @mtkcompile sys = System([D(x) ~ a * x], t)
    sys = SymbolicUtils.setmetadata(
        sys, ModelingToolkitBase.SymbolicADDisallowed, "test reason"
    )
    @test_throws ArgumentError calculate_paramjac(sys)
end
