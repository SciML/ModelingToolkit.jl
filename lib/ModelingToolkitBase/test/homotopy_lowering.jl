using ModelingToolkitBase
using ModelingToolkitBase: homotopy, has_homotopy, has_any_homotopy,
    _rewrite_with_lambda, lower_homotopy,
    t_nounits as t, D_nounits as D
using Symbolics
using Test

@testset "homotopy stays opaque through complete(): no λ parameter is injected" begin
    @variables x(t) y(t)
    @named sys = System([D(x) ~ homotopy(-x, -x + 1), y ~ homotopy(x^2, x)], t)
    csys = complete(sys)
    @test !any(p -> occursin("homotopy", string(p)), parameters(csys))
    @test has_any_homotopy(csys)                       # nodes survive untouched
end

@testset "runtime numeric fallback is the actual branch" begin
    @test ModelingToolkitBase.homotopy(2.0, 5.0) == 2.0
end

@testset "registered derivative is nodewise (Modelica-faithful)" begin
    @variables a b
    # d/da homotopy(a^2, b) = homotopy(1,0) * 2a  — the homotopy(1,0) factor is an
    # opaque node (NOT folded to 1) so the continuation rewrite can blend it to λ.
    d = Symbolics.derivative(homotopy(a^2, b), a)
    @test has_homotopy(d)
    # numeric endpoint: compiled evaluation folds homotopy(1,0) via the runtime
    # fallback, so the derivative at a=3 must equal actual's derivative 6.
    # (`substitute` does not re-fold a node whose arguments are already constant,
    # so the genuine runtime path — codegen + numeric fallback — is exercised.)
    fd = Symbolics.build_function(d, a, b; expression = Val(false))
    @test fd(3.0, 7.0) ≈ 6.0
end

@testset "registered derivative w.r.t. the simplified argument (∂₂)" begin
    @variables a b λ
    # d/db homotopy(a, b^2) = homotopy(0,1) * 2b — differentiation path through
    # the SECOND argument exercises the ∂₂ = homotopy(0,1) registration.
    d2 = Symbolics.derivative(homotopy(a, b^2), b)
    @test has_homotopy(d2)
    # runtime: actual = a does not depend on b, so the compiled derivative is 0.
    fd2 = Symbolics.build_function(d2, a, b; expression = Val(false))
    @test fd2(3.0, 2.0) == 0.0
    # blend side: lowering the factor homotopy(0,1) → (1-λ)*1 + λ*0 = 1-λ yields
    # the true derivative of the blended expression through the simplified
    # branch: d/db [(1-λ)b^2 + λa] = (1-λ)*2b = (1-0.3)*2*2 = 2.8.
    blended = _rewrite_with_lambda(Symbolics.unwrap(d2), Symbolics.unwrap(λ))
    fblend = Symbolics.build_function(blended, λ, b; expression = Val(false))
    @test fblend(0.3, 2.0) ≈ 2.8
end

@testset "symbolic jacobian of a homotopy system builds and matches actual" begin
    @variables x(t)
    @named sys = System([D(x) ~ homotopy(-x^2, -x)], t)
    csys = complete(sys)
    J = ModelingToolkitBase.calculate_jacobian(csys)
    # at runtime homotopy ≡ actual: J = d(-x^2)/dx = -2x; check numerically at x=2
    fJ = Symbolics.build_function(only(J), unknowns(csys)[1]; expression = Val(false))
    @test fJ(2.0) ≈ -4.0
end

@testset "lower_homotopy: one shared λ across equations and observed" begin
    @variables x(t) y(t) z(t)
    @named sys = System(
        [
            D(x) ~ homotopy(-x, -x + 1),
            0 ~ homotopy(atan(z - 3), z),
        ], t;
        observed = [y ~ homotopy(x^2, x)]
    )
    # un-`complete`d systems are rejected (the pass reconstructs toplevel fields)
    @test_throws ArgumentError lower_homotopy(sys)
    csys = complete(sys)
    shadow, λ = lower_homotopy(csys)
    λn = Symbolics.wrap(λ)
    # no homotopy nodes remain anywhere
    @test !has_any_homotopy(shadow)
    # single λ identity: substituting THE returned λ resolves every blend, at
    # BOTH endpoints, in EVERY equation and in observed — a per-equation λ
    # would leave at least one blend unresolved here.
    eq1 = equations(shadow)[1]
    @test isequal(Symbolics.substitute(eq1.rhs, Dict(λn => 1.0)), -x)
    @test isequal(Symbolics.substitute(eq1.rhs, Dict(λn => 0.0)), -x + 1)
    eq2 = equations(shadow)[2]
    @test isequal(Symbolics.substitute(eq2.rhs, Dict(λn => 1.0)), atan(z - 3))
    @test isequal(Symbolics.substitute(eq2.rhs, Dict(λn => 0.0)), z)
    obs = only(observed(shadow))
    @test isequal(Symbolics.substitute(obs.rhs, Dict(λn => 1.0)), x^2)
    @test isequal(Symbolics.substitute(obs.rhs, Dict(λn => 0.0)), x)
    # λ is NOT a parameter of the shadow
    @test !any(p -> isequal(Symbolics.unwrap(p), λ), Symbolics.unwrap.(parameters(shadow)))
end

@testset "fixed sentinel λ does not capture a near-name user symbol" begin
    @variables x(t)
    @parameters __homotopy_λ   # close to, but NOT, the `__homotopy_λₘₜₖ` sentinel
    @named sys = System([D(x) ~ homotopy(-x * __homotopy_λ, -x + 1)], t)
    csys = complete(sys)
    shadow, λ = lower_homotopy(csys)
    λn = Symbolics.wrap(λ)
    userλ = Symbolics.unwrap(__homotopy_λ)
    # the sentinel λ is a distinct symbol (Sym equality is name-keyed)
    @test !isequal(λ, userλ)
    @test isequal(λ, ModelingToolkitBase.HOMOTOPY_LAMBDA)
    # the user's parameter survives lowering; the sentinel λ is not a parameter
    @test any(
        p -> isequal(Symbolics.unwrap(p), userλ),
        Symbolics.unwrap.(parameters(shadow))
    )
    @test !any(p -> isequal(Symbolics.unwrap(p), λ), Symbolics.unwrap.(parameters(shadow)))
    rhs = only(equations(shadow)).rhs
    # substituting the USER symbol does NOT resolve the blend: the result still
    # depends on the sentinel λ (its endpoints differ)
    after_user = Symbolics.substitute(rhs, Dict(__homotopy_λ => 1.0))
    @test !isequal(
        Symbolics.substitute(after_user, Dict(λn => 1.0)),
        Symbolics.substitute(after_user, Dict(λn => 0.0))
    )
    # only the sentinel λ resolves it: λ=1 → actual, λ=0 → simplified
    @test isequal(Symbolics.substitute(rhs, Dict(λn => 1.0)), -x * __homotopy_λ)
    @test isequal(Symbolics.substitute(rhs, Dict(λn => 0.0)), -x + 1)
end

@testset "guard: a user symbol colliding with the sentinel name errors" begin
    # A user symbol spelled EXACTLY like the reserved sentinel would otherwise be
    # silently rewritten into a state/parameter access, dropping the λ-dependence
    # of the residual. The loud guard turns that into an error.
    @variables x(t)
    @parameters __homotopy_λₘₜₖ
    @named sys = System([D(x) ~ homotopy(-x * __homotopy_λₘₜₖ, -x + 1)], t)
    csys = complete(sys)
    @test_throws ArgumentError lower_homotopy(csys)
end

@testset "broadcast: homotopy.(actual, simplified) lowers elementwise, one shared λ" begin
    @variables x(t)[1:2]
    actual = [atan(x[1] - 3), atan(x[2] - 5)]
    simplified = [x[1], x[2]]
    eqs = collect(0 .~ homotopy.(actual, simplified))
    @named sys = System(eqs, t)
    csys = complete(sys)
    @test has_any_homotopy(csys)        # one homotopy node survives per element
    shadow, λ = lower_homotopy(csys)
    λn = Symbolics.wrap(λ)
    @test !has_any_homotopy(shadow)     # every element lowered
    # the SAME λ resolves both elementwise blends, at both endpoints
    for (i, eq) in enumerate(equations(shadow))
        @test isequal(Symbolics.substitute(eq.rhs, Dict(λn => 1.0)), actual[i])
        @test isequal(Symbolics.substitute(eq.rhs, Dict(λn => 0.0)), simplified[i])
    end
end

@testset "nested homotopy lowers through lower_homotopy with one λ" begin
    @variables w(t)
    @named sys = System([0 ~ homotopy(homotopy(w^3, w^2), w)], t)
    csys = complete(sys)
    shadow, λ = lower_homotopy(csys)
    λn = Symbolics.wrap(λ)
    @test !has_any_homotopy(shadow)
    rhs = only(equations(shadow)).rhs
    # the SAME returned λ resolves both nesting levels
    @test isequal(Symbolics.substitute(rhs, Dict(λn => 1.0)), w^3)
    @test isequal(Symbolics.substitute(rhs, Dict(λn => 0.0)), w)
    @test Symbolics.value(Symbolics.substitute(rhs, Dict(λn => 1.0, w => 2.0))) ≈ 8.0
    @test Symbolics.value(Symbolics.substitute(rhs, Dict(λn => 0.0, w => 2.0))) ≈ 2.0
end

@testset "nested homotopy rewrites recursively with the same λ" begin
    # worker-level isolation test for `_rewrite_with_lambda`; the system-level
    # nesting path through the public API is the testset above.
    @variables w λtest
    shadowex = ModelingToolkitBase._rewrite_with_lambda(
        Symbolics.unwrap(homotopy(homotopy(w^3, w^2), w)), Symbolics.unwrap(λtest)
    )
    @test Symbolics.value(
        Symbolics.substitute(
            Symbolics.wrap(shadowex), Dict(λtest => 1.0, w => 2.0)
        )
    ) ≈ 8.0
    @test Symbolics.value(
        Symbolics.substitute(
            Symbolics.wrap(shadowex), Dict(λtest => 0.0, w => 2.0)
        )
    ) ≈ 2.0
end

@testset "generate_homotopy_residual guards: raw nodes and time-dependence error" begin
    # Both misuse modes were probed to fail SILENTLY with wrong numerics rather than
    # erroring; the guards turn them loud. These lock that contract.
    #
    # Guard 1: a system still carrying raw homotopy(...) nodes (i.e. NOT run through
    # `lower_homotopy`) would codegen a λ-independent residual via the numeric
    # fallback — caught before any code is generated.
    @variables x
    @named raw = System([0 ~ homotopy(x^2 - 2, x)])
    craw = complete(raw)
    @test has_any_homotopy(craw)
    @test_throws ArgumentError ModelingToolkitBase.generate_homotopy_residual(
        craw, ModelingToolkitBase.HOMOTOPY_LAMBDA
    )
    # Guard 2: a correctly-lowered but time-dependent shadow — codegen would emit
    # `f(u, p, t)` and silently misuse the λ slot as `t`.
    @variables z(t)
    @named tdep = System([D(z) ~ homotopy(-z, -z + 1)], t)
    shadow, λ = lower_homotopy(complete(tdep))
    @test !has_any_homotopy(shadow)          # lowered: the raw-node guard would pass
    @test_throws ArgumentError ModelingToolkitBase.generate_homotopy_residual(shadow, λ)
end
