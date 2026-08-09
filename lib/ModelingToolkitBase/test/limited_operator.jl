using ModelingToolkitBase
using ModelingToolkitBase: limited, limitnew, limitold, has_limited, has_any_limited,
    lower_limited, strip_limited_system, LimitedCtx, StageLimitedCtx,
    attach_stage_limiters, bounds_limiter_specs, LIMIT_NEW, LIMIT_OLD,
    t_nounits as t, D_nounits as D
using NonlinearSolve
import NonlinearSolveBase
using OrdinaryDiffEqSDIRK
# `NonlinearSolveAlg` — the only stepper option that consumes `ODEFunction.nlstep_data` —
# is not exported by any package that declares it.
import OrdinaryDiffEqNonlinearSolve
using SciMLBase
using Symbolics
using SymbolicUtils: getmetadata
using Test

# The classic SPICE3 junction-voltage limiting rule, registered as an opaque function
# so its branches keep Julia short-circuit semantics.
function pnjlim(vnew, vold, vt, vcrit)
    if vnew > vcrit && abs(vnew - vold) > 2vt
        if vold > 0
            arg = 1 + (vnew - vold) / vt
            vnew = arg > 0 ? vold + vt * log(arg) : vcrit
        else
            vnew = vt * log(vnew / vt)
        end
    end
    return vnew
end
@register_symbolic pnjlim(vnew, vold, vt, vcrit)

@testset "runtime numeric fallback is the actual branch" begin
    @test ModelingToolkitBase.limited(2.0, 5.0) == 2.0
end

@testset "derivative treats the operator as actual" begin
    @variables a b
    d = Symbolics.derivative(limited(a^2, b), a)
    @test !has_limited(d)
    @test isequal(Symbolics.simplify(d), Symbolics.simplify(2a))
    @test isequal(Symbolics.derivative(limited(a^2, b), b), 0)
end

@testset "lowering: augmentation, irreducibility, guesses, registry" begin
    @variables v
    @parameters Vs R Is Vt vcrit
    eqs = [
        0 ~ (v - Vs) / R +
            Is * (exp(limited(v, pnjlim(limitnew, limitold, Vt, vcrit)) / Vt) - 1),
    ]
    @named rawsys = System(eqs)
    @test has_any_limited(rawsys)
    low = lower_limited(rawsys)
    @test !has_any_limited(low)
    dvs = unknowns(low)
    @test length(equations(low)) == 2
    lv = only(filter(u -> occursin("limited_1", string(u)), dvs))
    @test ModelingToolkitBase.isirreducible(lv)
    # sentinels are removed from the discovered variables
    @test !any(
        u -> occursin("limitnew", string(u)) || occursin("limitold", string(u)),
        string.(dvs)
    )
    @test !any(
        p -> occursin("limitnew", string(p)) || occursin("limitold", string(p)),
        string.(ModelingToolkitBase.get_ps(low))
    )
    specs = getmetadata(low, LimitedCtx, nothing)
    @test specs !== nothing && length(specs) == 1
end

@testset "nested limited operators are rejected" begin
    @variables x
    @named sys = System([0 ~ limited(limited(x, limitnew), limitnew) - 1])
    @test_throws ArgumentError lower_limited(sys)
end

@testset "compiled system solves with automatic PCNR limiting" begin
    @variables v
    @parameters Vs R Is Vt vcrit
    eqs = [
        0 ~ (v - Vs) / R +
            Is * (exp(limited(v, pnjlim(limitnew, limitold, Vt, vcrit)) / Vt) - 1),
    ]
    @mtkcompile sys = System(eqs)
    # the irreducible limited quantity survives structural simplification as the
    # representative of its alias class
    @test any(u -> occursin("limited_1", string(u)), string.(unknowns(sys)))
    pvals = [
        Vs => 5.0, R => 1.0e3, Is => 1.0e-14, Vt => 0.025,
        vcrit => 0.025 * log(0.025 / (sqrt(2) * 1.0e-14)),
    ]
    prob = NonlinearProblem(sys, [v => 0.0; pvals])
    @test haskey(prob.kwargs, :postcondition)
    sol = solve(prob, NewtonRaphson(); maxiters = 1000)
    @test SciMLBase.successful_retcode(sol)
    @test abs(sol[v] - 0.6698509496766559) < 1.0e-6

    # a corrector passed at solve time replaces the model-declared one
    calls = Ref(0)
    solve(
        prob, NewtonRaphson();
        postcondition = (up, uprev, p, cache) -> (calls[] += 1; nothing), maxiters = 1000
    )
    @test calls[] > 0

    eqs_plain = [0 ~ (v - Vs) / R + Is * (exp(v / Vt) - 1)]
    @mtkcompile sys_plain = System(eqs_plain)
    prob_plain = NonlinearProblem(
        sys_plain,
        [v => 0.0, Vs => 5.0, R => 1.0e3, Is => 1.0e-14, Vt => 0.025]
    )
    sol_plain = solve(prob_plain, NewtonRaphson(); maxiters = 1000)
    @test sol.stats.nsteps < sol_plain.stats.nsteps ÷ 4
end

@testset "hierarchical components: limiting composes through namespacing" begin
    function DCDiode(; name, Is = 1.0e-14, Vt = 0.025)
        @variables v i
        ps = @parameters begin
            (Is::Float64 = Is)
            (Vt::Float64 = Vt)
            (vcrit::Float64 = Vt * log(Vt / (sqrt(2) * Is)))
        end
        eqs = [i ~ Is * (exp(limited(v, pnjlim(limitnew, limitold, Vt, vcrit)) / Vt) - 1)]
        return System(eqs, [v, i], ps; name)
    end
    function DCResistor(; name, R = 1.0e3)
        @variables v i
        @parameters R = R
        return System([v ~ i * R], [v, i], [R]; name)
    end
    @named diode = DCDiode()
    @named res = DCResistor()
    @parameters Vs = 5.0
    conn = [res.i ~ diode.i, Vs ~ res.v + diode.v]
    @named circuit = System(conn, [], [Vs]; systems = [diode, res])
    csys = mtkcompile(circuit)
    prob = NonlinearProblem(
        csys, [diode.v => 0.0, diode.i => 0.0, res.v => 0.0, res.i => 0.0]
    )
    @test haskey(prob.kwargs, :postcondition)
    sol = solve(prob, NewtonRaphson(); maxiters = 1000)
    @test SciMLBase.successful_retcode(sol)
    @test abs(sol[diode.v] - 0.6698509496766559) < 1.0e-5
    @test abs(sol[diode.i] - 0.004330149050323345) < 1.0e-8
    @test sol.stats.nsteps < 40
end

# The transient RC-diode circuit used by the stage-limiting testsets below: a capacitor
# charged through a resistor and clamped by a diode, so the state is the junction voltage.
function rc_diode_system(; vcrit = 0.71, limiting = true)
    @variables v(t) = 0.0
    @parameters Is = 1.0e-14 Vt = 0.025 C = 1.0e-6 R = 1.0e3 Vsrc = 5.0
    vd = limiting ? limited(v, pnjlim(limitnew, limitold, Vt, vcrit)) : v
    eqs = [D(v) ~ ((Vsrc - v) / R - Is * (exp(vd / Vt) - 1)) / C]
    @named tsys = System(eqs, t)
    return v, tsys
end

@testset "time-dependent systems strip the operator but record it" begin
    v, tsys = rc_diode_system()
    ctsys = mtkcompile(tsys)
    @test !has_any_limited(ctsys)
    @test getmetadata(ctsys, LimitedCtx, nothing) === nothing
    specs = getmetadata(ctsys, StageLimitedCtx, nothing)
    @test specs !== nothing && length(specs) == 1
    @test isequal(only(specs)[1], ModelingToolkitBase.unwrap(v))
    # the limiter is canonicalized to the toplevel sentinels
    @test occursin("__limitnew_", string(only(specs)[2]))
    oprob = ODEProblem(ctsys, [], (0.0, 1.0e-5))
    du = [0.0]
    oprob.f(du, oprob.u0, oprob.p, 0.0)
    @test du[1] ≈ (5.0 / 1.0e3) / 1.0e-6 rtol = 1.0e-6
end

@testset "stage limiters: the limiter is conjugated onto the stage unknown" begin
    v, tsys = rc_diode_system()
    ctsys = mtkcompile(tsys)
    uv = ModelingToolkitBase.unwrap(v)

    # Mimic `inner_nlsystem`: the stage unknown `z` enters the right-hand side as the
    # physical value `γ₂ z + inner_tmp`, and the stage residual is that rhs minus `z`.
    @parameters g2 tmp
    subrules = Dict(uv => ModelingToolkitBase.unwrap(g2 * v + tmp))
    rhs = substitute(equations(ctsys)[1].rhs, subrules)
    stage = System(
        [0 ~ rhs - v], [v], [parameters(ctsys); g2; tmp]; name = :stage
    )
    stage = attach_stage_limiters(ctsys, complete(stage), subrules)

    specs = getmetadata(stage, LimitedCtx, nothing)
    @test specs !== nothing && length(specs) == 1
    @test isequal(only(specs)[1], uv)

    prob = NonlinearProblem(stage, [v => 0.0, g2 => 0.5, tmp => 0.1])
    @test haskey(prob.kwargs, :postcondition)
    post = prob.kwargs[:postcondition]

    # the corrector limits the physical quantity `0.5 z + 0.1` and writes back the `z`
    # that produces the limited value
    up = [3.0]
    post(up, [0.0], prob.p, nothing)
    @test up[1] ≈ (pnjlim(0.5 * 3.0 + 0.1, 0.5 * 0.0 + 0.1, 0.025, 0.71) - 0.1) / 0.5
    @test 0.5 * up[1] + 0.1 < 0.5 * 3.0 + 0.1

    # limiters are the identity at a fixed point, so converged iterates are untouched
    fixed = [0.9]
    post(fixed, [0.9], prob.p, nothing)
    @test fixed[1] ≈ 0.9
end

# The single-state stage system for a bounded transient model, built the way
# `inner_nlsystem` would: the physical quantity is `g2 * z + tmp`, the stage unknown is `z`.
function bounded_stage(; xbounds = (0.0, 1.0), limiting = false)
    @variables x(t) = 0.5 [bounds = xbounds]
    @parameters k = 1.0 g2 tmp
    xd = limiting ? limited(x, limitold + (limitnew - limitold) / 2) : x
    @named bsys = System([D(x) ~ -k * xd^2 * (1 - x)], t)
    cbsys = mtkcompile(bsys)
    subrules = Dict(
        ModelingToolkitBase.unwrap(x) => ModelingToolkitBase.unwrap(g2 * x + tmp)
    )
    rhs = substitute(equations(cbsys)[1].rhs, subrules)
    stage = System([0 ~ rhs - x], [x], [parameters(cbsys); g2; tmp]; name = :stage)
    return x, g2, tmp, cbsys, complete(stage), subrules
end

# `generate_ODENLStepData` builds the stage problem with the box explicitly disabled, since
# the stage unknowns are increments rather than states; mirror that here.
function stage_problem(stage, op)
    return NonlinearProblem(stage, op; lb = nothing, ub = nothing)
end

limiter_value(expr, new) = Symbolics.build_function(
    expr, [LIMIT_NEW]; expression = Val{false}
)([new])

@testset "bounds become a clamping limiter of the physical quantity" begin
    x, g2, tmp, cbsys, _, _ = bounded_stage()
    specs = bounds_limiter_specs(cbsys)
    @test length(specs) == 1
    @test isequal(only(specs)[1], ModelingToolkitBase.unwrap(x))
    clamp2 = only(specs)[2]
    @test limiter_value(clamp2, 3.0) == 1.0
    @test limiter_value(clamp2, -3.0) == 0.0
    @test limiter_value(clamp2, 0.5) == 0.5
    # an infinite side contributes no comparison at all
    _, _, _, half, _, _ = bounded_stage(; xbounds = (0.0, Inf))
    clamp1 = only(bounds_limiter_specs(half))[2]
    @test limiter_value(clamp1, 3.0) == 3.0
    @test limiter_value(clamp1, -3.0) == 0.0
    # and a wholly unbounded state contributes no corrector
    _, _, _, none, _, _ = bounded_stage(; xbounds = (-Inf, Inf))
    @test isempty(bounds_limiter_specs(none))
end

@testset "stage limiters: bounds are conjugated into stage coordinates" begin
    x, g2, tmp, cbsys, stage, subrules = bounded_stage()
    # the escape hatch turns the clamp off, but it is not what a user has to reach for
    @test getmetadata(
        attach_stage_limiters(cbsys, stage, subrules; limit_bounds = false),
        LimitedCtx, nothing
    ) === nothing

    stage = attach_stage_limiters(cbsys, stage, subrules)
    specs = getmetadata(stage, LimitedCtx, nothing)
    @test specs !== nothing && length(specs) == 1
    @test isequal(only(specs)[1], ModelingToolkitBase.unwrap(x))

    prob = stage_problem(stage, [x => 0.5, g2 => 0.5, tmp => 0.1])
    post = prob.kwargs[:postcondition]

    # the clamp acts on the physical quantity `0.5 z + 0.1`: z = 6 proposes 3.1, which is
    # projected to the upper bound 1.0, i.e. to the increment (1.0 - 0.1) / 0.5
    up = [6.0]
    post(up, [0.0], prob.p, nothing)
    @test up[1] ≈ 1.8
    @test 0.5 * up[1] + 0.1 ≈ 1.0
    # below the lower bound, symmetrically
    up = [-4.0]
    post(up, [0.0], prob.p, nothing)
    @test 0.5 * up[1] + 0.1 ≈ 0.0
    # and an in-range iterate is left exactly alone
    up = [1.0]
    post(up, [0.0], prob.p, nothing)
    @test up[1] == 1.0
end

@testset "stage limiters: a user limiter and bounds compose as clamp ∘ limiter" begin
    x, g2, tmp, cbsys, stage, subrules = bounded_stage(; limiting = true)
    stage = attach_stage_limiters(cbsys, stage, subrules)
    specs = getmetadata(stage, LimitedCtx, nothing)
    # both correctors are registered against the same stage unknown, in that order
    @test length(specs) == 2
    @test all(s -> isequal(s[1], ModelingToolkitBase.unwrap(x)), specs)

    prob = stage_problem(stage, [x => 0.5, g2 => 0.5, tmp => 0.1])
    post = prob.kwargs[:postcondition]

    # z = 6 proposes the physical value 3.1 against the previous 0.1. The halving limiter
    # takes it to 1.6, the clamp then projects that onto the upper bound 1.0. Applying the
    # two the other way round would give 0.55, so this pins the order.
    up = [6.0]
    post(up, [0.0], prob.p, nothing)
    @test 0.5 * up[1] + 0.1 ≈ 1.0
    # z = 2 proposes 1.1 -> halved to 0.6, which is in range and so survives the clamp;
    # the clamp alone would have produced 1.0.
    up = [2.0]
    post(up, [0.0], prob.p, nothing)
    @test 0.5 * up[1] + 0.1 ≈ 0.6
    # a converged iterate is a fixed point of both
    up = [0.8]
    post(up, [0.8], prob.p, nothing)
    @test up[1] ≈ 0.8
end

@testset "stage limiters: a box that cannot be conjugated is dropped, not an error" begin
    @variables x(t) = 0.5 [bounds = (0.0, 1.0)] y(t) = 0.5
    @parameters g2 tmp1 tmp2
    @named osys = System([D(x) ~ -x^2 * y, D(y) ~ -y^2], t)
    cosys = mtkcompile(osys)
    subrules = Dict(
        ModelingToolkitBase.unwrap(x) => ModelingToolkitBase.unwrap(g2 * x + tmp1),
        ModelingToolkitBase.unwrap(y) => ModelingToolkitBase.unwrap(g2 * y + tmp2)
    )
    # a stage system that has torn `x` away as a nonlinear function of `y`, leaving the box
    # on `x` with no single affine stage unknown to act on
    stage = complete(
        System(
            [0 ~ y - g2 * y^2 - tmp2], [y], [g2, tmp1, tmp2];
            observed = [x ~ y^2], name = :stage
        )
    )
    @test getmetadata(attach_stage_limiters(cosys, stage, subrules), LimitedCtx, nothing) ===
        nothing
end

@testset "explicit `lb`/`ub` of `nothing` suppress the metadata-derived box" begin
    @variables y [bounds = (0.0, 1.0)]
    @parameters a = 2.0
    @mtkcompile nsys = System([0 ~ y^2 - a * y + 0.5])
    derived = NonlinearProblem(nsys, [y => 0.3])
    @test derived.lb == [0.0] && derived.ub == [1.0]
    plain = NonlinearProblem(nsys, [y => 0.3]; lb = nothing, ub = nothing)
    @test plain.lb === nothing && plain.ub === nothing
    @test NonlinearProblem(nsys, [y => 0.3]; lb = [-1.0]).lb == [-1.0]
end

# Build the two-state stage system for `eqs`, the way `inner_nlsystem` would.
function two_state_stage(eqs, x, y, g2, tmp1, tmp2)
    @named bad = System(eqs, t)
    cbad = mtkcompile(bad)
    ux, uy = ModelingToolkitBase.unwrap(x), ModelingToolkitBase.unwrap(y)
    subrules = Dict(
        ux => ModelingToolkitBase.unwrap(g2 * x + tmp1),
        uy => ModelingToolkitBase.unwrap(g2 * y + tmp2)
    )
    rhss = [substitute(eq.rhs, subrules) for eq in equations(cbad)]
    stage = System(
        [0 ~ rhss[1] - x, 0 ~ rhss[2] - y], [x, y],
        [parameters(cbad); g2; tmp1; tmp2]; name = :stage
    )
    return cbad, complete(stage), subrules
end

@testset "stage limiters: quantities that cannot be conjugated are rejected" begin
    @parameters g2 tmp1 tmp2
    @variables x(t) = 1.0 y(t) = 1.0

    # spans two stage unknowns: there is no single iterate entry to correct
    cbad, stage, subrules = two_state_stage(
        [D(x) ~ -limited(x - y, limitnew / 2) - y, D(y) ~ -y - x], x, y, g2, tmp1, tmp2
    )
    err = try
        attach_stage_limiters(cbad, stage, subrules)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("exactly one unknown", sprint(showerror, err))

    # nonlinear in its stage unknown: correcting the quantity is not correcting the unknown
    cbad, stage, subrules = two_state_stage(
        [D(x) ~ -limited(x^2, limitnew / 2) - y, D(y) ~ -y - x], x, y, g2, tmp1, tmp2
    )
    err = try
        attach_stage_limiters(cbad, stage, subrules)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("affine", sprint(showerror, err))
end

# `nlstep = true` needs `generate_ODENLStepData`, which lives in ModelingToolkit.jl; this
# file is also run from the root package's suite, where it is loaded.
if @isdefined(ModelingToolkit)
    @testset "nlstep transient solve limits its stage Newton iterates" begin
        v, tsys = rc_diode_system()
        ctsys = mtkcompile(tsys)
        tspan = (0.0, 1.0e-4)

        prob = ODEProblem(ctsys, [], tspan; nlstep = true)
        nlstep = prob.f.nlstep_data
        @test nlstep !== nothing
        @test haskey(nlstep.nlprob.kwargs, :postcondition)
        # the stage system is still the ODE unknowns: no auxiliary unknown was introduced
        @test length(nlstep.u0perm) == length(unknowns(ctsys))
        @test nlstep.u0perm == [1]
        # and the corrector survives into the stage solver's cache, which is where the
        # implicit stepper drives it from
        stage_cache = init(nlstep.nlprob, NewtonRaphson())
        @test NonlinearSolveBase.get_postcondition(stage_cache) isa
            ModelingToolkitBase.LimitedPostcondition

        alg = ImplicitEuler(
            nlsolve = OrdinaryDiffEqNonlinearSolve.NonlinearSolveAlg()
        )
        sol = solve(prob, alg; dt = 1.0e-6, adaptive = false)
        @test SciMLBase.successful_retcode(sol)

        # limiting damps the stage iteration; it must not move the trajectory, so the
        # unannotated model integrated the same way has to agree
        vp, psys = rc_diode_system(; limiting = false)
        cpsys = mtkcompile(psys)
        @test getmetadata(cpsys, StageLimitedCtx, nothing) === nothing
        prob_plain = ODEProblem(cpsys, [], tspan; nlstep = true)
        @test !haskey(prob_plain.f.nlstep_data.nlprob.kwargs, :postcondition)
        sol_plain = solve(prob_plain, alg; dt = 1.0e-6, adaptive = false)
        @test SciMLBase.successful_retcode(sol_plain)
        @test sol[v, end] ≈ sol_plain[vp, end] rtol = 1.0e-6
    end

    @testset "nlstep rejects limiters it cannot attach" begin
        v, tsys = rc_diode_system()
        ctsys = mtkcompile(tsys)
        # the SCC decomposition splits the problem the corrector is attached to
        @test_throws ArgumentError ODEProblem(
            ctsys, [], (0.0, 1.0e-4); nlstep = true, nlstep_scc = true
        )
        # a `bounds` box, on the other hand, is nobody's request: it must not turn a model
        # that decomposes perfectly well into a construction error
        @variables x(t) = 0.5 [bounds = (0.0, 1.0)]
        @parameters k = 1.0
        @mtkcompile bsys = System([D(x) ~ -k * x^2 * (1 - x)], t)
        sccprob = ODEProblem(
            bsys, [], (0.0, 1.0); nlstep = true, nlstep_scc = true
        )
        @test sccprob.f.nlstep_data !== nothing
    end

    # A nonlinear right-hand side is essential here: `mtkcompile` solves a linear stage
    # system exactly and tears it to zero unknowns, so a linear model has no stage unknown
    # for a box to be attached to in the first place.
    function bounded_ode(; xbounds = (0.0, 1.0))
        @variables x(t) = 0.5 [bounds = xbounds]
        @parameters k = 1.0
        @named bsys = System([D(x) ~ -k * x^2 * (1 - x)], t)
        return x, mtkcompile(bsys)
    end

    @testset "nlstep does not put the state's box on the stage increments" begin
        x, bsys = bounded_ode()
        prob = ODEProblem(bsys, [], (0.0, 1.0); nlstep = true)
        nlprob = prob.f.nlstep_data.nlprob
        # the stage system really does iterate on a variable carrying the metadata
        @test length(unknowns(nlprob.f.sys)) == 1
        @test ModelingToolkitBase.getbounds(only(unknowns(nlprob.f.sys))) == (0.0, 1.0)
        # ... but the box describes `γ₂ z + inner_tmp`, not `z`, so it must not be here
        @test nlprob.lb === nothing
        @test nlprob.ub === nothing
        # which is what keeps the stage problem in the raw coordinates the implicit
        # stepper writes into `nlprob.u0`
        @test !NonlinearSolveBase.needs_bounds_transform(nlprob, NewtonRaphson())
        # the escape hatch is the only way to end up with neither
        noclamp = ODEProblem(
            bsys, [], (0.0, 1.0); nlstep = true, nlstep_limit_bounds = false
        )
        @test !haskey(noclamp.f.nlstep_data.nlprob.kwargs, :postcondition)
    end

    @testset "nlstep delivers the box as a per-stage clamp, automatically" begin
        x, bsys = bounded_ode()
        prob = ODEProblem(bsys, [], (0.0, 1.0); nlstep = true)
        nd = prob.f.nlstep_data
        nlprob = nd.nlprob
        @test nlprob.lb === nothing && nlprob.ub === nothing
        @test !NonlinearSolveBase.needs_bounds_transform(nlprob, NewtonRaphson())
        @test haskey(nlprob.kwargs, :postcondition)
        post = nlprob.kwargs[:postcondition]

        # `a` and `b` of the conjugation are the live stage parameters, so driving them the
        # way the stepper does moves the box the corrector enforces on the increment
        nd.set_γ_c(nlprob, (1.0, 0.5, 1.0, 0.0))
        nd.set_inner_tmp(nlprob, [0.1])
        up = [6.0]
        post(up, [0.0], nlprob.p, nothing)
        @test 0.5 * up[1] + 0.1 ≈ 1.0
        @test up[1] ≈ 1.8

        nd.set_γ_c(nlprob, (1.0, 0.25, 1.0, 0.0))
        nd.set_inner_tmp(nlprob, [-0.5])
        up = [8.0]
        post(up, [0.0], nlprob.p, nothing)
        @test 0.25 * up[1] - 0.5 ≈ 1.0
        @test up[1] ≈ 6.0
        up = [-6.0]
        post(up, [0.0], nlprob.p, nothing)
        @test 0.25 * up[1] - 0.5 ≈ 0.0

        alg = ImplicitEuler(nlsolve = OrdinaryDiffEqNonlinearSolve.NonlinearSolveAlg())
        sol = solve(prob, alg; dt = 0.05, adaptive = false)
        @test SciMLBase.successful_retcode(sol)
        # clamping intermediate iterates of a solve whose stage solutions are inside the
        # box changes how the step is found, not what it is
        _, psys = bounded_ode(; xbounds = (-Inf, Inf))
        prob_plain = ODEProblem(psys, [], (0.0, 1.0); nlstep = true)
        @test !haskey(prob_plain.f.nlstep_data.nlprob.kwargs, :postcondition)
        sol_plain = solve(prob_plain, alg; dt = 0.05, adaptive = false)
        @test SciMLBase.successful_retcode(sol_plain)
        # the two stage iterations take different paths to the same root, so they agree to
        # the stage solver's tolerance rather than to machine precision
        @test sol[x, end] ≈ sol_plain[unknowns(psys)[1], end] rtol = 1.0e-6
    end
end

@testset "guards" begin
    @variables x y
    # un-lowered limited nodes are rejected at NonlinearFunction construction
    @named gsys = System([0 ~ limited(x, limitnew) - 1])
    csys = complete(gsys)
    @test_throws ArgumentError NonlinearFunction(csys)
    # limiter referencing other unknowns is rejected at problem construction
    @named bsys = System([0 ~ limited(x, limitnew + y) - 1, 0 ~ y^2 + x - 3])
    sys2 = mtkcompile(bsys)
    err = try
        NonlinearProblem(sys2, [x => 1.0, y => 1.0])
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("unknowns", sprint(showerror, err))
end
