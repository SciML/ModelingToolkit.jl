using ModelingToolkitBase
using ModelingToolkitBase: limited, limitnew, limitold, has_limited, has_any_limited,
    lower_limited, strip_limited_system, LimitedCtx, LIMIT_NEW, LIMIT_OLD,
    t_nounits as t, D_nounits as D
using NonlinearSolve
import NonlinearSolveBase
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

# NonlinearSolve applies `NonlinearFunction.postcondition` only from versions carrying
# the iterate-commit hooks; older releases silently ignore the field, so the solver
# behavior assertions are gated (functional assertions still run everywhere).
const NLS_APPLIES_POSTCONDITION = isdefined(NonlinearSolveBase, :apply_postcondition!!)

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
    @test prob.f.postcondition !== nothing
    sol = solve(prob, NewtonRaphson(); maxiters = 1000)
    @test SciMLBase.successful_retcode(sol)
    @test abs(sol[v] - 0.6698509496766559) < 1.0e-6

    if NLS_APPLIES_POSTCONDITION
        eqs_plain = [0 ~ (v - Vs) / R + Is * (exp(v / Vt) - 1)]
        @mtkcompile sys_plain = System(eqs_plain)
        prob_plain = NonlinearProblem(
            sys_plain,
            [v => 0.0, Vs => 5.0, R => 1.0e3, Is => 1.0e-14, Vt => 0.025]
        )
        sol_plain = solve(prob_plain, NewtonRaphson(); maxiters = 1000)
        @test sol.stats.nsteps < sol_plain.stats.nsteps ÷ 4
    end
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
    @test prob.f.postcondition !== nothing
    sol = solve(prob, NewtonRaphson(); maxiters = 1000)
    @test SciMLBase.successful_retcode(sol)
    @test abs(sol[diode.v] - 0.6698509496766559) < 1.0e-5
    @test abs(sol[diode.i] - 0.004330149050323345) < 1.0e-8
    if NLS_APPLIES_POSTCONDITION
        @test sol.stats.nsteps < 40
    end
end

@testset "time-dependent systems strip the operator" begin
    @variables v(t) = 0.0
    @parameters Is = 1.0e-14 Vt = 0.025 C = 1.0e-6 R = 1.0e3 Vsrc = 5.0
    pvcrit = 0.71
    eqs = [
        D(v) ~ (
            (Vsrc - v) / R -
                Is * (exp(limited(v, pnjlim(limitnew, limitold, Vt, pvcrit)) / Vt) - 1)
        ) / C,
    ]
    @named tsys = System(eqs, t)
    ctsys = mtkcompile(tsys)
    @test !has_any_limited(ctsys)
    @test getmetadata(ctsys, LimitedCtx, nothing) === nothing
    oprob = ODEProblem(ctsys, [], (0.0, 1.0e-5))
    du = [0.0]
    oprob.f(du, oprob.u0, oprob.p, 0.0)
    @test du[1] ≈ (5.0 / 1.0e3) / 1.0e-6 rtol = 1.0e-6
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
