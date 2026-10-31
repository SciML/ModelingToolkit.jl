using ModelingToolkit, SciMLBase, SymbolicIndexingInterface, LinearAlgebra, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

const FWW = SciMLBase.FunctionWrappersWrappers.FunctionWrappersWrapper
const RGF = ModelingToolkit.RuntimeGeneratedFunctions.RuntimeGeneratedFunction

# Whether an object of type `T` is reachable from `x`. Modules and types are not entered.
function reaches(::Type{T}, x, seen = IdDict{Any, Nothing}()) where {T}
    x isa T && return true
    (x isa Module || x isa Type || x isa Symbol || x isa String || isbits(x)) && return false
    haskey(seen, x) && return false
    seen[x] = nothing
    if x isa AbstractArray
        return any(i -> isassigned(x, i) && reaches(T, x[i], seen), eachindex(x))
    elseif x isa AbstractDict
        return any(((k, v),) -> reaches(T, k, seen) || reaches(T, v, seen), x)
    end
    return any(i -> isdefined(x, i) && reaches(T, getfield(x, i), seen), 1:fieldcount(typeof(x)))
end

# Whether the optimized code of `f` for `argtypes` refers to a function whose name contains one of `names`.
function refers_to(f, argtypes, names)
    code = only(code_typed(f, argtypes; optimize = true)).first.code
    named(n) = any(s -> occursin(s, string(n)), names)
    function refers(x)
        x isa GlobalRef && return named(x.name)
        x isa Core.CodeInstance && return refers(x.def)
        x isa Core.MethodInstance && return x.def isa Method && named(x.def.name)
        x isa Function && return named(nameof(x))
        x isa Expr && return any(refers, x.args)
        return false
    end
    return any(refers, code)
end

# the prepared Jacobians and the observed function, which a call evaluates
prepared(lf) = (lf.h, lf.uf_jac, lf.h_jac, lf.pf_jac, lf.hp_jac)
# the Jacobians of a call, as the bits of their entries
jacobians(j) = map(m -> map(x -> reinterpret(UInt64, x), m), (j.f_x, j.f_z, j.g_x, j.g_z, j.f_u, j.g_u, j.h_x, j.h_z, j.h_u))

function plant(; name)
    @variables x(t) = 0.0 u(t) = 0.0 [input = true] y(t) [output = true]
    @parameters k = 1.0
    return System([D(x) ~ -k * x + u, y ~ x], t; name)
end
function pi_controller(; name)
    @variables I(t) = 0.0 y(t) = 0.0 [input = true] u(t) [output = true]
    @parameters kp = 1.0 ki = 0.5
    return System([D(I) ~ -y, u ~ -kp * y + ki * I], t; name)
end

@variables x(t) = 0.0 u(t) = 0.0
@parameters k = 1.0
@named scalar_plant = System([D(x) ~ -k * x + u], t)

@testset "linearization_function with specialize = FullSpecialize" begin
    kwargs = (; op = Dict(x => 0.0, u => 0.0), warn_empty_op = false)
    for initialize in (false, true)
        lf, sys = linearization_function(scalar_plant, [u], [x]; initialize, kwargs...)
        lfs, syss = linearization_function(
            scalar_plant, [u], [x]; initialize, specialize = SciMLBase.FullSpecialize, kwargs...
        )
        @test SciMLBase.specialization(lfs.prob.f) === SciMLBase.FullSpecialize
        @test any(f -> reaches(FWW, f), prepared(lf))
        @test !any(f -> reaches(FWW, f), prepared(lfs))
        set, sets = setp(sys, k), setp(syss, k)
        for kv in (0.5, 1.2, 0.2)
            set(lf.prob.p, kv)
            sets(lfs.prob.p, kv)
            j = lf(lf.prob.u0, lf.prob.p, 0.0)
            js = lfs(lfs.prob.u0, lfs.prob.p, 0.0)
            @test jacobians(js) == jacobians(j)
            @test js.f_x == fill(-kv, 1, 1)
        end
        T = Base.return_types(lfs, Tuple{typeof(lfs.prob.u0), typeof(lfs.prob.p), Float64})
        @test isconcretetype(only(T))
    end
end

@testset "eval_expression reaches the problem's functions" begin
    kwargs = (; op = Dict(x => 0.0, u => 0.0), initialize = false, warn_empty_op = false)
    lf, _ = linearization_function(scalar_plant, [u], [x]; specialize = SciMLBase.FullSpecialize, kwargs...)
    lfe, _ = linearization_function(
        scalar_plant, [u], [x]; specialize = SciMLBase.FullSpecialize, eval_expression = true,
        eval_module = @__MODULE__, kwargs...
    )
    @test reaches(RGF, lf.prob.f)
    @test !reaches(RGF, lfe.prob.f)
    @test !any(f -> reaches(RGF, f), prepared(lfe))
    @test jacobians(lfe(lfe.prob.u0, lfe.prob.p, 0.0)) == jacobians(lf(lf.prob.u0, lf.prob.p, 0.0))
    # Without `FullSpecialize`, `eval_expression` applies to the observed function only, so the
    # system function the problem wraps is still built from `RuntimeGeneratedFunction`s.
    unwrapped(f) = f isa FWW ? unwrapped(first(f.fw).obj[]) : f isa SciMLBase.Void ? unwrapped(f.f) : f
    has_rgf_type(f) = occursin("RuntimeGeneratedFunction", string(typeof(unwrapped(f))))
    @test !has_rgf_type(lfe.prob.f.f)
    for option in ((;), (; specialize = nothing))
        lfd, _ = linearization_function(
            scalar_plant, [u], [x]; eval_expression = true, eval_module = @__MODULE__,
            kwargs..., option...
        )
        @test has_rgf_type(lfd.prob.f.f)
    end
end

@testset "get_sensitivity_function with specialize = FullSpecialize" begin
    @named P = plant()
    @named C = pi_controller()
    model = System(
        [connect(P.y, C.y), connect(C.u, :plant_input, P.u)], t; systems = [P, C], name = :loop
    )
    kwargs = (; op = Dict(P.x => 0.0, C.I => 0.0), initialize = false, warn_empty_op = false)
    lf, sys = get_sensitivity_function(model, model.plant_input; kwargs...)
    lfs, syss = get_sensitivity_function(
        model, model.plant_input; specialize = SciMLBase.FullSpecialize, kwargs...
    )
    @test !any(f -> reaches(FWW, f), prepared(lfs))
    for (kp, ki) in ((1.0, 0.5), (2.5, 0.1), (0.3, 3.0))
        setp(sys, sys.C.kp)(lf.prob.p, kp)
        setp(sys, sys.C.ki)(lf.prob.p, ki)
        setp(syss, syss.C.kp)(lfs.prob.p, kp)
        setp(syss, syss.C.ki)(lfs.prob.p, ki)
        js = lfs(lfs.prob.u0, lfs.prob.p, 0.0)
        @test jacobians(js) == jacobians(lf(lf.prob.u0, lf.prob.p, 0.0))
        # Closed loop `ẋ = -x + u`, `İ = -x`, `u = -kp * x + ki * I`: `A` has trace `-(1 + kp)` and
        # determinant `ki` in any state order, and the sensitivity's feedthrough is 1.
        @test isempty(js.g_z)
        @test tr(js.f_x) ≈ -(1 + kp)
        @test det(js.f_x) ≈ ki
        @test js.h_u == ones(1, 1)
    end
end

# With `FullSpecialize` and `NoInit`, a call whose `u` and `p` have the problem's types does not
# remake the problem. With `OverrideInit` (`initialize = true`), or without `FullSpecialize`, the
# same call still does. All give the same Jacobians.
@testset "NoInit call without remake" begin
    kwargs = (; op = Dict(x => 0.5, u => 0.0), warn_empty_op = false)
    full = (; specialize = SciMLBase.FullSpecialize)
    lf, _ = linearization_function(scalar_plant, [u], [x]; initialize = false, kwargs..., full...)
    lfi, _ = linearization_function(scalar_plant, [u], [x]; initialize = true, kwargs..., full...)
    lfd, _ = linearization_function(scalar_plant, [u], [x]; initialize = false, kwargs...)
    sig(l) = Tuple{typeof(l.prob.u0), typeof(l.prob.p), Float64}
    names = ("remake", "get_initial_values")
    @test !refers_to(lf, sig(lf), names)
    @test refers_to(lfi, sig(lfi), names)
    @test refers_to(lfd, sig(lfd), names)
    j = lf(lf.prob.u0, lf.prob.p, 0.0)
    ji = lfi(lfi.prob.u0, lfi.prob.p, 0.0)
    jd = lfd(lfd.prob.u0, lfd.prob.p, 0.0)
    @test jacobians(j) == jacobians(ji)
    @test jacobians(j) == jacobians(jd)
    @test j.x == ji.x == jd.x
end

# The Jacobians of `D(x) ~ -k * x + u` with output `x` are `f_x = -k`, `f_u = 1`, `h_x = 1`, `h_u = 0`.
@testset "Analytic Jacobians for each specialization, initialization and evaluation setting" begin
    @variables x(t) = 0.0 u(t) = 0.0
    @parameters k = 1.0
    @named reference_plant = System([D(x) ~ -k * x + u], t)
    for initialize in (false, true), eval_expression in (false, true)
        eval_kwargs = eval_expression ? (; eval_expression, eval_module = @__MODULE__) : (;)
        common = (; op = Dict(x => 0.0, u => 0.0), initialize, warn_empty_op = false, eval_kwargs...)
        for option in ((;), (; specialize = nothing), (; specialize = SciMLBase.FullSpecialize))
            lf, sys = linearization_function(reference_plant, [u], [x]; common..., option...)
            setter = setp(sys, k)
            for kv in (0.5, 1.2, 0.2)
                setter(lf.prob.p, kv)
                # functions evaluated with `eval_expression` exist only in later worlds
                j = Base.invokelatest(lf, lf.prob.u0, lf.prob.p, 0.0)
                @test j.f_x == fill(-kv, 1, 1)
                @test j.f_u == ones(1, 1)
                @test j.h_x == ones(1, 1)
                @test j.h_u == zeros(1, 1)
            end
        end
    end
end
