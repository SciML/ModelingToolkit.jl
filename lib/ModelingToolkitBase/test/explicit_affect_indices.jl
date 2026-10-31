using ModelingToolkitBase, SciMLBase, OrdinaryDiffEqDefault, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using SymbolicUtils: BasicSymbolic
using SymbolicIndexingInterface: variable_index

# Whether `f`'s optimized code for `sig` makes a call that isn't resolved statically.
function has_dynamic_call(f, sig)
    code = only(code_typed(f, sig; optimize = true)).first.code
    late = (:invokelatest, :_call_latest, :_call_in_world)
    return any(code) do st
        st isa Expr || return false
        if st.head === :call
            g = st.args[1]
            g isa GlobalRef || return true
            g.name in late && return true
            v = getglobal(g.mod, g.name)
            return !(v isa Core.Builtin || v isa Core.IntrinsicFunction)
        elseif st.head === :invoke
            ci = st.args[1]
            mi = ci isa Core.CodeInstance ? ci.def : ci
            return mi.def isa Method && mi.def.name in late
        end
        return false
    end
end
# a dispatcher whose affect type is unknown to inference
dynamic_nth_affect(affects, k, integ) =
    ModelingToolkitBase.apply_affect_if_real(Base.inferencebarrier(affects[k]), integ)

# The problem's one continuous callback: alone, or in a `CallbackSet` with discrete ones.
function continuous_callback(prob)
    cb = prob.kwargs[:callback]
    return cb isa SciMLBase.CallbackSet ? only(cb.continuous_callbacks) : cb
end

# Every object reachable from `x` that is a symbolic term. Modules and types are not entered.
function symbolic_objects(x, seen = IdDict{Any, Nothing}(), found = Any[])
    (x isa Module || x isa Type || x isa Symbol || x isa String || isbits(x)) && return found
    haskey(seen, x) && return found
    seen[x] = nothing
    x isa BasicSymbolic && (push!(found, x); return found)
    if x isa AbstractArray
        for i in eachindex(x)
            isassigned(x, i) && symbolic_objects(x[i], seen, found)
        end
    elseif x isa AbstractDict
        for (k, v) in x
            symbolic_objects(k, seen, found)
            symbolic_objects(v, seen, found)
        end
    else
        for i in 1:fieldcount(typeof(x))
            isdefined(x, i) && symbolic_objects(getfield(x, i), seen, found)
        end
    end
    return found
end

# A compiled explicit affect records which unknowns and parameters it writes as indices, so the
# callback an integrator carries (and an ahead-of-time compiler serializes) holds no symbolic terms.
@testset "ExplicitAffect holds indices, not symbolic terms" begin
    @variables x(t) = 1.0 v(t) = 0.0
    @parameters g = 9.81 r = 0.9
    @named ball = System(
        [D(x) ~ v, D(v) ~ -g], t, [x, v], [g, r];
        continuous_events = [[t ~ 1.0] => [v ~ -r * Pre(v)]]
    )
    sys = mtkcompile(ball)
    prob = ODEProblem(sys, [], (0.0, 2.0); eval_expression = true, eval_module = @__MODULE__)
    cb = continuous_callback(prob)
    aff = cb.affect!
    @test aff isa ModelingToolkitBase.ExplicitAffect
    @test aff.dvs_to_update == [variable_index(sys, v)]
    @test isempty(aff.ps_to_update)
    @test isempty(symbolic_objects(cb))
    sol = solve(prob, DefaultODEAlgorithm(); abstol = 1.0e-8, reltol = 1.0e-8)
    @test sol.retcode == ReturnCode.Success
    # the event fires at t = 1 and resets the velocity: v(1⁺) = -r v(1⁻) = 0.9 g
    i = findfirst(τ -> abs(τ - 1.0) < 1.0e-8, sol.t)
    @test i !== nothing
    @test sol[v][i + 1] ≈ 0.9 * 9.81 atol = 1.0e-6
    @test sol[v][i] ≈ -9.81 atol = 1.0e-6
end

# Several continuous events become one `VectorContinuousCallback`; its affects are kept in tuples so
# that every affect call is a call of a concrete type, and the events fire as before.
@testset "VectorContinuousCallback affects are typed" begin
    @variables x(t) = 0.0
    @discretes c1(t) = 0.0 c2(t) = 0.0
    @named twoev = System(
        [D(x) ~ 1.0], t, [x], [c1, c2];
        continuous_events = [
            ModelingToolkitBase.SymbolicContinuousCallback([x ~ 0.5], [c1 ~ 1.0]; discrete_parameters = [c1]),
            ModelingToolkitBase.SymbolicContinuousCallback([x ~ 1.5], [c2 ~ 1.0]; discrete_parameters = [c2]),
        ]
    )
    sys = mtkcompile(twoev)
    prob = ODEProblem(sys, [], (0.0, 2.0); eval_expression = true, eval_module = @__MODULE__)
    cb = continuous_callback(prob)
    @test cb isa SciMLBase.VectorContinuousCallback
    @test cb.affect!.affects isa Tuple
    @test isconcretetype(typeof(cb.affect!.affects))
    # the affect for each event is called statically: the dispatcher's optimized code has no call through
    # a value, no generic-function `:call` and no `invokelatest` (checked against a deliberately dynamic one)
    integ = init(prob, DefaultODEAlgorithm())
    affects = cb.affect!.affects
    sig = (typeof(affects), Int, typeof(integ))
    @test has_dynamic_call(dynamic_nth_affect, sig)
    # asserted on the Julia versions `juliac --trim` supports; older versions' optimized code differs
    if VERSION >= v"1.12"
        @test !has_dynamic_call(ModelingToolkitBase.apply_nth_affect, sig)
    end
    @test isempty(symbolic_objects(cb))
    sol = solve(prob, DefaultODEAlgorithm(); abstol = 1.0e-8, reltol = 1.0e-8)
    @test sol.retcode == ReturnCode.Success
    @test sol.ps[c1][end] == 1.0 && sol.ps[c2][end] == 1.0
    @test sol(0.25)[1] ≈ 0.25 atol = 1.0e-8
end
