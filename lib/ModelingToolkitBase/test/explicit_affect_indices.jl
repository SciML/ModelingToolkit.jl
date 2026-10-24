using ModelingToolkitBase, SciMLBase, OrdinaryDiffEqDefault, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using SymbolicUtils: BasicSymbolic
using SymbolicIndexingInterface: variable_index

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
    @test isempty(symbolic_objects(cb))
    sol = solve(prob, DefaultODEAlgorithm(); abstol = 1.0e-8, reltol = 1.0e-8)
    @test sol.retcode == ReturnCode.Success
    @test sol.ps[c1][end] == 1.0 && sol.ps[c2][end] == 1.0
    @test sol(0.25)[1] ≈ 0.25 atol = 1.0e-8
end
