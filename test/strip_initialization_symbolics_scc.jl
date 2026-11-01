using ModelingToolkit, SciMLBase, OrdinaryDiffEqDefault, Test
using ModelingToolkit: t_nounits as t, D_nounits as D, ModelingToolkitBase
using SymbolicUtils: BasicSymbolic

# The initialization solve follows the integrator's tolerances.
const tols = (; abstol = 1.0e-12, reltol = 1.0e-12)

# Every symbolic term or system reachable from `x`. Modules and types are not entered.
function symbolic_objects(x, seen = IdDict{Any, Nothing}(), found = Any[])
    (x isa Module || x isa Type || x isa Symbol || x isa String || isbits(x)) && return found
    haskey(seen, x) && return found
    seen[x] = nothing
    if x isa BasicSymbolic || x isa ModelingToolkitBase.AbstractSystem
        push!(found, x)
        return found
    end
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

# With x = 1, the initialization solves y^2 = k x and then z^2 = y: two nonlinear SCCs with
# the positive roots y = √k and z = k^(1/4).
@testset "strip_initialization_symbolics, nonlinear SCC initialization" begin
    @parameters k = 2.0
    @variables x(t) y(t) z(t)
    @mtkcompile sys = System([D(x) ~ -x, 0 ~ y^2 - k * x, 0 ~ z^2 - y], t)
    kws = (; guesses = [y => 1.0, z => 1.0], use_scc = true)
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0); kws...)
    lean = ODEProblem(sys, [x => 1.0], (0.0, 1.0); kws..., strip_initialization_symbolics = true)

    full_init = prob.f.initialization_data.initializeprob
    initprob = lean.f.initialization_data.initializeprob
    @test full_init isa SCCNonlinearProblem
    @test initprob isa SCCNonlinearProblem
    @test length(initprob.probs) == length(full_init.probs) == 2
    @test all(p -> p isa NonlinearProblem, initprob.probs)
    @test !isempty(symbolic_objects(full_init))
    @test isempty(symbolic_objects(initprob))
    @test isempty(symbolic_objects(lean.f.initialization_data))
    for p in initprob.probs
        @test p.f.sys === nothing
        @test_throws ErrorException p.f.observed(y, p.u0, p.p)
    end

    integ = init(lean, DefaultODEAlgorithm(); tols...)
    @test integ[x] == 1.0
    @test integ[y] ≈ sqrt(2.0)
    @test integ[z] ≈ 2.0^(1 / 4)
    lean.ps[k] = 3.0
    integ = init(lean, DefaultODEAlgorithm(); tols...)
    @test integ[y] ≈ sqrt(3.0)
    @test integ[z] ≈ 3.0^(1 / 4)
end
