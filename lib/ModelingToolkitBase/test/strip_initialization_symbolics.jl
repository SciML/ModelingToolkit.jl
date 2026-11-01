using ModelingToolkitBase, SciMLBase, OrdinaryDiffEqDefault, NonlinearSolve, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D
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

function test_stripped_metadata(idata)
    @test isempty(idata.metadata.op)
    @test isempty(idata.metadata.guesses)
    @test isempty(idata.metadata.additional_initialization_eqs)
    return nothing
end

# With x = 1 the initialization is the linear y = k x.
@testset "strip_initialization_symbolics, linear initialization" begin
    @parameters k = 2.0
    @variables x(t) y(t)
    eqs = [D(x) ~ -k * x, 0 ~ y - k * x]
    @mtkcompile sys = System(eqs, t)
    op = [x => 1.0]
    prob = ODEProblem(sys, op, (0.0, 1.0))
    lean = ODEProblem(sys, op, (0.0, 1.0); strip_initialization_symbolics = true)
    idata = lean.f.initialization_data
    initprob = idata.initializeprob
    @test initprob isa SCCNonlinearProblem
    @test only(initprob.probs) isa LinearProblem
    @test initprob.f.sys === nothing
    @test_throws ErrorException initprob.f.observed(y, initprob.u0, initprob.p)
    @test only(initprob.probs).f.sys isa ModelingToolkitBase.StrippedLinearSystem
    @test !isempty(symbolic_objects(prob.f.initialization_data.initializeprob))
    @test isempty(symbolic_objects(initprob))
    @test isempty(symbolic_objects(lean.f.initialization_data))
    test_stripped_metadata(idata)

    integ = init(lean, DefaultODEAlgorithm(); tols...)
    @test integ[x] == 1.0
    @test integ[y] ≈ 2.0
    @test integ.u == init(prob, DefaultODEAlgorithm()).u
    lean.ps[k] = 3.0
    integ = init(lean, DefaultODEAlgorithm(); tols...)
    @test integ[x] == 1.0
    @test integ[y] ≈ 3.0
end

# With x = 1 the initialization is y^2 = k x, whose positive root from the guess y = 1 is √k.
@testset "strip_initialization_symbolics, nonlinear initialization" begin
    @parameters k = 2.0
    @variables x(t) y(t)
    @mtkcompile sys = System([D(x) ~ -x, 0 ~ y^2 - k * x], t)
    kws = (; guesses = [y => 1.0], use_scc = false)
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0); kws...)
    lean = ODEProblem(sys, [x => 1.0], (0.0, 1.0); kws..., strip_initialization_symbolics = true)
    idata = lean.f.initialization_data
    initprob = idata.initializeprob
    @test initprob isa NonlinearProblem
    @test initprob.f.sys === nothing
    @test_throws ErrorException initprob.f.observed(y, initprob.u0, initprob.p)
    @test !isempty(symbolic_objects(prob.f.initialization_data.initializeprob))
    @test isempty(symbolic_objects(initprob))
    @test isempty(symbolic_objects(lean.f.initialization_data))
    test_stripped_metadata(idata)

    integ = init(lean, DefaultODEAlgorithm(); tols...)
    @test integ[x] == 1.0
    @test integ[y] ≈ sqrt(2.0)
    lean.ps[k] = 3.0
    integ = init(lean, DefaultODEAlgorithm(); tols...)
    @test integ[y] ≈ sqrt(3.0)
end

# Three initialization equations in two unknowns, consistent at x = y = 2 (y^2 = k x with
# k = 2), make an overdetermined `NonlinearLeastSquaresProblem`.
@testset "strip_initialization_symbolics, least-squares initialization" begin
    @parameters k = 2.0
    @variables x(t) y(t)
    @mtkcompile sys = System([D(x) ~ -x, 0 ~ y^2 - k * x], t)
    kws = (; initialization_eqs = [x ~ 2.0, y ~ 2.0], guesses = [x => 1.0, y => 1.0])
    prob = ODEProblem(sys, [], (0.0, 1.0); kws...)
    lean = ODEProblem(sys, [], (0.0, 1.0); kws..., strip_initialization_symbolics = true)
    idata = lean.f.initialization_data
    initprob = idata.initializeprob
    @test initprob isa NonlinearLeastSquaresProblem
    @test initprob.f.sys === nothing
    @test !isempty(symbolic_objects(prob.f.initialization_data.initializeprob))
    @test isempty(symbolic_objects(initprob))
    @test isempty(symbolic_objects(lean.f.initialization_data))
    test_stripped_metadata(idata)

    integ = init(lean, DefaultODEAlgorithm(); tols...)
    @test integ[x] ≈ 2.0
    @test integ[y] ≈ 2.0
end

# With x = 1 the initialization is atan(y - 3) = 0, reached from the out-of-basin guess
# y = 12 by continuation from the simplified equation y = 0; its root is y = 3.
@testset "strip_initialization_symbolics, homotopy initialization" begin
    @variables x(t) y(t)
    @mtkcompile sys = System(
        [D(x) ~ -x, 0 ~ ModelingToolkitBase.homotopy(atan(y - 3), y)], t;
        guesses = [y => 12.0]
    )
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0))
    lean = ODEProblem(sys, [x => 1.0], (0.0, 1.0); strip_initialization_symbolics = true)
    full_init = prob.f.initialization_data.initializeprob
    idata = lean.f.initialization_data
    initprob = idata.initializeprob
    @test initprob isa SciMLBase.HomotopyProblem
    @test initprob.λspan == full_init.λspan
    @test initprob.u0 == full_init.u0
    @test initprob.f.sys === nothing
    @test_throws ErrorException initprob.f.observed(y, initprob.u0, initprob.p)
    @test !isempty(symbolic_objects(full_init))
    @test isempty(symbolic_objects(initprob))
    @test isempty(symbolic_objects(lean.f.initialization_data))
    test_stripped_metadata(idata)

    integ = init(lean, DefaultODEAlgorithm(); tols...)
    @test integ[x] == 1.0
    @test integ[y] ≈ 3.0
end

# The helper keeps the problem's numerical data and keywords; the residual u - p, 2(u - p) has
# its least-squares solution at u = p.
@testset "strip_problem_symbolics(::NonlinearLeastSquaresProblem)" begin
    @variables u
    @named nsys = System([0 ~ u - 1], [u], [])
    f = NonlinearFunction(
        (r, u, p) -> (r[1] = u[1] - p; r[2] = 2 * (u[1] - p); nothing);
        sys = nsys, resid_prototype = zeros(2)
    )
    q = NonlinearLeastSquaresProblem(f, [0.0], 2.0; abstol = 1.0e-12)
    r = ModelingToolkitBase.strip_problem_symbolics(q)
    @test r isa NonlinearLeastSquaresProblem
    @test SciMLBase.isinplace(r)
    @test r.f.sys === nothing
    @test r.f.resid_prototype == zeros(2)
    @test r.u0 == [0.0]
    @test r.p == 2.0
    @test r.kwargs[:abstol] == 1.0e-12
    @test solve(r).u ≈ [2.0]
end

# A numeric `u0` passed to `remake` is written into the `Initial` parameters. With
# x = 4 the initialization is y^2 = k x with k = 2, whose positive root from y = 1 is √8.
@testset "strip_initialization_symbolics, numeric remake, $S, split = $split" for S in (SciMLBase.FullSpecialize, SciMLBase.AutoSpecialize),
        split in (true, false)

    @parameters k = 2.0
    @variables x(t) y(t)
    @mtkcompile sys = System([D(x) ~ -x, 0 ~ y^2 - k * x], t) split = split
    kws = (; guesses = [y => 1.0], use_scc = false)
    prob = ODEProblem{true, S}(sys, [x => 1.0], (0.0, 1.0); kws...)
    lean = ODEProblem{true, S}(
        sys, [x => 1.0], (0.0, 1.0); kws..., strip_initialization_symbolics = true
    )
    @test !isempty(symbolic_objects(prob.f.initialization_data))
    @test isempty(symbolic_objects(lean.f.initialization_data))

    u0 = map(v -> isequal(v, x) ? 4.0 : 1.0, unknowns(sys))
    lean2 = remake(lean; u0)
    integ = init(lean2, DefaultODEAlgorithm(); tols...)
    @test integ[x] == 4.0
    @test integ[y] ≈ sqrt(8.0)
end
