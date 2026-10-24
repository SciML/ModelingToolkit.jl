using ModelingToolkitBase, SciMLBase, OrdinaryDiffEqDefault, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@testset "strip_initialization_symbolics" begin
    @parameters k = 2.0
    @variables x(t) y(t)
    eqs = [D(x) ~ -k * x, 0 ~ y - k * x]
    @mtkcompile sys = System(eqs, t)
    op = [x => 1.0]
    prob = ODEProblem(sys, op, (0.0, 1.0))
    lean = ODEProblem(sys, op, (0.0, 1.0); strip_initialization_symbolics = true)
    idata = lean.f.initialization_data
    @test idata.initializeprob.f.sys === nothing
    @test idata.initializeprob.f.observed === SciMLBase.DEFAULT_OBSERVED_NO_TIME
    @test isempty(idata.metadata.op) && isempty(idata.metadata.guesses)
    @test isempty(idata.metadata.additional_initialization_eqs)
    i1 = init(prob, DefaultODEAlgorithm())
    i2 = init(lean, DefaultODEAlgorithm())
    @test i1.u == i2.u
    prob.ps[k] = 3.0
    lean.ps[k] = 3.0
    @test init(prob, DefaultODEAlgorithm()).u == init(lean, DefaultODEAlgorithm()).u == [1.0, 3.0]
end

@testset "strip_initialization_symbolics, nonlinear initialization" begin
    @parameters k = 2.0
    @variables x(t) y(t)
    @mtkcompile sys = System([D(x) ~ -x, 0 ~ y^2 - k * x], t)
    kws = (; guesses = [y => 1.0], use_scc = false)
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0); kws...)
    lean = ODEProblem(sys, [x => 1.0], (0.0, 1.0); kws..., strip_initialization_symbolics = true)
    initprob = lean.f.initialization_data.initializeprob
    @test initprob isa NonlinearProblem
    @test initprob.f.sys === nothing
    @test init(prob, DefaultODEAlgorithm()).u ≈ init(lean, DefaultODEAlgorithm()).u
end
