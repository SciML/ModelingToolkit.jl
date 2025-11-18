using ModelingToolkit
using Test
using ModelingToolkit: t_nounits as t, D_nounits as D, SymbolicDiscreteCallback,
                       SymbolicContinuousCallback
using OrdinaryDiffEq
using StochasticDiffEq
using JumpProcesses
using StableRNGs
using SciMLStructures: canonicalize, Tunable, replace, replace!
using SymbolicIndexingInterface
using NonlinearSolve
import DiffEqNoiseProcess

@testset "ODESystem with callbacks" begin
    @discretes p1(t) = 1.0
    @parameters p2 = 2p1
    @variables x(t)
    cb1 = SymbolicContinuousCallback([x ~ 2.0] => [p1 ~ 2.0], discrete_parameters = [p1]) # triggers at t=-2+√6
    function affect1!(mod, obs, ctx, integ)
        return (; p1 = obs.p2)
    end
    cb2 = [x ~ 4.0] => (f = affect1!, observed = (; p2), modified = (; p1)) # triggers at t=-2+√7
    cb3 = SymbolicDiscreteCallback([1.0] => [p1 ~ 5.0], discrete_parameters = [p1])

    @mtkcompile sys = System(
        [D(x) ~ p1 * t + p2],
        t;
        continuous_events = [cb1, cb2],
        discrete_events = [cb3]
    )
    @test !(p2 in Set(parameters(sys)))
    @test p2 in Set(bound_parameters(sys))
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.5), jac = true)
    @test prob.ps[p1] == 1.0
    @test prob.ps[p2] == 2.0
    @test SciMLBase.successful_retcode(solve(prob, Tsit5()))
    prob = ODEProblem(sys, [x => 1.0, p1 => 1.0], (0.0, 1.5), jac = true)
    @test prob.ps[p1] == 1.0
    @test prob.ps[p2] == 2.0
    integ = init(prob, Tsit5())
    @test integ.ps[p1] == 1.0
    @test integ.ps[p2] == 2.0
    step!(integ, 0.5, true) # after cb1, before cb2
    @test integ.ps[p1] == 2.0
    @test integ.ps[p2] == 4.0
    step!(integ, 0.4, true) # after cb2, before cb3
    @test integ.ps[p1] == 4.0
    @test integ.ps[p2] == 8.0
    step!(integ, 0.2, true) # after cb3
    @test integ.ps[p1] == 5.0
    @test integ.ps[p2] == 10.0
end

@testset "vector parameter bindings" begin
    @parameters p1[1:2]=[1.0, 2.0] p2[1:2]=2p1
    @variables x(t) = 0

    @named sys = System(
        [D(x) ~ sum(p1) * t + sum(p2)],
        t
    )
    prob = ODEProblem(complete(sys), [], (0.0, 1.0))
    setp1! = setp(prob, p1)
    get_p1 = getp(prob, p1)
    get_p2 = getp(prob, p2)
    setp1!(prob, [1.5, 2.5])

    @test get_p1(prob) == [1.5, 2.5]
    @test get_p2(prob) == [3.0, 5.0]
end

@testset "extend" begin
    @parameters p1=1.0 p2
    @variables x(t) = 0

    @mtkcompile sys1 = System(
        [D(x) ~ p1 * t + p2],
        t
    )
    @named sys2 = System(Equation[],
        t, [], [];
        bindings = [p2 => 2p1]
    )
    sys = complete(extend(sys2, sys1))
    @test !(p2 in Set(parameters(sys)))
    @test p2 in Set(bound_parameters(sys))
    prob = ODEProblem(complete(sys), nothing, (0.0, 1.0))
    get_dep = getu(prob, 2p2)
    @test get_dep(prob) == 4
end

@testset "getu with parameter bindings" begin
    @parameters p1=1.0 p2=2p1
    @variables x(t) = 0

    @named sys = System(
        [D(x) ~ p1 * t + p2],
        t
    )
    prob = ODEProblem(complete(sys), nothing, (0.0, 1.0))
    get_dep = getu(prob, 2p2)
    @test get_dep(prob) == 4
end

@testset "getu with vector parameter bindings" begin
    @parameters p1[1:2]=[1.0, 2.0] p2[1:2]=2p1
    @variables x(t) = 0

    @named sys = System(
        [D(x) ~ sum(p1) * t + sum(p2)],
        t
    )
    prob = ODEProblem(complete(sys), [], (0.0, 1.0))
    get_dep = getu(prob, 2p1)
    @test get_dep(prob) == [2.0, 4.0]
end

@testset "composing systems with parameter bindings" begin
    @parameters p1=1.0 p2
    @variables x(t) = 0

    @named sys1 = System(
        [D(x) ~ p1 * t + p2],
        t; initial_conditions = [p2 => 1.0]
    )
    @named sys2 = System(
        [D(x) ~ p1 * t - p2],
        t, bindings = [p2 => 2p1]
    )
    sys = complete(System(Equation[], t, systems = [sys1, sys2], name = :sys))

    prob = ODEProblem(sys, [], (0.0, 1.0))
    v1 = sys.sys2.p2
    v2 = 2 * v1
    @test is_observed(prob, v1)
    @test is_observed(prob, v2)
    get_v1 = getu(prob, v1)
    get_v2 = getu(prob, v2)
    @test get_v1(prob) == 2
    @test get_v2(prob) == 4

    setp1! = setp(prob, sys2.p1)
    setp1!(prob, 2.5)
    @test prob.ps[sys2.p2] == 5.0

    new_prob = remake(prob, p = [sys2.p1 => 1.5])

    @test !isempty(ModelingToolkit.bindings(sys))
    @test new_prob.ps[sys2.p1] == 1.5
    @test new_prob.ps[sys2.p2] == 3.0
end

@testset "parameter bindings across model hierarchy" begin
    sys2 = let name = :sys2
        @parameters p2
        @variables x(t) = 1.0
        eqs = [D(x) ~ p2]
        System(eqs, t, [x], [p2]; name)
    end

    @parameters p1 = 1.0
    sys1 = System(
        Equation[], t, [], [p1];
        bindings = [sys2.p2 => 2p1], name = :sys1, systems = [sys2])

    sys = mtkcompile(sys1)

    prob = ODEProblem(sys, [], (0.0, 1.0))
    sol = solve(prob)
    @test SciMLBase.successful_retcode(sol)
end

@testset "Change Tunables" begin
    @variables θ(t)=π/6 ω(t)=0.
    @parameters g=9.81 L=1.0 b=0.1 errp=1
    eqs = [
        D(θ) ~ ω,
        D(ω) ~ -(g/L)*sin(θ) - b*ω
    ]
    @named pendulum_sys = System(eqs, t, [θ, ω], [g, L, b])
    sys = mtkcompile(pendulum_sys)

    new_tunables = [L, b]
    old_tunables = copy(ModelingToolkit.tunable_parameters(sys, ModelingToolkit.parameters(sys)))
    sys2 = ModelingToolkit.subset_tunables(sys, new_tunables)
    sys2_tunables = ModelingToolkit.tunable_parameters(sys2, ModelingToolkit.parameters(sys2))
    @test length(sys2_tunables) == 2
    @test isempty(setdiff(sys2_tunables, new_tunables))
    @test_throws ArgumentError ModelingToolkit.subset_tunables(sys, [errp])
    @test_throws ArgumentError ModelingToolkit.subset_tunables(sys, [θ, L])
    sys3 = ModelingToolkit.subset_tunables(sys, [])
    sys3_tunables = ModelingToolkit.tunable_parameters(sys3, ModelingToolkit.parameters(sys3))
    @test length(sys3_tunables) == 0

    sys_incomplete = pendulum_sys
    @test_throws ArgumentError ModelingToolkit.subset_tunables(sys_incomplete, new_tunables)
    sys_nonsplit = mtkcompile(pendulum_sys; split = false)
    @test_throws ArgumentError ModelingToolkit.subset_tunables(sys_nonsplit, new_tunables)

    @test length(ModelingToolkit.tunable_parameters(sys, ModelingToolkit.parameters(sys))) == length(old_tunables)
end

struct CallableFoo
    p::Any
end

@register_symbolic CallableFoo(x)

(f::CallableFoo)(x) = f.p + x

@testset "callable parameters" begin
    @variables y(t) = 1
    @parameters p=2 (i::CallableFoo)(..)

    eqs = [D(y) ~ i(t) + p]
    @named model = System(eqs, t, [y], [p, i]; bindings = [i => CallableFoo(p)])
    sys = mtkcompile(model)

    prob = ODEProblem(sys, [], (0.0, 1.0))
    sol = solve(prob, Tsit5())

    @test SciMLBase.successful_retcode(sol)
end

@testset "SDESystem" begin
    @parameters σ ρ β
    @variables x(t) y(t) z(t)

    eqs = [D(x) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z]

    noiseeqs = [0.1 * x,
        0.1 * y,
        0.1 * z]

    @named sys = System(eqs, t)
    @named sdesys = SDESystem(sys, noiseeqs; bindings = [ρ => 2σ])
    sdesys = complete(sdesys)
    @test !(ρ in Set(parameters(sdesys)))
    @test ρ in Set(bound_parameters(sdesys))

    prob = SDEProblem(
        sdesys, [x => 1.0, y => 0.0, z => 0.0, σ => 10.0, β => 2.33], (0.0, 100.0))
    @test prob.ps[ρ] == 2prob.ps[σ]
    @test_nowarn solve(prob, SRIW1())
end

@testset "JumpSystem" begin
    rng = StableRNG(12345)
    @parameters β γ
    @constants h = 1
    @variables S(t) I(t) R(t)
    rate₁ = β * S * I * h
    affect₁ = [S ~ Pre(S) - 1 * h, I ~ Pre(I) + 1]
    rate₃ = γ * I * h
    affect₃ = [I ~ Pre(I) * h - 1, R ~ Pre(R) + 1]
    j₁ = ConstantRateJump(rate₁, affect₁)
    j₃ = ConstantRateJump(rate₃, affect₃)
    @named js2 = JumpSystem(
        [j₃], t, [S, I, R], [β, γ, h]; bindings = [β => 0.01γ])
    @test issetequal(parameters(js2), [β, γ, h])
    js2 = complete(js2)
    @test issetequal(parameters(js2), [γ, h])
    @test Set(bound_parameters(js2)) == Set([β])
    tspan = (0.0, 250.0)
    u₀map = [S => 999, I => 1, R => 0]
    parammap = [γ => 0.01]
    jprob = JumpProblem(js2, [u₀map; parammap], tspan; aggregator = Direct(),
        save_positions = (false, false), rng = rng)
    @test jprob.ps[γ] == 0.01
    @test jprob.ps[β] == 0.0001
    @test_nowarn solve(jprob, SSAStepper())
end

@testset "NonlinearSystem" begin
    @parameters p1=1.0 p2
    @variables x(t)
    eqs = [0 ~ p1 * x * exp(x) + p2]
    @mtkcompile sys = System(eqs; bindings = [p2 => 2p1])
    @test isequal(only(parameters(sys)), p1)
    @test Set(bound_parameters(sys)) == Set([p2])
    prob = NonlinearProblem(sys, [x => 1.0])
    @test prob.ps[p1] == 1.0
    @test prob.ps[p2] == 2.0
    @test_nowarn solve(prob, NewtonRaphson())
    prob = NonlinearProblem(sys, [x => 1.0, p1 => 2.0])
    @test prob.ps[p1] == 2.0
    @test prob.ps[p2] == 4.0
end

@testset "SciMLStructures interface" begin
    @parameters p1=1.0 p2
    @variables x(t)
    cb1 = [x ~ 2.0] => [p1 ~ 2.0] # triggers at t=-2+√6
    function affect1!(integ, u, p, ctx)
        integ.ps[p[1]] = integ.ps[p[2]]
    end
    cb2 = [x ~ 4.0] => (affect1!, [], [p1, p2], [p1]) # triggers at t=-2+√7
    cb3 = [1.0] => [p1 ~ 5.0]

    @mtkcompile sys = System(
        [D(x) ~ p1 * t + p2],
        t; bindings = [p2 => 2p1]
    )
    prob = ODEProblem(sys, [x => 1.0, p1 => 1.0], (0.0, 1.5), jac = true)
    prob.ps[p1] = 3.0
    @test prob.ps[p1] == 3.0
    @test prob.ps[p2] == 6.0

    ps = prob.p
    buffer, repack, _ = canonicalize(Tunable(), ps)
    idx = parameter_index(sys, p1)
    @test buffer[idx.idx] == 3.0
    buffer[idx.idx] = 4.0
    ps = repack(buffer)
    @test getp(sys, p1)(ps) == 4.0
    @test getp(sys, p2)(ps) == 8.0

    replace!(Tunable(), ps, ones(length(ps.tunable)))
    @test getp(sys, p1)(ps) == 1.0
    @test getp(sys, p2)(ps) == 2.0

    ps2 = replace(Tunable(), ps, 2 .* ps.tunable)
    @test getp(sys, p1)(ps2) == 2.0
    @test getp(sys, p2)(ps2) == 4.0
end

@testset "Scalarized array as RHS of parameter binding" begin
    @parameters p[1:2] p1 = p[1] p2 = p[2]
    @variables x(t)
    @named sys = System([D(x) ~ x], t, [x], [p, p1, p2])
    sys = complete(sys)
    @test any(isequal(p), ModelingToolkit.get_ps(sys))
    sys = mtkcompile(sys)
    @test length(ModelingToolkit.bound_parameters(sys)) == 2
end
