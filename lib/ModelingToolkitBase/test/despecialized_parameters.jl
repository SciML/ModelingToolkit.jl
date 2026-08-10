using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D, MTKParameters
using OrdinaryDiffEq
using Random
using SciMLStructures
using SymbolicIndexingInterface
using Test

@testset "despecialized AutoSpecialize parameters" begin
    @parameters a = 2.0
    @variables x(t) = 1.0 z(t)
    sys = mtkcompile(
        System(
            [D(x) ~ -a * x], t;
            observed = [z ~ a * x], name = :despecialized_parameters
        )
    )

    auto_prob = ODEProblem(sys, [], (0.0, 1.0); jac = true)
    full_prob = ODEProblem{true, SciMLBase.FullSpecialize}(sys, [], (0.0, 1.0))
    expression_prob = eval(ODEProblem(sys, [], (0.0, 1.0); expression = Val{true}))

    @test auto_prob.p isa SciMLBase.DespecializedParameters
    @test auto_prob.p.params isa MTKParameters
    @test auto_prob.p.tunable === auto_prob.p.params.tunable
    @test eltype(typeof(auto_prob.p)) === Any
    @test SciMLStructures.isscimlstructure(auto_prob.p)
    @test SciMLStructures.ismutablescimlstructure(auto_prob.p)
    @test full_prob.p isa MTKParameters
    @test expression_prob.p isa SciMLBase.DespecializedParameters

    @parameters b = 2.0
    @variables y = 0.0
    nonlinear_sys = mtkcompile(
        System([0 ~ y - b]; name = :despecialized_nonlinear_parameters)
    )
    nonlinear_prob = NonlinearProblem(nonlinear_sys, [y => 1.0])
    full_nonlinear_prob = NonlinearProblem{true, SciMLBase.FullSpecialize}(
        nonlinear_sys, [y => 1.0]
    )
    @test nonlinear_prob.p isa SciMLBase.DespecializedParameters
    @test full_nonlinear_prob.p isa MTKParameters
    @test nonlinear_prob.f(nonlinear_prob.u0, nonlinear_prob.p) == [-1.0]

    du = similar(auto_prob.u0)
    auto_prob.f(du, auto_prob.u0, auto_prob.p, first(auto_prob.tspan))
    @test du == [-2.0]
    jac = zeros(1, 1)
    auto_prob.f.jac(jac, auto_prob.u0, auto_prob.p, first(auto_prob.tspan))
    @test jac == [-2.0;;]

    get_a = getp(sys, a)
    set_a! = setp(sys, a)
    @test get_a(auto_prob.p) == 2.0
    set_a!(auto_prob.p, 3.0)
    @test get_a(auto_prob.p) == 3.0

    tunables, repack, aliases = SciMLStructures.canonicalize(
        SciMLStructures.Tunable(), auto_prob.p
    )
    @test aliases
    @test tunables == [3.0]
    new_p = repack([4.0])
    @test new_p isa SciMLBase.DespecializedParameters
    @test get_a(new_p) == 4.0
    @test copy(new_p) == new_p
    @test SciMLStructures.replace!(SciMLStructures.Tunable(), new_p, [5.0]) === nothing
    @test get_a(new_p) == 5.0
    replaced_p = SciMLStructures.replace(SciMLStructures.Tunable(), new_p, [4.0])
    @test replaced_p isa SciMLBase.DespecializedParameters
    @test get_a(replaced_p) == 4.0

    other_inner = SciMLStructures.replace(
        SciMLStructures.Tunable(), new_p.params, Float32[5]
    )
    @test typeof(new_p) === typeof(SciMLBase.DespecializedParameters(other_inner))
    @test typeof(new_p.params) !== typeof(other_inner)
    other_prob = remake(auto_prob; p = SciMLBase.DespecializedParameters(other_inner))
    @test typeof(other_prob) === typeof(auto_prob)

    remade = remake(auto_prob; p = [a => 4.0])
    @test remade.p isa SciMLBase.DespecializedParameters
    sol = solve(remade, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10)
    @test sol[x][end] ≈ exp(-4)
    @test sol[z][end] ≈ 4exp(-4)

    @parameters drift = 0.0 rate = 1.0
    @variables population(t) = 10.0
    jump = SymbolicMassActionJump(rate, [population => 1], [population => -1])
    hybrid_sys = mtkcompile(
        System(
            [D(population) ~ drift], t, [population], [drift, rate]; jumps = [jump],
            name = :despecialized_jump_parameters
        )
    )
    hybrid_prob = JumpProblem(hybrid_sys, [], (0.0, 0.1); rng = Xoshiro(1))
    @test hybrid_prob.prob.p isa SciMLBase.DespecializedParameters
    @test SciMLBase.successful_retcode(solve(hybrid_prob, Tsit5()))
end
