using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D, MTKParameters
using DiffEqBase
using OrdinaryDiffEq
using Random
using SciMLStructures
using SymbolicIndexingInterface
using Test

@testset "AutoDespecialize parameters" begin
    @parameters a = 2.0
    @variables x(t) = 1.0 z(t)
    sys = mtkcompile(
        System(
            [D(x) ~ -a * x], t;
            observed = [z ~ a * x], name = :despecialized_parameters
        )
    )

    despecialized_prob = ODEProblem(sys, [], (0.0, 1.0); jac = true)
    auto_prob = ODEProblem{true, SciMLBase.AutoSpecialize}(sys, [], (0.0, 1.0))
    respecialized_prob = ODEProblem{true, SciMLBase.AutoRespecialize}(
        sys, [], (0.0, 1.0)
    )
    full_prob = ODEProblem{true, SciMLBase.FullSpecialize}(sys, [], (0.0, 1.0))
    function_wrapper_prob = ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(
        sys, [], (0.0, 1.0)
    )
    expression_prob = eval(ODEProblem(sys, [], (0.0, 1.0); expression = Val{true}))

    @test SciMLBase.specialization(despecialized_prob.f) === SciMLBase.AutoDespecialize
    @test SciMLBase.specialization(auto_prob.f) === SciMLBase.AutoSpecialize
    @test despecialized_prob.p isa MTKParameters
    @test auto_prob.p isa MTKParameters
    @test respecialized_prob.p isa MTKParameters
    @test full_prob.p isa MTKParameters
    @test despecialized_prob.f.initialization_data.initializeprob.p isa MTKParameters
    @test respecialized_prob.f.initialization_data.initializeprob.p isa MTKParameters
    @test full_prob.f.initialization_data.initializeprob.p isa MTKParameters
    @test SciMLBase.specialization(
        despecialized_prob.f.initialization_data.initializeprob.f
    ) === SciMLBase.AutoDespecialize
    @test SciMLBase.specialization(auto_prob.f.initialization_data.initializeprob.f) ===
        SciMLBase.AutoSpecialize
    @test SciMLBase.specialization(
        respecialized_prob.f.initialization_data.initializeprob.f
    ) === SciMLBase.AutoSpecialize
    @test SciMLBase.specialization(full_prob.f.initialization_data.initializeprob.f) ===
        SciMLBase.AutoSpecialize
    @test SciMLBase.specialization(
        function_wrapper_prob.f.initialization_data.initializeprob.f
    ) === SciMLBase.AutoSpecialize
    @test expression_prob.p isa MTKParameters

    concrete_prob = DiffEqBase.get_concrete_problem(
        despecialized_prob, true; alg = Tsit5()
    )
    concrete_parameters = SciMLBase.unwrap_parameters(concrete_prob.p)
    @test concrete_prob.p isa SciMLBase.DespecializedParameters
    @test concrete_parameters isa MTKParameters
    @test concrete_prob.p.tunable === concrete_parameters.tunable
    @test eltype(typeof(concrete_prob.p)) === Any
    @test SciMLStructures.isscimlstructure(concrete_prob.p)
    @test SciMLStructures.ismutablescimlstructure(concrete_prob.p)

    seen_parameter_type = Ref{DataType}()
    gfw = ModelingToolkitBase.GeneratedFunctionWrapper{Tuple{2, 3, false}}(
        (u, p, t) -> (seen_parameter_type[] = typeof(p); u),
        (du, u, p, t) -> (seen_parameter_type[] = typeof(p); copyto!(du, u))
    )
    @test gfw(concrete_prob.u0, concrete_prob.p, 0.0) === concrete_prob.u0
    @test seen_parameter_type[] === typeof(concrete_parameters)
    du = similar(concrete_prob.u0)
    gfw(du, concrete_prob.u0, concrete_prob.p, 0.0)
    @test seen_parameter_type[] === typeof(concrete_parameters)

    @parameters b = 2.0
    @variables y = 0.0
    nonlinear_sys = mtkcompile(
        System([0 ~ y - b]; name = :despecialized_nonlinear_parameters)
    )
    nonlinear_prob = NonlinearProblem(nonlinear_sys, [y => 1.0])
    full_nonlinear_prob = NonlinearProblem{true, SciMLBase.FullSpecialize}(
        nonlinear_sys, [y => 1.0]
    )
    @test SciMLBase.specialization(nonlinear_prob.f) === SciMLBase.AutoDespecialize
    @test nonlinear_prob.p isa MTKParameters
    @test full_nonlinear_prob.p isa MTKParameters
    @test nonlinear_prob.f(nonlinear_prob.u0, nonlinear_prob.p) == [-1.0]

    du = similar(concrete_prob.u0)
    concrete_prob.f(du, concrete_prob.u0, concrete_prob.p, first(concrete_prob.tspan))
    @test du == [-2.0]
    jac = zeros(1, 1)
    concrete_prob.f.jac(
        jac, concrete_prob.u0, concrete_prob.p, first(concrete_prob.tspan)
    )
    @test jac == [-2.0;;]

    get_a = getp(sys, a)
    set_a! = setp(sys, a)
    @test get_a(concrete_prob.p) == 2.0
    set_a!(concrete_prob.p, 3.0)
    @test get_a(concrete_prob.p) == 3.0

    tunables, repack, aliases = SciMLStructures.canonicalize(
        SciMLStructures.Tunable(), concrete_prob.p
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
        SciMLStructures.Tunable(), SciMLBase.unwrap_parameters(new_p), Float32[5]
    )
    @test typeof(new_p) === typeof(SciMLBase.DespecializedParameters(other_inner))
    @test typeof(SciMLBase.unwrap_parameters(new_p)) !== typeof(other_inner)
    other_prob = remake(concrete_prob; p = SciMLBase.DespecializedParameters(other_inner))
    @test typeof(other_prob) === typeof(concrete_prob)

    remade = remake(despecialized_prob; p = [a => 4.0])
    @test remade.p isa MTKParameters
    sol = solve(remade, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10)
    @test sol.prob.p isa SciMLBase.DespecializedParameters
    @test sol[x][end] ≈ exp(-4)
    @test sol[z][end] ≈ 4exp(-4)

    @variables event_x(t) = 0.0
    event_sys = complete(
        System(
            [D(event_x) ~ 1], t; continuous_events = [event_x ~ 1],
            name = :despecialized_nonsplit_event
        ); split = false
    )
    event_prob = ODEProblem(event_sys, [], (0.0, 2.0))
    concrete_event_prob = DiffEqBase.get_concrete_problem(event_prob, true; alg = Tsit5())
    @test concrete_event_prob.p isa SciMLBase.DespecializedParameters
    event_sol = solve(event_prob, Tsit5())
    @test event_sol.prob.p isa SciMLBase.DespecializedParameters
    @test minimum(t -> abs(t - 1), event_sol.t) < 1.0e-10

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
    @test hybrid_prob.prob.p isa MTKParameters
    @test SciMLBase.successful_retcode(solve(hybrid_prob, Tsit5()))
end
