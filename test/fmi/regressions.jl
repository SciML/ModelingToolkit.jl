using Test, ModelingToolkit, FMI, FMIZoo, OrdinaryDiffEq, SciMLBase
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit as MTK
import ModelingToolkitBase as MTKBase
import FiniteDiff, ForwardDiff

const FIXTURE_FMUS = joinpath(@__DIR__, "fmus")

function feedthrough_step_hook(step)
    # Feedthrough has no state to advance. Supply a communication clock for arbitrary intervals.
    function do_step(instance, point, h, no_reset, event, terminate, early, last_time)
        unsafe_store!(event, UInt8(0))
        unsafe_store!(terminate, UInt8(0))
        unsafe_store!(early, UInt8(0))
        unsafe_store!(last_time, step(point, h))
        return Cuint(FMI.fmi3StatusOK)
    end
    return @cfunction(
        $do_step, Cuint,
        (Ptr{Cvoid}, Cdouble, Cdouble, UInt8, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{Cdouble})
    )
end

function adder_system(fmu, type; initialization_eq = false, kwargs...)
    step = type == :CS ? (; communication_step_size = 0.01) : (;)
    @named adder = MTK.FMIComponent(Val(2); fmu, type, step..., kwargs...)
    @variables z(t) = 0.0
    @mtkcompile sys = System(
        [adder.a ~ 1.0, adder.b ~ 2.0, D(z) ~ 1.0], t; systems = [adder],
        initialization_eqs = initialization_eq ? [adder.c ~ 3.0] : []
    )
    return sys
end

@testset "Co-Simulation holds the left input sample" begin
    for jump in (0.01, 0.025)
        fmu = loadFMU(joinpath(FIXTURE_FMUS, "SimpleAdder.fmu"); type = :CS)
        @named adder = MTK.FMIComponent(Val(2); fmu, type = :CS, communication_step_size = 0.01)
        @variables z(t) = 0.0
        @mtkcompile sys = System(
            [adder.a ~ ifelse(t >= jump, 100.0, 0.0), adder.b ~ 2.0, D(z) ~ 1.0], t;
            systems = [adder]
        )
        sol = solve(ODEProblem(sys, [sys.adder.c => 0.0], (0.0, 0.025)), Tsit5())
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test sol.t[end] == 0.025
        @test sol[sys.adder.c][end] ≈ 3 * 0.025 + 100 * (0.025 - jump)
    end
end

@testset "FMI 3 reports the achieved Co-Simulation time" begin
    for achieved in (0.0025, 0.005)
        fmu = loadFMU("Feedthrough", "ModelicaReferenceFMUs", "0.0.39", "3.0"; type = :CS)
        hook = feedthrough_step_hook((point, h) -> achieved)
        fmu.cDoStep = Base.unsafe_convert(Ptr{Cvoid}, hook)
        GC.@preserve hook begin
            @named ports = MTK.FMIComponent(Val(3); fmu, type = :CS, communication_step_size = 0.01)
            eqs = [
                getproperty(ports, name) ~ (name == :Enumeration_input ? 1.0 : 0.0)
                    for name in (MTK.getname(v) for v in unknowns(ports) if MTK.isinput(v))
            ]
            @variables z(t) = 0.0
            @mtkcompile sys = System([eqs; D(z) ~ 1.0], t; systems = [ports])
            prob = ODEProblem(sys, [], (0.0, 0.005))
            if achieved == 0.005
                sol = solve(prob, Tsit5())
                @test sol.retcode == SciMLBase.ReturnCode.Success
                @test sol.t[end] == achieved
            else
                @test_throws r"FMI 3.*0\.005.*0\.0025" solve(prob, Tsit5())
            end
        end
    end
end

@testset "Explicit FMU initialization equations" begin
    @testset "$type, $source" for type in (:ME, :CS), source in (:system, :problem)
        fmu = loadFMU(joinpath(FIXTURE_FMUS, "SimpleAdder.fmu"); type)
        sys = adder_system(fmu, type; initialization_eq = source == :system)
        eqs = source == :problem ? [sys.adder.c ~ 3.0] : []
        sol = solve(
            ODEProblem(sys, Dict(), (0.0, 0.1); fully_determined = true, initialization_eqs = eqs), Tsit5()
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol[sys.adder.c][1] == 3.0
    end
end

@testset "FMI 2 references in different base types" begin
    fmu = loadFMU("Feedthrough", "ModelicaReferenceFMUs", "0.0.39", "2.0"; type = :ME)
    md = fmu.modelDescription
    ri, bi = map(("Float64_continuous_input", "Boolean_input")) do name
        only(filter(mv -> mv.name == name, md.modelVariables))
    end
    real_ref, bool_ref = ri.valueReference, bi.valueReference
    bi.valueReference = real_ref
    md.stringValueReferences[bi.name] = real_ref
    replace!(md.inputValueReferences, bool_ref => real_ref)
    original = fmu.cSetBoolean
    function set_boolean(instance, refs, n, values)
        refs = [ref == real_ref ? bool_ref : ref for ref in unsafe_wrap(Vector{UInt32}, refs, Int(n))]
        return ccall(original, Cuint, (Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{Cint}), instance, refs, n, values)
    end
    hook = @cfunction($set_boolean, Cuint, (Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{Cint}))
    fmu.cSetBoolean = Base.unsafe_convert(Ptr{Cvoid}, hook)
    GC.@preserve hook begin
        @named ports = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @test MTK.isinput(ports.Float64_continuous_input)
        @test MTK.isinput(ports.Boolean_input)
        values = Dict(:Float64_continuous_input => 1.5, :Boolean_input => Float64(false), :Enumeration_input => 1.0)
        eqs = [
            getproperty(ports, name) ~ get(values, name, 0.0)
                for name in (MTK.getname(v) for v in unknowns(ports) if MTK.isinput(v))
        ]
        @mtkcompile sys = System(eqs, t; systems = [ports])
        sol = solve(ODEProblem(sys, [], (0.0, 0.02)), Tsit5())
        @test SciMLBase.successful_retcode(sol)
        @test sol[sys.ports.Float64_continuous_output][end] == 1.5
        @test sol[sys.ports.Boolean_output][end] == 0.0
    end
end

@testset "FMI 3 wide integer ports and range diagnostics" begin
    @testset "$type" for type in (:ME, :CS)
        fmu = loadFMU("Feedthrough", "ModelicaReferenceFMUs", "0.0.39", "3.0"; type)
        step = type == :CS ? (; communication_step_size = 0.1) : (;)
        @named ports = MTK.FMIComponent(Val(3); fmu, type, step...)
        @test MTK.isinput(ports.Int64_input)
        @test MTK.isinput(ports.UInt64_input)
        @variables z(t) = 0.0
        function port_problem(i64, u64)
            values = Dict{Symbol, Integer}(
                :Int64_input => i64, :UInt64_input => u64, :Int32_input => typemax(Int32),
                :UInt32_input => typemax(UInt32), :Enumeration_input => 1
            )
            eqs = [
                getproperty(ports, name) ~ get(values, name, 0)
                    for name in (MTK.getname(v) for v in unknowns(ports) if MTK.isinput(v))
            ]
            @mtkcompile sys = System([eqs; D(z) ~ 1.0], t; systems = [ports])
            return sys, ODEProblem(sys, [], (0.0, 0.2))
        end
        sys, prob = port_problem(-137, UInt64(251))
        sol = solve(prob, Tsit5())
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test sol[sys.ports.Int64_output][end] == -137
        @test sol[sys.ports.UInt64_output][end] == 251
        @test sol[sys.ports.Int32_output][end] == typemax(Int32)
        @test sol[sys.ports.UInt32_output][end] == typemax(UInt32)
        @test_throws r"Int64.*Int64_input.*9007199254740993" solve(last(port_problem(9007199254740993, 0)), Tsit5())
        @test_throws r"UInt64.*UInt64_input.*9007199254740993" solve(last(port_problem(0, UInt64(9007199254740993))), Tsit5())
    end
end

@testset "Model Exchange snapshots event-dependent continuous outputs" begin
    fmu = loadFMU(joinpath(FIXTURE_FMUS, "SimpleAdder.fmu"); type = :ME)
    completed = Ref(0)
    original_get = fmu.cGetReal
    original_enter = fmu.cEnterEventMode
    output_ref = fmu.modelDescription.stringValueReferences["out2"]
    offsets = Dict{Ptr{Cvoid}, Float64}()
    function complete(instance, no_reset, event, terminate)
        completed[] += 1
        first_step = completed[] == 1
        unsafe_store!(event, Cint(first_step))
        unsafe_store!(terminate, Cint(0))
        return zero(Cuint)
    end
    function enter(instance)
        offsets[instance] = 100.0
        return ccall(original_enter, Cuint, (Ptr{Cvoid},), instance)
    end
    function get_real(instance, refs, n, values)
        status = ccall(original_get, Cuint, (Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{Cdouble}), instance, refs, n, values)
        for i in 1:Int(n)
            if unsafe_load(refs, i) == output_ref
                unsafe_store!(values, unsafe_load(values, i) + get(offsets, instance, 0.0), i)
            end
        end
        return status
    end
    complete_hook = @cfunction($complete, Cuint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}))
    enter_hook = @cfunction($enter, Cuint, (Ptr{Cvoid},))
    get_hook = @cfunction($get_real, Cuint, (Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{Cdouble}))
    fmu.cCompletedIntegratorStep = Base.unsafe_convert(Ptr{Cvoid}, complete_hook)
    fmu.cEnterEventMode = Base.unsafe_convert(Ptr{Cvoid}, enter_hook)
    fmu.cGetReal = Base.unsafe_convert(Ptr{Cvoid}, get_hook)
    GC.@preserve complete_hook enter_hook get_hook begin
        sys = adder_system(fmu, :ME)
        prob = ODEProblem(sys, [sys.adder.c => 0.0], (0.0, 0.05))
        sol = solve(prob, Tsit5(); dt = 0.01, adaptive = false, saveat = [0.0, 0.005, 0.025, 0.05])
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test sol(0.005; idxs = sys.adder.out2) ≈ 0.04
        @test sol(0.01; idxs = sys.adder.out2, continuity = :left) ≈ 0.08
        @test sol(0.01; idxs = sys.adder.out2, continuity = :right) ≈ 100.08
        @test sol(0.025; idxs = sys.adder.out2) ≈ 100.2
        @test sol(0.05; idxs = sys.adder.out2) ≈ 100.4
        @test_throws r"FMI.*snapshot.*0\.035" sol(0.035; idxs = sys.adder.out2)
        empty!(offsets)
        again = solve(remake(prob; u0 = [sys.adder.c => 2.0]), Tsit5(); dt = 0.01, adaptive = false)
        @test again.retcode == SciMLBase.ReturnCode.Success
        @test again(0.05; idxs = sys.adder.out2) ≈ 4.4
        @test sol(0.05; idxs = sys.adder.out2) ≈ 100.4
    end
end

@testset "Co-Simulation communication reaches the solve endpoint" begin
    @testset "FMI 2, $start to $finish" for (start, finish) in ((0.0, 0.025), (0.0, 1.0e-9), (1.0e8, 1.0e8 + 0.025))
        fmu = loadFMU(joinpath(FIXTURE_FMUS, "SimpleAdder.fmu"); type = :CS)
        sys = adder_system(fmu, :CS; stop_time = start + 0.05)
        sol = solve(ODEProblem(sys, [sys.adder.c => 0.0], (start, finish)), Tsit5())
        @test SciMLBase.successful_retcode(sol)
        @test sol[sys.adder.c][end] ≈ 4 * (finish - start)
    end
    @testset "FMI 3, $start to $finish, variable steps: $variable" for (start, finish, variable) in (
            (0.0, 0.025, true), (0.0, 1.0e-9, true), (1.0e8, 1.0e8 + 0.025, true), (0.0, 0.02, false),
        )
        fmu = loadFMU("Feedthrough", "ModelicaReferenceFMUs", "0.0.39", "3.0"; type = :CS)
        fmu.modelDescription.coSimulation.canHandleVariableCommunicationStepSize = variable
        requests = Tuple{Float64, Float64}[]
        hook = feedthrough_step_hook() do point, h
            push!(requests, (point, h))
            return point + h
        end
        fmu.cDoStep = Base.unsafe_convert(Ptr{Cvoid}, hook)
        GC.@preserve hook begin
            @named ports = MTK.FMIComponent(
                Val(3); fmu, type = :CS, communication_step_size = 0.01, stop_time = start + 0.05
            )
            values = Dict(:Float64_continuous_input => 4 * (t - start), :Enumeration_input => 1.0)
            eqs = [
                getproperty(ports, name) ~ get(values, name, 0.0)
                    for name in (MTK.getname(v) for v in unknowns(ports) if MTK.isinput(v))
            ]
            @mtkcompile sys = System(eqs, t; systems = [ports])
            prob = ODEProblem(sys, [], (start, finish))
            sol = solve(prob, Tsit5())
            @test SciMLBase.successful_retcode(sol)
            @test sol[sys.ports.Float64_continuous_output][end] ≈ 4 * (first(last(requests)) - start)
            @test sum(last, requests; init = 0.0) ≈ finish - start
            if !variable
                @test requests == [(0.0, 0.01), (0.01, 0.01)]
                offgrid = solve(remake(prob; tspan = (0.0, 0.025)), Tsit5())
                @test SciMLBase.successful_retcode(offgrid)
            end
        end
    end
end

@testset "Model Exchange accepted-step event iteration" begin
    @testset "$cause" for cause in (:step, :input, :terminate)
        fmu = loadFMU(joinpath(FIXTURE_FMUS, "SimpleAdder.fmu"); type = :ME)
        notified = Ref(0)
        events = Ref(0)
        reset_pending = Ref(false)
        original_enter = fmu.cEnterEventMode
        original_update = fmu.cNewDiscreteStates
        original_get = fmu.cGetContinuousStates
        function completed_step(instance, no_reset, event, terminate)
            notified[] += 1
            unsafe_store!(event, Cint(cause == :step && notified[] == 1))
            unsafe_store!(terminate, Cint(cause == :terminate && notified[] == 1))
            return zero(Cuint)
        end
        function enter_event(instance)
            events[] += 1
            reset_pending[] = true
            return ccall(original_enter, Cuint, (Ptr{Cvoid},), instance)
        end
        function discrete_states(instance, info)
            status = ccall(original_update, Cuint, (Ptr{Cvoid}, Ptr{Cvoid}), instance, info)
            if reset_pending[]
                # valuesOfContinuousStatesChanged is the fourth fmi2Boolean in fmi2EventInfo.
                unsafe_store!(Ptr{Cint}(info), 1, 4)
            end
            return status
        end
        function get_states(instance, values, n)
            status = ccall(original_get, Cuint, (Ptr{Cvoid}, Ptr{Cdouble}, Csize_t), instance, values, n)
            if reset_pending[]
                unsafe_store!(values, 10.0)
                reset_pending[] = false
            end
            return status
        end
        complete_hook = @cfunction($completed_step, Cuint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}))
        enter_hook = @cfunction($enter_event, Cuint, (Ptr{Cvoid},))
        discrete_hook = @cfunction($discrete_states, Cuint, (Ptr{Cvoid}, Ptr{Cvoid}))
        get_hook = @cfunction($get_states, Cuint, (Ptr{Cvoid}, Ptr{Cdouble}, Csize_t))
        fmu.cCompletedIntegratorStep = Base.unsafe_convert(Ptr{Cvoid}, complete_hook)
        fmu.cEnterEventMode = Base.unsafe_convert(Ptr{Cvoid}, enter_hook)
        fmu.cNewDiscreteStates = Base.unsafe_convert(Ptr{Cvoid}, discrete_hook)
        fmu.cGetContinuousStates = Base.unsafe_convert(Ptr{Cvoid}, get_hook)
        if cause == :input
            only(filter(v -> v.name == "a", fmu.modelDescription.modelVariables)).variability = FMI.fmi2VariabilityDiscrete
        end
        GC.@preserve complete_hook enter_hook discrete_hook get_hook begin
            @named adder = MTK.FMIComponent(Val(2); fmu, type = :ME)
            input = cause == :input ? ifelse(t >= 0.015, 1.0, 0.0) : 1.0
            @mtkcompile sys = System([adder.a ~ input, adder.b ~ 2.0], t; systems = [adder])
            prob = ODEProblem(sys, [sys.adder.c => 0.0], (0.0, 0.04))
            sol = solve(prob, Tsit5(); dt = 0.01, adaptive = false)
            @test notified[] == sol.stats.naccept
            if cause == :terminate
                @test sol.retcode == SciMLBase.ReturnCode.Terminated
                @test sol.t[end] == 0.01
                @test events[] == 0
            else
                event_time = cause == :step ? 0.01 : 0.02
                @test SciMLBase.successful_retcode(sol)
                @test events[] == 1
                @test minimum(sol.t[findall(==(10.0), sol[sys.adder.c])]; init = Inf) == event_time
                @test sol[sys.adder.c][end] ≈ 10.0 + 4.0 * (0.04 - event_time)
            end
        end
    end
end


@testset "Fixed-step Co-Simulation endpoint grid" begin
    function instrumented_adder(finish)
        fmu = loadFMU(joinpath(FIXTURE_FMUS, "SimpleAdder.fmu"); type = :CS)
        fmu.modelDescription.coSimulation.canHandleVariableCommunicationStepSize = false
        requests = Tuple{Float64, Float64}[]
        original = fmu.cDoStep
        function do_step(instance, point, h, no_reset)
            push!(requests, (point, h))
            return ccall(original, Cuint, (Ptr{Cvoid}, Cdouble, Cdouble, Cint), instance, point, h, no_reset)
        end
        hook = @cfunction($do_step, Cuint, (Ptr{Cvoid}, Cdouble, Cdouble, Cint))
        fmu.cDoStep = Base.unsafe_convert(Ptr{Cvoid}, hook)
        return fmu, requests, hook
    end
    @testset "endpoint $finish" for finish in (0.02, nextfloat(0.02), 0.025)
        fmu, requests, hook = instrumented_adder(finish)
        GC.@preserve hook begin
            sys = adder_system(fmu, :CS; stop_time = 0.05)
            sol = solve(ODEProblem(sys, [sys.adder.c => 0.0], (0.0, finish)), Tsit5())
            @test SciMLBase.successful_retcode(sol)
            @test requests == [(0.0, 0.01), (0.01, 0.01)]
            @test sol[sys.adder.c][end] ≈ 0.08
        end
    end
    @testset "off-grid endpoint with experiment stop at the last communication" begin
        fmu, requests, hook = instrumented_adder(0.025)
        GC.@preserve hook begin
            sys = adder_system(fmu, :CS; stop_time = 0.02)
            sol = solve(ODEProblem(sys, [sys.adder.c => 0.0], (0.0, 0.025)), Tsit5())
            @test SciMLBase.successful_retcode(sol)
            @test sol.t[end] == 0.025
            @test requests == [(0.0, 0.01), (0.01, 0.01)]
            @test sol[sys.adder.c][end] ≈ 0.08
        end
    end
    @testset "large nonzero start" begin
        start = 1.0e8
        finish = start + 2 * 0.01
        fmu, requests, hook = instrumented_adder(finish)
        GC.@preserve hook begin
            sys = adder_system(fmu, :CS; stop_time = start + 0.05)
            sol = solve(ODEProblem(sys, [sys.adder.c => 0.0], (start, finish)), Tsit5())
            @test SciMLBase.successful_retcode(sol)
            @test length(requests) == 2
            @test first.(requests) == [start, start + 0.01]
            @test sum(last, requests) ≈ finish - start
            @test sol[sys.adder.c][end] ≈ 4 * (finish - start) atol = 8eps(start)
        end
    end
end

@testset "Coupled Co-Simulation FMUs publish before sampling" begin
    function coupled_solution(order)
        fmus = [loadFMU(joinpath(FIXTURE_FMUS, "SimpleAdder.fmu"); type = :CS) for _ in 1:2]
        @named adder1 = MTK.FMIComponent(Val(2); fmu = fmus[1], type = :CS, communication_step_size = 0.01)
        @named adder2 = MTK.FMIComponent(Val(2); fmu = fmus[2], type = :CS, communication_step_size = 0.01)
        systems = order == :forward ? [adder1, adder2] : [adder2, adder1]
        @mtkcompile sys = System(
            [adder1.a ~ adder2.c, adder2.a ~ adder1.c, adder1.b ~ 2.0, adder2.b ~ 2.0],
            t; systems
        )
        prob = ODEProblem(sys, [sys.adder1.c => 0.0, sys.adder2.c => 0.0], (0.0, 0.02))
        return sys, solve(prob, Tsit5())
    end
    for order in (:forward, :reverse)
        sys, sol = coupled_solution(order)
        @test SciMLBase.successful_retcode(sol)
        @test sol[sys.adder1.c][end] ≈ 0.0603
        @test sol[sys.adder2.c][end] ≈ 0.0603
        @test sol(0.01; idxs = sys.adder1.a) ≈ 0.03
        @test sol(0.01; idxs = sys.adder2.a) ≈ 0.03
    end
end

@testset "Model Exchange snapshots use the saved parameter values" begin
    fmu = loadFMU(joinpath(FIXTURE_FMUS, "SimpleAdder.fmu"); type = :ME)
    only(filter(v -> v.name == "a", fmu.modelDescription.modelVariables)).variability = FMI.fmi2VariabilityDiscrete
    completed = Ref(0)
    original_complete = fmu.cCompletedIntegratorStep
    function complete_step(instance, no_reset, event, terminate)
        completed[] += 1
        status = ccall(
            original_complete, Cuint,
            (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}), instance, no_reset, event, terminate
        )
        unsafe_store!(event, Cint(completed[] == 2))
        return status
    end
    complete_hook = @cfunction($complete_step, Cuint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}))
    fmu.cCompletedIntegratorStep = Base.unsafe_convert(Ptr{Cvoid}, complete_hook)
    GC.@preserve complete_hook begin
        @named adder = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @discretes host_a(t) = 1.0
        set_host_a(modified, observed, ctx, integrator) = (; host_a = 5.0)
        host_event = MTKBase.SymbolicDiscreteCallback(
            [0.02], MTKBase.ImperativeAffect(set_host_a; modified = (; host_a));
            reinitializealg = SciMLBase.NoInit()
        )
        @mtkcompile sys = System(
            [adder.a ~ host_a, adder.b ~ 0.0], t; systems = [adder], discrete_events = [host_event]
        )
        prob = ODEProblem(
            sys, [sys.adder.c => 0.0, sys.adder.value => 0.0], (0.0, 0.04)
        )
        sol = solve(prob, Tsit5(); dt = 0.01, adaptive = false, saveat = [0.015, 0.025])
        @test SciMLBase.successful_retcode(sol)
        event_indices = findall(==(0.02), sol.t)
        @test sol[sys.adder.out][first(event_indices)] == 1.0
        @test sol(0.02; idxs = sys.adder.out, continuity = :left) == 1.0
        @test sol(0.02; idxs = sys.adder.out, continuity = :right) == 5.0
        @test sol(0.015; idxs = sys.host_a) == 1.0
        @test sol(0.015; idxs = sys.adder.out) == 1.0
        @test sol(0.025; idxs = sys.host_a) == 5.0
        @test sol(0.025; idxs = sys.adder.out) == 5.0
    end
end

@testset "VanDerPol finite-difference Jacobian" begin
    fmu = loadFMU("VanDerPol", "ModelicaReferenceFMUs", "0.0.39", "3.0"; type = :ME)
    @mtkcompile sys = MTK.FMIComponent(Val(3); fmu, type = :ME)
    order = Int[findfirst(isequal(var), unknowns(sys)) for var in (sys.x0, sys.x1)]
    for start in (0.0, 1.0e8)
        prob = ODEProblem(sys, [sys.x0 => 1.0, sys.x1 => 2.0], (start, start + 1))
        rhs(u) = (du = similar(u); prob.f(du, u, prob.p, start); du)
        jac = FiniteDiff.finite_difference_jacobian(rhs, prob.u0)
        @test jac[order, order] ≈ [0.0 1.0; -5.0 0.0] atol = 1.0e-6
        @test_throws MethodError ForwardDiff.jacobian(rhs, prob.u0)
    end
end
