using Test, ModelingToolkit, FMI, FMIZoo, OrdinaryDiffEq, SciMLBase
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit as MTK
import FiniteDiff, ForwardDiff

const FIXTURE_FMUS = joinpath(@__DIR__, "fmus")

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

@testset "FMI 3 unsupported wide integer ports" begin
    @testset "$type" for type in (:ME, :CS)
        fmu = loadFMU("Feedthrough", "ModelicaReferenceFMUs", "0.0.39", "3.0"; type)
        step = type == :CS ? (; communication_step_size = 0.01) : (;)
        ports = @test_logs match_mode = :any (:warn, r"FMI 3 Int64 variables Int64_input") (
            :warn, r"FMI 3 UInt64 variables UInt64_input",
        ) MTK.FMIComponent(Val(3); fmu, type, name = :ports, step...)
        names = Set(MTK.getname.(unknowns(ports)))
        @test isempty(intersect(names, (:Int64_input, :UInt64_input, :Int64_output, :UInt64_output)))
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
        original = fmu.cDoStep
        requests = Tuple{Float64, Float64}[]
        function record_step(instance, point, h, no_reset, event, terminate, early, last_time)
            push!(requests, (point, h))
            return ccall(
                original, Cuint,
                (Ptr{Cvoid}, Cdouble, Cdouble, UInt8, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{Cdouble}),
                instance, point, h, no_reset, event, terminate, early, last_time
            )
        end
        hook = @cfunction(
            $record_step, Cuint,
            (Ptr{Cvoid}, Cdouble, Cdouble, UInt8, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{Cdouble})
        )
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
            @test sol[sys.ports.Float64_continuous_output][end] ≈ 4 * (finish - start)
            @test sum(last, requests; init = 0.0) ≈ finish - start
            if !variable
                @test requests == [(0.0, 0.01), (0.01, 0.01)]
                @test_throws ArgumentError solve(remake(prob; tspan = (0.0, 0.025)), Tsit5())
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
