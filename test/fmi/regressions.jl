using Test, ModelingToolkit, FMI, FMIZoo, OrdinaryDiffEq
using ModelingToolkit: SciMLBase, t_nounits as t, D_nounits as D
import ModelingToolkit as MTK

const FIXTURE_FMUS = joinpath(@__DIR__, "fmus")

function load_feedthrough()
    fmu = loadFMU("Feedthrough", "ModelicaReferenceFMUs", "0.0.39", "3.0"; type = :CS)
    md = fmu.modelDescription
    # Use the Float64 ports to exercise communication independently of typed ports.
    for refs in (md.inputValueReferences, md.outputValueReferences, md.parameterValueReferences)
        filter!(ref -> FMI.dataTypeForValueReference(fmu, ref) === Float64, refs)
    end
    return fmu
end

function feedthrough_step_hook(step; stop_time = nothing)
    # Feedthrough has no state to advance. Supply a communication clock for arbitrary intervals.
    function do_step(instance, point, h, no_reset, event, terminate, early, last_time)
        unsafe_store!(event, UInt8(0))
        unsafe_store!(terminate, UInt8(0))
        unsafe_store!(early, UInt8(0))
        achieved = step(point, h)
        if stop_time !== nothing && point + h > stop_time[]
            return Cuint(FMI.fmi3StatusError)
        end
        unsafe_store!(last_time, achieved)
        return Cuint(FMI.fmi3StatusOK)
    end
    return @cfunction(
        $do_step, Cuint,
        (Ptr{Cvoid}, Cdouble, Cdouble, UInt8, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{Cdouble})
    )
end

function feedthrough_init_hook(fmu, limit)
    original_enter = fmu.cEnterInitializationMode
    function enter_init(instance, tolerance_defined, tolerance, t0, stop_defined, stop)
        limit[] = iszero(stop_defined) ? Inf : stop
        return ccall(
            original_enter, Cuint, (Ptr{Cvoid}, UInt8, Cdouble, Cdouble, UInt8, Cdouble),
            instance, tolerance_defined, tolerance, t0, stop_defined, stop
        )
    end
    return @cfunction($enter_init, Cuint, (Ptr{Cvoid}, UInt8, Cdouble, Cdouble, UInt8, Cdouble))
end

function adder_system(fmu; stop_time)
    @named adder = MTK.FMIComponent(Val(2); fmu, type = :CS, communication_step_size = 0.01, stop_time)
    @variables z(t) = 0.0
    @mtkcompile sys = System([adder.a ~ 1.0, adder.b ~ 2.0, D(z) ~ 1.0], t; systems = [adder])
    return sys
end

function feedthrough_problem(fmu, tspan; communication_step_size = 0.01, stop_time = nothing)
    @named ports = MTK.FMIComponent(Val(3); fmu, type = :CS, communication_step_size, stop_time)
    inputs = [
        getproperty(ports, MTK.getname(v)) ~ 0.0
            for v in unknowns(ports) if MTK.isinput(v)
    ]
    @variables z(t) = 0.0
    @mtkcompile sys = System([inputs; D(z) ~ 1.0], t; systems = [ports])
    return ODEProblem(sys, [], tspan)
end

@testset "Co-Simulation scheduling and limits" begin
    @testset "Co-Simulation communication reaches the solve endpoint" begin
        @testset "FMI 2, $start to $finish" for (start, finish) in ((0.0, 0.025), (0.0, 1.0e-9), (1.0e8, 1.0e8 + 0.025))
            fmu = loadFMU(joinpath(FIXTURE_FMUS, "SimpleAdder.fmu"); type = :CS)
            sys = adder_system(fmu; stop_time = start + 0.05)
            sol = solve(ODEProblem(sys, [sys.adder.c => 0.0], (start, finish)), Tsit5())
            @test SciMLBase.successful_retcode(sol)
            @test sol[sys.adder.c][end] ≈ 4 * (finish - start)
        end
        @testset "FMI 3, $start to $finish, variable steps: $variable" for (start, finish, variable) in (
                (0.0, 0.025, true), (0.0, 1.0e-9, true), (1.0e8, 1.0e8 + 0.025, true), (0.0, 0.02, false),
            )
            fmu = load_feedthrough()
            fmu.modelDescription.coSimulation.canHandleVariableCommunicationStepSize = variable
            requests = Tuple{Float64, Float64}[]
            hook = feedthrough_step_hook() do point, h
                push!(requests, (point, h))
                return point + h
            end
            fmu.cDoStep = Base.unsafe_convert(Ptr{Cvoid}, hook)
            GC.@preserve hook begin
                prob = feedthrough_problem(fmu, (start, finish); stop_time = start + 0.05)
                sol = solve(prob, Tsit5())
                @test SciMLBase.successful_retcode(sol)
                @test sol.t[end] == finish
                @test sum(last, requests; init = 0.0) ≈ finish - start
                if !variable
                    @test requests == [(0.0, 0.01), (0.01, 0.01)]
                    offgrid = solve(remake(prob; tspan = (0.0, 0.025)), Tsit5())
                    @test SciMLBase.successful_retcode(offgrid)
                end
            end
        end
    end

    @testset "Variable communication rounding at the experiment limit" begin
        start, finish = 0.001, 0.01
        @test start + (finish - start) > finish
        @testset "stop_time = $stop_time" for stop_time in (nothing, finish)
            fmu = load_feedthrough()
            fmu.modelDescription.coSimulation.canHandleVariableCommunicationStepSize = true
            limit, short = Ref(Inf), Ref(false)
            requests = Tuple{Float64, Float64}[]
            init_hook = feedthrough_init_hook(fmu, limit)
            step_hook = feedthrough_step_hook(; stop_time = limit) do point, h
                push!(requests, (point, h))
                return short[] ? point + h / 2 : point + h
            end
            fmu.cEnterInitializationMode = Base.unsafe_convert(Ptr{Cvoid}, init_hook)
            fmu.cDoStep = Base.unsafe_convert(Ptr{Cvoid}, step_hook)
            GC.@preserve init_hook step_hook begin
                prob = feedthrough_problem(fmu, (start, finish); stop_time)
                sol = solve(prob, Tsit5())
                @test SciMLBase.successful_retcode(sol)
                @test sol.t[end] == finish
                @test limit[] == finish
                @test length(requests) == 1
                point, h = only(requests)
                @test point == start && h > 0
                @test prevfloat(finish) <= point + h <= limit[]
                short[] = true
                @test_throws r"FMI 3.*lastSuccessfulTime" solve(prob, Tsit5())
            end
        end
    end

    @testset "Fixed-step Co-Simulation endpoint grid" begin
        function instrumented_adder(; strict = false)
            fmu = loadFMU(joinpath(FIXTURE_FMUS, "SimpleAdder.fmu"); type = :CS)
            fmu.modelDescription.coSimulation.canHandleVariableCommunicationStepSize = false
            requests = Tuple{Float64, Float64}[]
            original = fmu.cDoStep
            function do_step(instance, point, h, no_reset)
                push!(requests, (point, h))
                strict && h != 0.01 && return Cuint(FMI.fmi2StatusError)
                return ccall(original, Cuint, (Ptr{Cvoid}, Cdouble, Cdouble, Cint), instance, point, h, no_reset)
            end
            hook = @cfunction($do_step, Cuint, (Ptr{Cvoid}, Cdouble, Cdouble, Cint))
            fmu.cDoStep = Base.unsafe_convert(Ptr{Cvoid}, hook)
            return fmu, requests, hook
        end
        @testset "endpoint $finish" for finish in (0.02, nextfloat(0.02), 0.025)
            fmu, requests, hook = instrumented_adder()
            GC.@preserve hook begin
                sys = adder_system(fmu; stop_time = 0.05)
                sol = solve(ODEProblem(sys, [sys.adder.c => 0.0], (0.0, finish)), Tsit5())
                @test SciMLBase.successful_retcode(sol)
                @test requests == [(0.0, 0.01), (0.01, 0.01)]
                @test sol[sys.adder.c][end] ≈ 0.08
            end
        end
        @testset "off-grid endpoint with experiment stop at the last communication" begin
            fmu, requests, hook = instrumented_adder()
            GC.@preserve hook begin
                sys = adder_system(fmu; stop_time = 0.02)
                sol = solve(ODEProblem(sys, [sys.adder.c => 0.0], (0.0, 0.025)), Tsit5())
                @test SciMLBase.successful_retcode(sol)
                @test sol.t[end] == 0.025
                @test requests == [(0.0, 0.01), (0.01, 0.01)]
                @test sol[sys.adder.c][end] ≈ 0.08
            end
        end
        @testset "strict fixed-step native requests" begin
            fmu, requests, hook = instrumented_adder(; strict = true)
            GC.@preserve hook begin
                sys = adder_system(fmu; stop_time = 0.05)
                sol = solve(ODEProblem(sys, [sys.adder.c => 0.0], (0.0, 0.035)), Tsit5())
                @test SciMLBase.successful_retcode(sol)
                @test requests == [(0.0, 0.01), (0.01, 0.01), (0.02, 0.01)]
                @test sol[sys.adder.c][end] ≈ 0.12
            end
            start = 1.0e8
            finish = start + 2 * 0.01
            fmu, requests, hook = instrumented_adder(; strict = true)
            GC.@preserve hook begin
                sys = adder_system(fmu; stop_time = start + 0.05)
                sol = solve(ODEProblem(sys, [sys.adder.c => 0.0], (start, finish)), Tsit5())
                @test SciMLBase.successful_retcode(sol)
                @test requests == [(start, 0.01), (start + 0.01, 0.01)]
                @test sol[sys.adder.c][end] ≈ 4 * (finish - start) atol = 8eps(start)
            end
            fmu3 = load_feedthrough()
            fmu3.modelDescription.coSimulation.canHandleVariableCommunicationStepSize = false
            requests3 = Tuple{Float64, Float64}[]
            hook3 = feedthrough_step_hook() do point, h
                push!(requests3, (point, h))
                return point + h
            end
            fmu3.cDoStep = Base.unsafe_convert(Ptr{Cvoid}, hook3)
            GC.@preserve hook3 begin
                prob = feedthrough_problem(fmu3, (start, finish); stop_time = start + 0.05)
                sol = solve(prob, Tsit5())
                @test SciMLBase.successful_retcode(sol)
                @test requests3 == [(start, 0.01), (start + 0.01, 0.01)]
            end
            coarse_h = 3.5 * eps(start)
            coarse_finish = start + 2 * coarse_h
            fmu3_coarse = load_feedthrough()
            fmu3_coarse.modelDescription.coSimulation.canHandleVariableCommunicationStepSize = false
            coarse_requests = Tuple{Float64, Float64}[]
            coarse_short = Ref(false)
            coarse_limit = Ref(Inf)
            init_hook = feedthrough_init_hook(fmu3_coarse, coarse_limit)
            coarse_hook = feedthrough_step_hook(; stop_time = coarse_limit) do point, step_h
                push!(coarse_requests, (point, step_h))
                return coarse_short[] ? point + 0.5 * step_h : point + step_h
            end
            fmu3_coarse.cDoStep = Base.unsafe_convert(Ptr{Cvoid}, coarse_hook)
            fmu3_coarse.cEnterInitializationMode = Base.unsafe_convert(Ptr{Cvoid}, init_hook)
            native_finish = (start + coarse_h) + coarse_h
            GC.@preserve coarse_hook init_hook begin
                @testset "coarse fixed steps (stop_time = $stop_time)" for stop_time in (nothing, coarse_finish)
                    empty!(coarse_requests)
                    coarse_prob = feedthrough_problem(
                        fmu3_coarse, (start, coarse_finish); communication_step_size = coarse_h, stop_time
                    )
                    if stop_time === nothing
                        coarse_sol = solve(coarse_prob, Tsit5())
                        @test SciMLBase.successful_retcode(coarse_sol)
                        @test coarse_sol.t[end] == coarse_finish
                        @test coarse_limit[] == native_finish > coarse_finish
                        @test coarse_requests == [(start, coarse_h), (start + coarse_h, coarse_h)]
                        coarse_short[] = true
                        @test_throws r"FMI 3.*lastSuccessfulTime" solve(coarse_prob, Tsit5())
                    else
                        coarse_short[] = false
                        @test_throws "stop_time must be at least $native_finish" solve(coarse_prob, Tsit5())
                        @test isempty(coarse_requests)
                    end
                end
            end
        end
    end
end
