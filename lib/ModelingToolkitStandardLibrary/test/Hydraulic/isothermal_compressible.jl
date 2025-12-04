using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC
import ModelingToolkitStandardLibrary.Blocks as B
import ModelingToolkitStandardLibrary.Mechanical.Translational as T

using ModelingToolkitStandardLibrary.Blocks: Parameter

NEWTON = NLNewton(
    check_div = false, always_new = true, max_iter = 100, relax = 9 // 10, κ = 1e-6)

@testset "Fluid Domain and Tube" begin
    function FluidSystem(N; bulk_modulus, name)
        pars = @parameters begin
            bulk_modulus = bulk_modulus
            p_int = 0
        end

        systems = @named begin
            fluid = IC.HydraulicFluid(; bulk_modulus)
            stp = B.Step(; height = 10e5, offset = 0, start_time = 0.005,
                duration = Inf, smooth = true)
            src = IC.Pressure()
            vol = IC.FixedVolume(; vol = 10.0, p_int)
            res = IC.Tube(N; area = 0.01, length = 50.0, p_int)
        end

        eqs = [connect(stp.output, src.p)
               connect(fluid, src.port)
               connect(src.port, res.port_a)
               connect(res.port_b, vol.port)]

        System(eqs, t, [], pars; name, systems)
    end

    @mtkcompile s1_1 = FluidSystem(1; bulk_modulus = 1e9)
    @mtkcompile s1_2 = FluidSystem(1; bulk_modulus = 2e9)
    @mtkcompile s5_1 = FluidSystem(5; bulk_modulus = 1e9)

    p1_1 = ODEProblem(s1_1, [], (0, 0.05))
    p1_2 = ODEProblem(s1_2, [], (0, 0.05))
    p5_1 = ODEProblem(s5_1, [], (0, 0.05))

    sol1_1 = solve(p1_1, Rodas5P())
    sol1_2 = solve(p1_2, Rodas5P())
    sol5_1 = solve(p5_1, Rodas5P())

    # fig = Figure()
    # tm = 0:0.001:0.05 |> collect
    # ax = Axis(fig[1,1])
    # lines!(ax, tm, sol1_1.(tm; idxs=s1_2.vol.port.p)); fig
    # lines!(ax, tm, sol1_2.(tm; idxs=s1_1.vol.port.p)); fig
    # lines!(ax, tm, sol5_1.(tm; idxs=s5_1.vol.port.p)); fig
    # fig

    # higher stiffness should compress more quickly and give a higher pressure
    @test sol1_2[s1_2.vol.port.p][end] > sol1_1[s1_1.vol.port.p][end]

    # N=5 pipe is compressible, will pressurize more slowly
    @test sol1_1[s1_1.vol.port.p][end] > sol5_1[s5_1.vol.port.p][end]
end

@testset "Valve" begin
    function ValveSystem(; name)
        pars = []

        systems = @named begin
            fluid = IC.HydraulicFluid()
            sink = IC.FixedPressure(; p = 10e5)
            vol = IC.FixedVolume(; vol = 5, p_int = 1e5)
            valve = IC.Valve(; Cd = 1e5, minimum_area = 0)
            ramp = B.Ramp(;
                height = 0.1, duration = 0.1, offset = 0, start_time = 0.1, smooth = true)
        end

        eqs = [connect(fluid, sink.port)
               connect(sink.port, valve.port_a)
               connect(valve.port_b, vol.port)
               connect(valve.area, ramp.output)]

        System(eqs, t, [], pars; name, systems)
    end

    @named valve_system = ValveSystem()
    sys = mtkcompile(valve_system)
    prob = ODEProblem(sys, [], (0, 1))
    sol = solve(prob, Rodas5P(); abstol = 1e-6, reltol = 1e-9)
    s = complete(valve_system)

    # the volume should discharge to 10bar
    @test sol[s.vol.port.p][end]≈10e5 atol=1e5

    # fig = Figure()
    # tm = 0:0.01:1 |> collect
    # ax = Axis(fig[1,1])
    # lines!(ax, tm, sol.(tm; idxs=sys.vol.port.p));  
    # fig
end

@testset "DynamicVolume and minimum_volume feature" begin # Need help here
    function TestSystem(; name, area = 0.01, length = 0.1, damping_volume = length * area *
                                                                            0.1)
        pars = []

        # DynamicVolume values
        systems = @named begin
            fluid = IC.HydraulicFluid(; bulk_modulus = 1e9)

            src1 = IC.Pressure(;)
            src2 = IC.Pressure(;)

            vol1 = IC.DynamicVolume(; direction = +1,
                area,
                x_int = length,
                x_max = length * 2,
                x_min = length * 0.1,
                x_damp = damping_volume / area + length * 0.1,
                d = 1e3,
                p_int = 10e5)
            # vol1 = IC.Volume(;area, direction = +1, x_int=length)

            vol2 = IC.DynamicVolume(; direction = -1,
                area,
                x_int = length,
                x_max = length * 2,
                x_min = length * 0.1,
                x_damp = damping_volume / area + length * 0.1,
                d = 1e3,
                p_int = 10e5)
            # vol2 = IC.Volume(;area, direction = -1, x_int=length)

            mass = T.Mass(; m = 10)

            sin1 = B.Sine(; frequency = 0.5, amplitude = +1e5, offset = 10e5)
            sin2 = B.Sine(; frequency = 0.5, amplitude = -1e5, offset = 10e5)
        end

        eqs = [connect(fluid, src1.port)
               connect(fluid, src2.port)
               connect(src1.port, vol1.port)
               connect(src2.port, vol2.port)
               connect(vol1.flange, mass.flange, vol2.flange)
               connect(src1.p, sin1.output)
               connect(src2.p, sin2.output)]

        initialization_eqs = [mass.s ~ 0.0
                              mass.v ~ 0.0]

        System(eqs, t, [], pars; name, systems, initialization_eqs)
    end

    @named sys = TestSystem()
    sys = mtkcompile(sys; allow_symbolic = true)
    prob = ODEProblem(sys, [], (0, 5))
    sol = solve(prob, Rodas5P(); abstol = 1e-6, reltol = 1e-9)
    # begin
    #     fig = Figure()

    #     ax = Axis(fig[1,1], ylabel="position [m]", xlabel="time [s]")
    #     lines!(ax, sol.t, sol[sys.vol1.x]; label="vol1")
    #     lines!(ax, sol.t, sol[sys.vol2.x]; label="vol2")
    #     Legend(fig[1,2], ax)

    #     ax = Axis(fig[2,1], ylabel="pressure [bar]", xlabel="time [s]")
    #     lines!(ax, sol.t, sol[sys.vol1.damper.port_a.p]/1e5; label="vol1")
    #     lines!(ax, sol.t, sol[sys.vol2.damper.port_a.p]/1e5; label="vol2")
    #     ylims!(ax, 10-2, 10+2)

    #     ax = Axis(fig[3,1], ylabel="area", xlabel="time [s]")
    #     lines!(ax, sol.t, sol[sys.vol1.damper.area]; label="area 1")
    #     lines!(ax, sol.t, sol[sys.vol2.damper.area]; label="area 2")

    #     display(fig)
    # end

    # volume/mass should stop moving at opposite ends
    @test sol(0; idxs = sys.vol1.x) == 0.1
    @test sol(0; idxs = sys.vol2.x) == 0.1

    @test round(sol(1; idxs = sys.vol1.x); digits = 2) == 0.19
    @test round(sol(1; idxs = sys.vol2.x); digits = 2) == 0.01

    @test round(sol(2; idxs = sys.vol1.x); digits = 2) == 0.01
    @test round(sol(2; idxs = sys.vol2.x); digits = 2) == 0.19

    @test round(sol(3; idxs = sys.vol1.x); digits = 2) == 0.19
    @test round(sol(3; idxs = sys.vol2.x); digits = 2) == 0.01

    @test round(sol(4; idxs = sys.vol1.x); digits = 2) == 0.01
    @test round(sol(4; idxs = sys.vol2.x); digits = 2) == 0.19
end

@testset "Actuator System" begin
    function ActuatorSystem(use_input; name)
        pars = @parameters begin
            p_s = 200e5
            p_r = 5e5

            A_1 = 360e-4
            A_2 = 360e-4

            p_1 = 45e5
            p_2 = 45e5

            l_1 = 1.5
            l_2 = 1.5
            m_f = 250
            g = 0

            d = 100e-3

            Cd = 0.01

            m_piston = 880
        end

        vars = @variables begin
            ddx(t) = 0
        end

        systems = @named begin
            src = IC.FixedPressure(; p = p_s)
            valve = IC.SpoolValve2Way(; g, m = m_f, d, Cd, x_int = 0)
            piston = IC.Actuator(;
                length_a_int = l_1,
                length_b_int = l_2,
                area_a = A_1,
                area_b = A_2,
                m = m_piston,
                g = 0,
                minimum_volume_a = A_1 * 1e-3,
                minimum_volume_b = A_2 * 1e-3,
                damping_volume_a = A_1 * 5e-3,
                damping_volume_b = A_2 * 5e-3,
                p_a_int = p_1,
                p_b_int = p_2)
            # body = T.Mass(; m = 1500)
            # pipe = IC.Tube(1; area = A_2, length = 2.0, p_int = p_2)
            snk = IC.FixedPressure(; p = p_r)
            pos = T.Position()

            # m1 = IC.FlowDivider(; n = 3)
            # m2 = IC.FlowDivider(; n = 3)

            fluid = IC.HydraulicFluid()
        end

        if use_input
            @named input = B.SampledData(Float64)
        else
            #@named input = B.TimeVaryingFunction(f)
            @named input = B.Constant(k = 0)
        end

        push!(systems, input)

        eqs = [connect(input.output, pos.s)
               connect(valve.flange, pos.flange)
               connect(valve.port_a, piston.port_a)
               #    connect(piston.flange, body.flange)

               connect(piston.port_b, valve.port_b)

               #    connect(piston.port_b, pipe.port_b)
               # #    connect(piston.port_b, m1.port_a)
               # #    connect(m1.port_b, pipe.port_b)

               # connect(pipe.port_a, valve.port_b)
               # #    connect(pipe.port_a, m2.port_b)
               # #    connect(m2.port_a, valve.port_b)

               connect(src.port, valve.port_s)
               connect(snk.port, valve.port_r)
               connect(fluid, src.port, snk.port)
               D(piston.mass.v) ~ ddx]

        initialization_eqs = [
        # body.s ~ 0
        ]

        System(eqs, t, vars, pars; name, systems, initialization_eqs)
    end

    @mtkcompile initsys = ActuatorSystem(false)

    initprob = ODEProblem(initsys, [], (0, 0))
    initsol = solve(initprob, Rodas5P())

    @mtkcompile sys = ActuatorSystem(true)

    dt = 1e-4
    time = 0:dt:0.1

    x = @. (time - 0.015)^2 - 10 * (time - 0.02)^3
    x[1:150] = zeros(150)

    defs = ModelingToolkit.defaults(sys)
    defs[sys.input.buffer] = Parameter(0.5 * x, dt)

    # NOTE: bypassing initialization system: https://github.com/SciML/ModelingToolkit.jl/issues/3312
    prob = ODEProblem(sys, initsol[1], (0, 0.1); build_initializeprob = false)

    #TODO: Implement proper initialization system after issue is resolved
    #TODO: How to bring the body back and not have an overdetermined system?

    # check the fluid domain
    @test Symbol(defs[sys.src.port.ρ]) == Symbol(sys.fluid.ρ)
    @test Symbol(defs[sys.valve.port_s.ρ]) == Symbol(sys.fluid.ρ)
    @test Symbol(defs[sys.valve.port_a.ρ]) == Symbol(sys.fluid.ρ)
    @test Symbol(defs[sys.valve.port_b.ρ]) == Symbol(sys.fluid.ρ)
    @test Symbol(defs[sys.valve.port_r.ρ]) == Symbol(sys.fluid.ρ)
    @test Symbol(defs[sys.snk.port.ρ]) == Symbol(sys.fluid.ρ)

    @time sol = solve(prob, Rodas5P(); initializealg = NoInit())

    @test sol[sys.ddx][1] == 0.0
    @test maximum(sol[sys.ddx]) > 200
    @test sol[sys.piston.x][end] > 0.6
end

@testset "Prevent Negative Pressure" begin
    @component function HydraulicSystem(; name)
        pars = @parameters let_gas = 1

        systems = @named begin
            fluid = IC.HydraulicFluid(; let_gas)
            vol = IC.DynamicVolume(; area = 0.001, x_int = 0.05,
                x_max = 0.1, x_damp = 0.02, x_min = 0.01, direction = +1, p_int = 100e5)
            mass = T.Mass(; m = 100, g = -9.807) # s = 0.05
            cap = IC.Cap()
        end

        eqs = [connect(fluid, cap.port, vol.port)
               connect(vol.flange, mass.flange)]

        initialization_eqs = [mass.s ~ 0.05
                              mass.v ~ 0]

        return System(eqs, t, [], pars; name, systems, initialization_eqs)
    end

    @mtkcompile sys = HydraulicSystem()

    prob1 = ODEProblem(sys, [], (0, 0.05))
    # prob1 = remake(prob1; u0 = BigFloat.(prob1.u0))
    prob2 = ODEProblem(sys, [sys.let_gas => 0], (0, 0.05))

    # @time sol1 = solve(prob1, Rodas5P(); abstol=1e-9, reltol=1e-9) #BUG: Using BigFloat gives... ERROR: MethodError: no method matching getindex(::Missing, ::Int64)
    @time sol1 = solve(prob1, Rodas5P(); adaptive = false, dt = 1e-6) #TODO: fix BigFloat to implement abstol=1e-9, reltol=1e-9
    @time sol2 = solve(prob2, Rodas5P())

    # case 1: no negative pressure will only have gravity pulling mass back down
    # case 2: with negative pressure, added force pulling mass back down
    # - case 1 should push the mass higher
    @test sol1[sys.mass.s][end] > sol2[sys.mass.s][end]

    # case 1 should prevent negative pressure less than -1000
    @test minimum(sol1[sys.vol.port.p]) > -5000
    @test minimum(sol2[sys.vol.port.p]) < -5000

    # fig = Figure()
    # ax = Axis(fig[1,1])
    # lines!(ax, sol1.t, sol1[sys.vol.port.p]); fig
    # lines!(ax, sol2.t, sol2[sys.vol.port.p]); fig

    # ax = Axis(fig[1,2])
    # lines!(ax, sol1.t, sol1[sys.mass.s])
    # lines!(ax, sol2.t, sol2[sys.mass.s])
    # fig
end

#TODO
# @testset "Component Flow Reversals" begin
# # Check Component Flow Reversals
#     function System(; name)
#         pars = []

#         systems = @named begin
#             fluid = IC.HydraulicFluid()
#             source = IC.Pressure()
#             sink = IC.FixedPressure(; p = 101325)
#             pipe = IC.Tube(1, false; area = 0.1, length =.1, head_factor = 1)
#             osc = Sine(; frequency = 0.01, amplitude = 100, offset = 101325)
#         end

#         eqs = [connect(fluid, pipe.port_a)
#             connect(source.port, pipe.port_a)
#             connect(pipe.port_b, sink.port)
#             connect(osc.output, source.p)]

#         System(eqs, t, [], []; systems)
#     end

#     @named sys = System()

#     syss = mtkcompile.([sys])
#     tspan = (0.0, 1000.0)
#     prob = ODEProblem(sys, tspan)  # u0 guess can be supplied or not
#     @time sol = solve(prob)

# end

#TODO
# @testset "Tube Discretization" begin
#     # Check Tube Discretization
# end

#TODO
# @testset "Pressure BC" begin
#     # Ensure Pressure Boundary Condition Works
# end

#TODO
# @testset "Massflow BC" begin
#     # Ensure Massflow Boundary Condition Works
# end

#TODO
# @testset "Splitter Flow Test" begin
#     # Ensure FlowDivider Splits Flow Properly
#     # 1) Set flow into port A, expect reduction in port B

#     # 2) Set flow into port B, expect increase in port B
# end

#TODO: Test Valve Inversion
