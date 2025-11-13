using ModelingToolkitStandardLibrary.Magnetic, ModelingToolkit, OrdinaryDiffEq, Test

import ModelingToolkitStandardLibrary.Electrical
import ModelingToolkitStandardLibrary.Blocks
import ModelingToolkitStandardLibrary.Magnetic
using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq: ReturnCode.Success
# using Plots

@testset "Inductor" begin
    mu_air = 1
    l_air = 0.0001
    mu_Fe = 1000
    l_Fe = 4 * 0.065
    a = b = 0.25

    @named source = Blocks.Sine(amplitude = 230 * sqrt(2), frequency = 50, phase = pi / 2)
    @named voltage = Electrical.Voltage()
    @named r = Electrical.Resistor(R = 7.5)
    @named ground = Electrical.Ground()
    @named coil = Magnetic.FluxTubes.ElectroMagneticConverter(N = 600, Phi = 0.0)
    @named ground_m = Magnetic.FluxTubes.Ground()
    @named r_mAirPar = Magnetic.FluxTubes.ConstantReluctance(R_m = a * b * l_air * mu_air)
    @named r_mFe = Magnetic.FluxTubes.ConstantReluctance(R_m = a * b * l_Fe * mu_Fe)
    @named r_mLeak = Magnetic.FluxTubes.ConstantReluctance(R_m = 1.2e6)
    connections = [connect(source.output, voltage.V)
                   connect(voltage.p, r.p)
                   connect(r.n, coil.p)
                   connect(voltage.n, coil.n)
                   connect(coil.port_p, r_mLeak.port_p)
                   connect(r_mLeak.port_p, r_mAirPar.port_p)
                   connect(r_mAirPar.port_n, r_mFe.port_p)
                   connect(r_mFe.port_n, r_mLeak.port_n)
                   connect(r_mFe.port_n, coil.port_n)
                   connect(ground.g, voltage.n)
                   connect(ground_m.port, r_mFe.port_n)]
    @named model = System(connections, t,
        systems = [
            source,
            r,
            ground,
            coil,
            ground_m,
            r_mAirPar,
            r_mFe,
            r_mLeak,
            voltage
        ])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[], (0, 0.1))
    sol = solve(prob, Rodas4())

    # Plots.plot(sol; vars=[r.i])
    # Plots.plot(sol; vars=[r_mFe.V_m, r_mFe.Phi])

    @test SciMLBase.successful_retcode(sol)
    @test sol[r_mFe.Phi] == sol[r_mAirPar.Phi]
    @test all(sol[coil.port_p.Phi] + sol[r_mLeak.Phi] + sol[r_mAirPar.Phi] .== 0)
end
