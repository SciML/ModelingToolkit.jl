using ModelingToolkit, BenchmarkTools
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Thermal
using OrdinaryDiffEqDefault

const SUITE = BenchmarkGroup()

@mtkmodel DCMotor begin
    @structural_parameters begin
        R = 0.5
        L = 4.5e-3
        k = 0.5
        J = 0.02
        f = 0.01
        V_step = 10
        tau_L_step = -3
    end
    @components begin
        ground = Ground()
        source = Voltage()
        voltage_step = Blocks.Step(height = V_step, start_time = 0)
        R1 = Resistor(R = R)
        L1 = Inductor(L = L, i = 0.0)
        emf = EMF(k = k)
        fixed = Fixed()
        load = Torque()
        load_step = Blocks.Step(height = tau_L_step, start_time = 3)
        inertia = Inertia(J = J)
        friction = Damper(d = f)
    end
    @equations begin
        connect(fixed.flange, emf.support, friction.flange_b)
        connect(emf.flange, friction.flange_a, inertia.flange_a)
        connect(inertia.flange_b, load.flange)
        connect(load_step.output, load.tau)
        connect(voltage_step.output, source.V)
        connect(source.p, R1.p)
        connect(R1.n, L1.p)
        connect(L1.n, emf.p)
        connect(emf.n, source.n, ground.g)
    end
end

@named model = DCMotor()

SUITE["mtkcompile"] = @benchmarkable mtkcompile($model)

model = mtkcompile(model)
u0 = unknowns(model) .=> 0.0
tspan = (0.0, 6.0)
SUITE["ODEProblem"] = @benchmarkable ODEProblem($model, $u0, $tspan)

prob = ODEProblem(model, u0, tspan)
SUITE["init"] = init($prob)
