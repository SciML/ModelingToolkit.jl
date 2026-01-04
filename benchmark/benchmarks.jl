using ModelingToolkit, BenchmarkTools
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEqDefault
using ModelingToolkit: t_nounits as t, D_nounits as D

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

# first call
mtkcompile(model)
SUITE["mtkcompile"] = @benchmarkable mtkcompile($model)

model = mtkcompile(model)
u0 = unknowns(model) .=> 0.0
tspan = (0.0, 6.0)

prob = ODEProblem(model, u0, tspan)
SUITE["ODEProblem"] = @benchmarkable ODEProblem($model, $u0, $tspan)

# first call
init(prob)
SUITE["init"] = @benchmarkable init($prob)

large_param_init = SUITE["large_parameter_init"] = BenchmarkGroup()

N = 25
@variables x(t)[1:N]
@parameters A[1:N, 1:N]

defval = collect(x) * collect(x)'
@mtkcompile model = System(
    [D(x) ~ x], t, [x], [A]; defaults = [A => defval], guesses = [A => fill(NaN, N, N)]
)

u0 = [x => rand(N)]
prob = ODEProblem(model, u0, tspan)
large_param_init["ODEProblem"] = @benchmarkable ODEProblem($model, $u0, $tspan)

large_param_init["init"] = @benchmarkable init($prob)

sparse_analytical_jacobian = SUITE["sparse_analytical_jacobian"]

eqs = [D(x[i]) ~ prod(x[j] for j in 1:N if (i + j) % 3 == 0) for i in 1:N]
@mtkcompile model = System(eqs, t)
u0 = collect(x .=> 1.0)
tspan = (0.0, 1.0)
jac = true
sparse = true
prob = ODEProblem(model, u0, tspan; jac, sparse)
out = similar(prob.f.jac_prototype)

sparse_analytical_jacobian["ODEProblem"] = @benchmarkable ODEProblem($model, $u0, $tspan; jac, sparse)
sparse_analytical_jacobian["f_oop"] = @benchmarkable $(prob.f.jac.f_oop)($(prob.u0), $(prob.p), $(first(tspan)))
sparse_analytical_jacobian["f_iip"] = @benchmarkable $(prob.f.jac.f_iip)($out, $(prob.u0), $(prob.p), $(first(tspan)))
