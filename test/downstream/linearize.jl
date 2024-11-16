using ModelingToolkit, Test

# r is an input, and y is an output.
@independent_variables t
@variables x(t)=0 y(t)=0 u(t)=0 r(t)=0
@variables x(t)=0 y(t)=0 u(t)=0 r(t)=0 [input = true]
@parameters kp = 1
D = Differential(t)

eqs = [u ~ kp * (r - y)
       D(x) ~ -x + u
       y ~ x]

@named sys = ODESystem(eqs, t)

lsys, ssys = linearize(sys, [r], [y])

@test lsys.A[] == -2
@test lsys.B[] == 1
@test lsys.C[] == 1
@test lsys.D[] == 0

lsys, ssys = linearize(sys, [r], [r])

@test lsys.A[] == -2
@test lsys.B[] == 1
@test lsys.C[] == 0
@test lsys.D[] == 1

lsys, ssys = linearize(sys, r, r) # Test allow scalars

@test lsys.A[] == -2
@test lsys.B[] == 1
@test lsys.C[] == 0
@test lsys.D[] == 1

##
```

  r ┌─────┐       ┌─────┐     ┌─────┐
───►│     ├──────►│     │  u  │     │
    │  F  │       │  C  ├────►│  P  │ y
    └─────┘     ┌►│     │     │     ├─┬─►
                │ └─────┘     └─────┘ │
                │                     │
                └─────────────────────┘
```

function plant(; name)
    @variables x(t) = 1
    @variables u(t)=0 y(t)=0
    D = Differential(t)
    eqs = [D(x) ~ -x + u
           y ~ x]
    ODESystem(eqs, t; name = name)
end

function filt_(; name)
    @variables x(t)=0 y(t)=0
    @variables u(t)=0 [input = true]
    D = Differential(t)
    eqs = [D(x) ~ -2 * x + u
           y ~ x]
    ODESystem(eqs, t, name = name)
end

function controller(kp; name)
    @variables y(t)=0 r(t)=0 u(t)=0
    @parameters kp = kp
    eqs = [
        u ~ kp * (r - y)
    ]
    ODESystem(eqs, t; name = name)
end

@named f = filt_()
@named c = controller(1)
@named p = plant()

connections = [f.y ~ c.r # filtered reference to controller reference
               c.u ~ p.u # controller output to plant input
               p.y ~ c.y]

@named cl = ODESystem(connections, t, systems = [f, c, p])

lsys0, ssys = linearize(cl, [f.u], [p.x])
desired_order = [f.x, p.x]
lsys = ModelingToolkit.reorder_unknowns(lsys0, unknowns(ssys), desired_order)

@test lsys.A == [-2 0; 1 -2]
@test lsys.B == reshape([1, 0], 2, 1)
@test lsys.C == [0 1]
@test lsys.D[] == 0

## Symbolic linearization
lsyss, _ = ModelingToolkit.linearize_symbolic(cl, [f.u], [p.x])

@test ModelingToolkit.fixpoint_sub(lsyss.A, ModelingToolkit.defaults(cl)) == lsys.A
@test ModelingToolkit.fixpoint_sub(lsyss.B, ModelingToolkit.defaults(cl)) == lsys.B
@test ModelingToolkit.fixpoint_sub(lsyss.C, ModelingToolkit.defaults(cl)) == lsys.C
@test ModelingToolkit.fixpoint_sub(lsyss.D, ModelingToolkit.defaults(cl)) == lsys.D
##
using ModelingToolkitStandardLibrary.Blocks: LimPID
k = 400
Ti = 0.5
Td = 1
Nd = 10
@named pid = LimPID(; k, Ti, Td, Nd)

@unpack reference, measurement, ctr_output = pid
lsys0, ssys = linearize(pid, [reference.u, measurement.u], [ctr_output.u])
@unpack int, der = pid
desired_order = [int.x, der.x]
lsys = ModelingToolkit.reorder_unknowns(lsys0, unknowns(ssys), desired_order)

@test lsys.A == [0 0; 0 -10]
@test lsys.B == [2 -2; 10 -10]
@test lsys.C == [400 -4000]
@test lsys.D == [4400 -4400]

lsyss, _ = ModelingToolkit.linearize_symbolic(pid, [reference.u, measurement.u],
    [ctr_output.u])

@test ModelingToolkit.fixpoint_sub(
    lsyss.A, ModelingToolkit.defaults_and_guesses(pid)) == lsys.A
@test ModelingToolkit.fixpoint_sub(
    lsyss.B, ModelingToolkit.defaults_and_guesses(pid)) == lsys.B
@test ModelingToolkit.fixpoint_sub(
    lsyss.C, ModelingToolkit.defaults_and_guesses(pid)) == lsys.C
@test ModelingToolkit.fixpoint_sub(
    lsyss.D, ModelingToolkit.defaults_and_guesses(pid)) == lsys.D

# Test with the reverse desired unknown order as well to verify that similarity transform and reoreder_unknowns really works
lsys = ModelingToolkit.reorder_unknowns(lsys, unknowns(ssys), reverse(desired_order))

@test lsys.A == [-10 0; 0 0]
@test lsys.B == [10 -10; 2 -2]
@test lsys.C == [-4000 400]
@test lsys.D == [4400 -4400]

## Test that there is a warning when input is misspecified
if VERSION >= v"1.8"
    @test_throws "Some specified inputs were not found" linearize(pid,
        [
            pid.reference.u,
            pid.measurement.u
        ], [ctr_output.u])
    @test_throws "Some specified outputs were not found" linearize(pid,
        [
            reference.u,
            measurement.u
        ],
        [pid.ctr_output.u])
else # v1.6 does not have the feature to match error message
    @test_throws ErrorException linearize(pid,
        [
            pid.reference.u,
            pid.measurement.u
        ], [ctr_output.u])
    @test_throws ErrorException linearize(pid,
        [reference.u, measurement.u],
        [pid.ctr_output.u])
end

## Test operating points

# The saturation has no dynamics
function saturation(; y_max, y_min = y_max > 0 ? -y_max : -Inf, name)
    @variables u(t)=0 y(t)=0
    @parameters y_max=y_max y_min=y_min
    ie = ifelse
    eqs = [
    # The equation below is equivalent to y ~ clamp(u, y_min, y_max)
        y ~ ie(u > y_max, y_max, ie((y_min < u) & (u < y_max), u, y_min))
    ]
    ODESystem(eqs, t, name = name)
end
@named sat = saturation(; y_max = 1)
# inside the linear region, the function is identity
@unpack u, y = sat
lsys, ssys = linearize(sat, [u], [y])
@test isempty(lsys.A) # there are no differential variables in this system
@test isempty(lsys.B)
@test isempty(lsys.C)
@test lsys.D[] == 1

@test_skip lsyss, _ = ModelingToolkit.linearize_symbolic(sat, [u], [y]) # Code gen replaces ifelse with if statements causing symbolic evaluation to fail
# @test substitute(lsyss.A, ModelingToolkit.defaults(sat)) == lsys.A
# @test substitute(lsyss.B, ModelingToolkit.defaults(sat)) == lsys.B
# @test substitute(lsyss.C, ModelingToolkit.defaults(sat)) == lsys.C
# @test substitute(lsyss.D, ModelingToolkit.defaults(sat)) == lsys.D

# outside the linear region the derivative is 0
lsys, ssys = linearize(sat, [u], [y]; op = Dict(u => 2))
@test isempty(lsys.A) # there are no differential variables in this system
@test isempty(lsys.B)
@test isempty(lsys.C)
@test lsys.D[] == 0

# Test case when unknowns in system do not have equations in initialization system
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks: Add, Sine, PID, SecondOrder, Step, RealOutput
using ModelingToolkit: connect

# Parameters
m1 = 1
m2 = 1
k = 1000 # Spring stiffness
c = 10   # Damping coefficient
@named inertia1 = Inertia(; J = m1)
@named inertia2 = Inertia(; J = m2)
@named spring = Spring(; c = k)
@named damper = Damper(; d = c)
@named torque = Torque()

function SystemModel(u = nothing; name = :model)
    eqs = [connect(torque.flange, inertia1.flange_a)
           connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
           connect(inertia2.flange_a, spring.flange_b, damper.flange_b)]
    if u !== nothing
        push!(eqs, connect(torque.tau, u.output))
        return ODESystem(eqs, t;
            systems = [
                torque,
                inertia1,
                inertia2,
                spring,
                damper,
                u
            ],
            name)
    end
    ODESystem(eqs, t; systems = [torque, inertia1, inertia2, spring, damper], name)
end

@named r = Step(start_time = 0)
model = SystemModel()
@named pid = PID(k = 100, Ti = 0.5, Td = 1)
@named filt = SecondOrder(d = 0.9, w = 10)
@named sensor = AngleSensor()
@named er = Add(k2 = -1)

connections = [connect(r.output, :r, filt.input)
               connect(filt.output, er.input1)
               connect(pid.ctr_output, :u, model.torque.tau)
               connect(model.inertia2.flange_b, sensor.flange)
               connect(sensor.phi, :y, er.input2)
               connect(er.output, :e, pid.err_input)]

closed_loop = ODESystem(connections, t, systems = [model, pid, filt, sensor, r, er],
    name = :closed_loop, defaults = [
        model.inertia1.phi => 0.0,
        model.inertia2.phi => 0.0,
        model.inertia1.w => 0.0,
        model.inertia2.w => 0.0,
        filt.x => 0.0,
        filt.xd => 0.0
    ])

@test_nowarn linearize(closed_loop, :r, :y)

# https://discourse.julialang.org/t/mtk-change-in-linearize/115760/3
@mtkmodel Tank_noi begin
    # Model parameters
    @parameters begin
        ρ = 1, [description = "Liquid density"]
        A = 5, [description = "Cross sectional tank area"]
        K = 5, [description = "Effluent valve constant"]
        h_ς = 3, [description = "Scaling level in valve model"]
    end
    # Model variables, with initial values needed
    @variables begin
        m(t) = 1.5 * ρ * A, [description = "Liquid mass"]
        md_i(t), [description = "Influent mass flow rate"]
        md_e(t), [description = "Effluent mass flow rate"]
        V(t), [description = "Liquid volume"]
        h(t), [description = "level"]
    end
    # Providing model equations
    @equations begin
        D(m) ~ md_i - md_e
        m ~ ρ * V
        V ~ A * h
        md_e ~ K * sqrt(h / h_ς)
    end
end

@named tank_noi = Tank_noi()
@unpack md_i, h, m = tank_noi
m_ss = 2.4000000003229878
@test_nowarn linearize(tank_noi, [md_i], [h]; op = Dict(m => m_ss, md_i => 2))

# Test initialization
@variables x(t) y(t) u(t)=1.0
@parameters p = 1.0
eqs = [D(x) ~ p * u, x ~ y]
@named sys = ODESystem(eqs, t)

matrices1, _ = linearize(sys, [u], []; op = Dict(x => 2.0))
matrices2, _ = linearize(sys, [u], []; op = Dict(y => 2.0))
@test matrices1 == matrices2

# Ensure parameter values passed as `Dict` are respected
linfun, _ = linearization_function(sys, [u], []; op = Dict(x => 2.0))
matrices = linfun([1.0], Dict(p => 3.0), 1.0)
# this would be 1 if the parameter value isn't respected
@test matrices.f_u[] == 3.0

@testset "linearization_function handles dependent values" begin
    @parameters q
    matrices = @test_nowarn linfun([1.0], Dict(p => 3q, q => 1.0), 1.0)
    @test matrices.f_u[] == 3.0
end

@testset "Issue #2941" begin
    @variables x(t) y(t)
    @parameters p
    eqs = [0 ~ x * log(y) - p]
    @named sys = ODESystem(eqs, t; defaults = [p => 1.0])
    sys = complete(sys)
    @test_nowarn linearize(
        sys, [x], []; op = Dict(x => 1.0), allow_input_derivatives = true)
end
