using ModelingToolkit, Test

# r is an input, and y is an output.
@variables t x(t)=0 y(t)=0 u(t)=0 r(t)=0
@variables t x(t)=0 y(t)=0 u(t)=0 r(t)=0 [input = true]
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

@test substitute(lsyss.A, ModelingToolkit.defaults(cl)) == lsys.A
@test substitute(lsyss.B, ModelingToolkit.defaults(cl)) == lsys.B
@test substitute(lsyss.C, ModelingToolkit.defaults(cl)) == lsys.C
@test substitute(lsyss.D, ModelingToolkit.defaults(cl)) == lsys.D
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

@test substitute(lsyss.A, ModelingToolkit.defaults(pid)) == lsys.A
@test substitute(lsyss.B, ModelingToolkit.defaults(pid)) == lsys.B
@test substitute(lsyss.C, ModelingToolkit.defaults(pid)) == lsys.C
@test substitute(lsyss.D, ModelingToolkit.defaults(pid)) == lsys.D

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

## Test that dummy_derivatives can be set to zero
if VERSION >= v"1.8"
    # The call to Link(; m = 0.2, l = 10, I = 1, g = -9.807) hangs forever on Julia v1.6
    using LinearAlgebra
    using ModelingToolkit
    using ModelingToolkitStandardLibrary
    using ModelingToolkitStandardLibrary.Blocks
    using ModelingToolkitStandardLibrary.Mechanical.MultiBody2D
    using ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition

    using ControlSystemsMTK
    using ControlSystemsMTK.ControlSystemsBase: sminreal, minreal, poles
    connect = ModelingToolkit.connect

    @parameters t
    D = Differential(t)

    @named link1 = Link(; m = 0.2, l = 10, I = 1, g = -9.807)
    @named cart = TranslationalPosition.Mass(; m = 1, s = 0)
    @named fixed = Fixed()
    @named force = Force(use_support = false)

    eqs = [connect(link1.TX1, cart.flange)
           connect(cart.flange, force.flange)
           connect(link1.TY1, fixed.flange)]

    @named model = ODESystem(eqs, t, [], []; systems = [link1, cart, force, fixed])
    def = ModelingToolkit.defaults(model)
    def[cart.s] = 10
    def[cart.v] = 0
    def[link1.A] = -pi / 2
    def[link1.dA] = 0
    lin_outputs = [cart.s, cart.v, link1.A, link1.dA]
    lin_inputs = [force.f.u]

    @info "named_ss"
    G = named_ss(model, lin_inputs, lin_outputs, allow_symbolic = true, op = def,
        allow_input_derivatives = true, zero_dummy_der = true)
    G = sminreal(G)
    @info "minreal"
    G = minreal(G)
    @info "poles"
    ps = poles(G)

    @test minimum(abs, ps) < 1e-6
    @test minimum(abs, complex(0, 1.3777260367206716) .- ps) < 1e-10

    lsys, syss = linearize(model, lin_inputs, lin_outputs, allow_symbolic = true, op = def,
        allow_input_derivatives = true, zero_dummy_der = true)
    lsyss, sysss = ModelingToolkit.linearize_symbolic(model, lin_inputs, lin_outputs;
        allow_input_derivatives = true)

    dummyder = setdiff(unknowns(sysss), unknowns(model))
    def = merge(def, Dict(x => 0.0 for x in dummyder))
    def[link1.fy1] = -def[link1.g] * def[link1.m]

    @test substitute(lsyss.A, def) ≈ lsys.A
    # We cannot pivot symbolically, so the part where a linear solve is required
    # is not reliable.
    @test substitute(lsyss.B, def)[1:6, 1] ≈ lsys.B[1:6, 1]
    @test substitute(lsyss.C, def) ≈ lsys.C
    @test substitute(lsyss.D, def) ≈ lsys.D
end
