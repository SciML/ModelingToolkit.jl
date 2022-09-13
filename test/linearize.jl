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
        u ~ kp * (r - y),
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

lsys, ssys = linearize(cl, [f.u], [p.x])
desired_order = [f.x, p.x]
lsys = ModelingToolkit.reorder_states(lsys, states(ssys), desired_order)

@test lsys.A == [-2 0; 1 -2]
@test lsys.B == reshape([1, 0], 2, 1)
@test lsys.C == [0 1]
@test lsys.D[] == 0

##
using ModelingToolkitStandardLibrary.Blocks: LimPID
k = 400
Ti = 0.5
Td = 1
Nd = 10
@named pid = LimPID(; k, Ti, Td, Nd)

@unpack reference, measurement, ctr_output = pid
lsys, ssys = linearize(pid, [reference.u, measurement.u], [ctr_output.u])
@unpack int, der = pid
desired_order = [int.x, der.x]
lsys = ModelingToolkit.reorder_states(lsys, states(ssys), desired_order)

@test lsys.A == [0 0; 0 -10]
@test lsys.B == [2 -2; 10 -10]
@test lsys.C == [400 -4000]
@test lsys.D == [4400 -4400]

# Test with the reverse desired state order as well to verify that similarity transform and reoreder_states really works
lsys = ModelingToolkit.reorder_states(lsys, states(ssys), reverse(desired_order))

@test lsys.A == [-10 0; 0 0]
@test lsys.B == [10 -10; 2 -2]
@test lsys.C == [-4000 400]
@test lsys.D == [4400 -4400]

## Test that there is a warning when input is misspecified
if VERSION >= v"1.8"
    @test_throws "Some specified inputs were not found" linearize(pid,
                                                                  [
                                                                      pid.reference.u,
                                                                      pid.measurement.u,
                                                                  ], [ctr_output.u])
    @test_throws "Some specified outputs were not found" linearize(pid,
                                                                   [
                                                                       reference.u,
                                                                       measurement.u,
                                                                   ],
                                                                   [pid.ctr_output.u])
else # v1.6 does not have the feature to match error message
    @test_throws ErrorException linearize(pid,
                                          [
                                              pid.reference.u,
                                              pid.measurement.u,
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
    ie = ModelingToolkit.IfElse.ifelse
    eqs = [
        # The equation below is equivalent to y ~ clamp(u, y_min, y_max)
        y ~ ie(u > y_max, y_max, ie((y_min < u) & (u < y_max), u, y_min)),
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

# outside the linear region the derivative is 0
lsys, ssys = linearize(sat, [u], [y]; op = Dict(u => 2))
@test isempty(lsys.A) # there are no differential variables in this system
@test isempty(lsys.B)
@test isempty(lsys.C)
@test lsys.D[] == 0
