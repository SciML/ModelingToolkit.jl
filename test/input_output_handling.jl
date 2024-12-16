using ModelingToolkit, Symbolics, Test
using ModelingToolkit: get_namespace, has_var, inputs, outputs, is_bound, bound_inputs,
                       unbound_inputs, bound_outputs, unbound_outputs, isinput, isoutput,
                       ExtraVariablesSystemException
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables xx(t) some_input(t) [input = true]
eqs = [D(xx) ~ some_input]
@named model = ODESystem(eqs, t)
@test_throws ExtraVariablesSystemException structural_simplify(model, ((), ()))
if VERSION >= v"1.8"
    err = "In particular, the unset input(s) are:\n some_input(t)"
    @test_throws err structural_simplify(model, ((), ()))
end

# Test input handling
@variables x(t) u(t) [input = true] v(t)[1:2] [input = true]
@test isinput(u)

@named sys = ODESystem([D(x) ~ -x + u], t) # both u and x are unbound
@named sys1 = ODESystem([D(x) ~ -x + v[1] + v[2]], t) # both v and x are unbound
@named sys2 = ODESystem([D(x) ~ -sys.x], t, systems = [sys]) # this binds sys.x in the context of sys2, sys2.x is still unbound
@named sys21 = ODESystem([D(x) ~ -sys1.x], t, systems = [sys1]) # this binds sys.x in the context of sys2, sys2.x is still unbound
@named sys3 = ODESystem([D(x) ~ -sys.x + sys.u], t, systems = [sys]) # This binds both sys.x and sys.u
@named sys31 = ODESystem([D(x) ~ -sys1.x + sys1.v[1]], t, systems = [sys1]) # This binds both sys.x and sys1.v[1]

@named sys4 = ODESystem([D(x) ~ -sys.x, u ~ sys.u], t, systems = [sys]) # This binds both sys.x and sys3.u, this system is one layer deeper than the previous. u is directly forwarded to sys.u, and in this case sys.u is bound while u is not

@test has_var(x ~ 1, x)
@test has_var(1 ~ x, x)
@test has_var(x + x, x)
@test !has_var(2 ~ 1, x)

@test get_namespace(x) == ""
@test get_namespace(sys.x) == "sys"
@test get_namespace(sys2.x) == "sys2"
@test get_namespace(sys2.sys.x) == "sys2₊sys"
@test get_namespace(sys21.sys1.v) == "sys21₊sys1"

@test !is_bound(sys, u)
@test !is_bound(sys, x)
@test !is_bound(sys, sys.u)
@test is_bound(sys2, sys.x)
@test !is_bound(sys2, sys.u)
@test !is_bound(sys2, sys2.sys.u)
@test is_bound(sys21, sys1.x)
@test !is_bound(sys21, sys1.v[1])
@test !is_bound(sys21, sys1.v[2])
@test is_bound(sys31, sys1.v[1])
@test !is_bound(sys31, sys1.v[2])

# simplification turns input variables into parameters
ssys, _ = structural_simplify(sys, ([u], []))
@test ModelingToolkit.isparameter(unbound_inputs(ssys)[])
@test !is_bound(ssys, u)
@test u ∈ Set(unbound_inputs(ssys))

fsys2 = flatten(sys2)
@test is_bound(fsys2, sys.x)
@test !is_bound(fsys2, sys.u)
@test !is_bound(fsys2, sys2.sys.u)

@test is_bound(sys3, sys.u) # I would like to write sys3.sys.u here but that's not how the variable is stored in the equations
@test is_bound(sys3, sys.x)

@test is_bound(sys4, sys.u)
@test !is_bound(sys4, u)

fsys4 = flatten(sys4)
@test is_bound(fsys4, sys.u)
@test !is_bound(fsys4, u)

@test isequal(inputs(sys), [u])
@test isequal(inputs(sys2), [sys.u])

@test isempty(bound_inputs(sys))
@test isequal(unbound_inputs(sys), [u])

@test isempty(bound_inputs(sys2))
@test isempty(bound_inputs(fsys2))
@test isequal(unbound_inputs(sys2), [sys.u])
@test isequal(unbound_inputs(fsys2), [sys.u])

@test isequal(bound_inputs(sys3), [sys.u])
@test isempty(unbound_inputs(sys3))

# Test output handling
@variables x(t) y(t) [output = true]
@test isoutput(y)
@named sys = ODESystem([D(x) ~ -x, y ~ x], t) # both y and x are unbound
syss = structural_simplify(sys) # This makes y an observed variable

@named sys2 = ODESystem([D(x) ~ -sys.x, y ~ sys.y], t, systems = [sys])

@test !is_bound(sys, y)
@test !is_bound(sys, x)
@test !is_bound(sys, sys.y)

@test !is_bound(syss, y)
@test !is_bound(syss, x)
@test !is_bound(syss, sys.y)

@test isequal(unbound_outputs(sys), [y])
@test isequal(unbound_outputs(syss), [y])

@test isequal(unbound_outputs(sys2), [y])
@test isequal(bound_outputs(sys2), [sys.y])

syss = structural_simplify(sys2)

@test !is_bound(syss, y)
@test !is_bound(syss, x)
@test is_bound(syss, sys.y)

#@test isequal(unbound_outputs(syss), [y])
@test isequal(bound_outputs(syss), [sys.y])

using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Mechanical.Rotational
@named inertia1 = Inertia(; J = 1)
@named inertia2 = Inertia(; J = 1)
@named spring = Rotational.Spring(; c = 10)
@named damper = Rotational.Damper(; d = 3)
@named torque = Torque(; use_support = false)
@variables y(t) = 0
eqs = [connect(torque.flange, inertia1.flange_a)
       connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
       connect(inertia2.flange_a, spring.flange_b, damper.flange_b)
       y ~ inertia2.w + torque.tau.u]
model = ODESystem(eqs, t; systems = [torque, inertia1, inertia2, spring, damper],
    name = :name)
model_outputs = [inertia1.w, inertia2.w, inertia1.phi, inertia2.phi]
model_inputs = [torque.tau.u]
matrices, ssys = linearize(model, model_inputs, model_outputs);
@test length(ModelingToolkit.outputs(ssys)) == 4

if VERSION >= v"1.8" # :opaque_closure not supported before
    let # Just to have a local scope for D
        matrices, ssys = linearize(model, model_inputs, [y])
        A, B, C, D = matrices
        obsf = ModelingToolkit.build_explicit_observed_function(ssys,
            [y],
            inputs = [torque.tau.u])
        x = randn(size(A, 1))
        u = randn(size(B, 2))
        p = (getindex.(
            Ref(ModelingToolkit.defaults_and_guesses(ssys)),
            parameters(ssys)),)
        y1 = obsf(x, u, p, 0)
        y2 = C * x + D * u
        @test y1[] ≈ y2[]
    end
end

## Code generation with unbound inputs
@testset "generate_control_function with disturbance inputs" begin
    for split in [true, false]
        simplify = true

        @variables x(t)=0 u(t)=0 [input = true]
        eqs = [
            D(x) ~ -x + u
        ]

        @named sys = ODESystem(eqs, t)
        f, dvs, ps, io_sys = ModelingToolkit.generate_control_function(sys; simplify, split)

        @test isequal(dvs[], x)
        @test isempty(ps)

        p = nothing
        x = [rand()]
        u = [rand()]
        @test f[1](x, u, p, 1) == -x + u

        # With disturbance inputs
        @variables x(t)=0 u(t)=0 [input = true] d(t)=0
        eqs = [
            D(x) ~ -x + u + d^2
        ]

        @named sys = ODESystem(eqs, t)
        f, dvs, ps, io_sys = ModelingToolkit.generate_control_function(
            sys, [u], [d]; simplify, split)

        @test isequal(dvs[], x)
        @test isempty(ps)

        p = nothing
        x = [rand()]
        u = [rand()]
        @test f[1](x, u, p, 1) == -x + u

        ## With added d argument
        @variables x(t)=0 u(t)=0 [input = true] d(t)=0
        eqs = [
            D(x) ~ -x + u + d^2
        ]

        @named sys = ODESystem(eqs, t)
        f, dvs, ps, io_sys = ModelingToolkit.generate_control_function(
            sys, [u], [d]; simplify, split, disturbance_argument = true)

        @test isequal(dvs[], x)
        @test isempty(ps)

        p = nothing
        x = [rand()]
        u = [rand()]
        d = [rand()]
        @test f[1](x, u, p, t, d) == -x + u + [d[]^2]
    end
end

## more complicated system

@variables u(t) [input = true]

function Mass(; name, m = 1.0, p = 0, v = 0)
    @variables y(t)=0 [output = true]
    ps = @parameters m = m
    sts = @variables pos(t)=p vel(t)=v
    eqs = [D(pos) ~ vel
           y ~ pos]
    ODESystem(eqs, t, [pos, vel, y], ps; name)
end

function MySpring(; name, k = 1e4)
    ps = @parameters k = k
    @variables x(t) = 0 # Spring deflection
    ODESystem(Equation[], t, [x], ps; name)
end

function MyDamper(; name, c = 10)
    ps = @parameters c = c
    @variables vel(t) = 0
    ODESystem(Equation[], t, [vel], ps; name)
end

function SpringDamper(; name, k = false, c = false)
    spring = MySpring(; name = :spring, k)
    damper = MyDamper(; name = :damper, c)
    compose(ODESystem(Equation[], t; name),
        spring, damper)
end

connect_sd(sd, m1, m2) = [sd.spring.x ~ m1.pos - m2.pos, sd.damper.vel ~ m1.vel - m2.vel]
sd_force(sd) = -sd.spring.k * sd.spring.x - sd.damper.c * sd.damper.vel

# Parameters
m1 = 1
m2 = 1
k = 1000
c = 10

@named mass1 = Mass(; m = m1)
@named mass2 = Mass(; m = m2)
@named sd = SpringDamper(; k, c)

eqs = [connect_sd(sd, mass1, mass2)
       D(mass1.vel) ~ (sd_force(sd) + u) / mass1.m
       D(mass2.vel) ~ (-sd_force(sd)) / mass2.m]
@named _model = ODESystem(eqs, t)
@named model = compose(_model, mass1, mass2, sd);

f, dvs, ps, io_sys = ModelingToolkit.generate_control_function(model, simplify = true)
@test length(dvs) == 4
@test length(ps) == length(parameters(model))
p = MTKParameters(io_sys, [io_sys.u => NaN])
x = ModelingToolkit.varmap_to_vars(
    merge(ModelingToolkit.defaults(model),
        Dict(D.(unknowns(model)) .=> 0.0)), dvs)
u = [rand()]
out = f[1](x, u, p, 1)
i = findfirst(isequal(u[1]), out)
@test i isa Int
@test iszero(out[[1:(i - 1); (i + 1):end]])

@variables x(t) u(t) [input = true]
eqs = [D(x) ~ u]
@named sys = ODESystem(eqs, t)
@test_nowarn structural_simplify(sys, ([u], []))

#=
## Disturbance input handling
We test that the generated disturbance dynamics is correct by calling the dynamics in two different points that differ in the disturbance state, and check that we get the same result as when we call the linearized dynamics in the same two points. The true system is linear so the linearized dynamics are exact.

The test below builds a double-mass model and adds an integrating disturbance to the input
=#

using ModelingToolkit
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks

# Parameters
m1 = 1
m2 = 1
k = 1000 # Spring stiffness
c = 10   # Damping coefficient

@named inertia1 = Rotational.Inertia(; J = m1)
@named inertia2 = Rotational.Inertia(; J = m2)
@named spring = Rotational.Spring(; c = k)
@named damper = Rotational.Damper(; d = c)
@named torque = Rotational.Torque(; use_support = false)

function SystemModel(u = nothing; name = :model)
    eqs = [connect(torque.flange, inertia1.flange_a)
           connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
           connect(inertia2.flange_a, spring.flange_b, damper.flange_b)]
    if u !== nothing
        push!(eqs, connect(torque.tau, u.output))
        return @named model = ODESystem(eqs, t;
            systems = [
                torque,
                inertia1,
                inertia2,
                spring,
                damper,
                u
            ])
    end
    ODESystem(eqs, t; systems = [torque, inertia1, inertia2, spring, damper], name)
end

model = SystemModel() # Model with load disturbance
model = complete(model)
model_outputs = [model.inertia1.w, model.inertia2.w, model.inertia1.phi, model.inertia2.phi]

@named dmodel = Blocks.StateSpace([0.0], [1.0], [1.0], [0.0]) # An integrating disturbance

@named dist = ModelingToolkit.DisturbanceModel(model.torque.tau.u, dmodel)
(f_oop, f_ip), outersys, dvs, p, io_sys = ModelingToolkit.add_input_disturbance(model, dist)

@unpack u, d = outersys
matrices, ssys = linearize(outersys, [u, d], model_outputs)

def = ModelingToolkit.defaults(outersys)

# Create a perturbation in the disturbance state
dstate = setdiff(dvs, model_outputs)[]
x_add = ModelingToolkit.varmap_to_vars(merge(Dict(dvs .=> 0), Dict(dstate => 1)), dvs)

x0 = randn(5)
x1 = copy(x0) + x_add # add disturbance state perturbation
u = randn(1)
pn = MTKParameters(io_sys, [])
xp0 = f_oop(x0, u, pn, 0)
xp1 = f_oop(x1, u, pn, 0)

@test xp0 ≈ matrices.A * x0 + matrices.B * [u; 0]
@test xp1 ≈ matrices.A * x1 + matrices.B * [u; 0]

@variables x(t)[1:3] = 0
@variables u(t)[1:2]
y₁, y₂, y₃ = x
u1, u2 = u
k₁, k₂, k₃ = 1, 1, 1
eqs = [D(y₁) ~ -k₁ * y₁ + k₃ * y₂ * y₃ + u1
       D(y₂) ~ k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2 + u2
       y₁ + y₂ + y₃ ~ 1]

@named sys = ODESystem(eqs, t)
m_inputs = [u[1], u[2]]
m_outputs = [y₂]
sys_simp, input_idxs = structural_simplify(sys, (; inputs = m_inputs, outputs = m_outputs))
@test isequal(unknowns(sys_simp), collect(x[1:2]))
@test length(input_idxs) == 2

# https://github.com/SciML/ModelingToolkit.jl/issues/1577
@named c = Constant(; k = 2)
@named gain = Gain(1;)
@named int = Integrator(; k = 1)
@named fb = Feedback(;)
@named model = ODESystem(
    [
        connect(c.output, fb.input1),
        connect(fb.input2, int.output),
        connect(fb.output, gain.input),
        connect(gain.output, int.input)
    ],
    t,
    systems = [int, gain, c, fb])
sys = structural_simplify(model)
@test length(unknowns(sys)) == length(equations(sys)) == 1

## Disturbance models when plant has multiple inputs
using ModelingToolkit, LinearAlgebra
using ModelingToolkit: DisturbanceModel, io_preprocessing, get_iv, get_disturbance_system
using ModelingToolkitStandardLibrary.Blocks
A, C = [randn(2, 2) for i in 1:2]
B = [1.0 0; 0 1.0]
@named model = Blocks.StateSpace(A, B, C)
@named integrator = Blocks.StateSpace([-0.001;;], [1.0;;], [1.0;;], [0.0;;])

ins = collect(complete(model).input.u)
outs = collect(complete(model).output.u)

disturbed_input = ins[1]
@named dist_integ = DisturbanceModel(disturbed_input, integrator)

(f_oop, f_ip), augmented_sys, dvs, p = ModelingToolkit.add_input_disturbance(model,
    dist_integ,
    ins)

augmented_sys = complete(augmented_sys)
matrices, ssys = linearize(augmented_sys,
    [
        augmented_sys.u,
        augmented_sys.input.u[2],
        augmented_sys.d
    ], outs)
@test matrices.A ≈ [A [1; 0]; zeros(1, 2) -0.001]
@test matrices.B == I
@test matrices.C == [C zeros(2)]
@test matrices.D == zeros(2, 3)

# Verify using ControlSystemsBase
# P = ss(A,B,C,0)
# G = ss(matrices...)
# @test sminreal(G[1, 3]) ≈ sminreal(P[1,1])*dist

@testset "Observed functions with inputs" begin
    @variables x(t)=0 u(t)=0 [input = true]
    eqs = [
        D(x) ~ -x + u
    ]

    @named sys = ODESystem(eqs, t)
    (; io_sys,) = ModelingToolkit.generate_control_function(sys, simplify = true)
    obsfn = ModelingToolkit.build_explicit_observed_function(
        io_sys, [x + u * t]; inputs = [u])
    @test obsfn([1.0], [2.0], nothing, 3.0) == [7.0]
end

# https://github.com/SciML/ModelingToolkit.jl/issues/2896
@testset "Constants substitution" begin
    @constants c = 2.0
    @variables x(t)
    eqs = [D(x) ~ c * x]
    @named sys = ODESystem(eqs, t, [x], [])

    f, dvs, ps = ModelingToolkit.generate_control_function(sys, simplify = true)
    @test f[1]([0.5], nothing, nothing, 0.0) == [1.0]
end

@testset "With callable symbolic" begin
    @variables x(t)=0 u(t)=0 [input = true]
    @parameters p(::Real) = (x -> 2x)
    eqs = [D(x) ~ -x + p(u)]
    @named sys = ODESystem(eqs, t)
    f, dvs, ps, io_sys = ModelingToolkit.generate_control_function(sys, simplify = true)
    p = MTKParameters(io_sys, [])
    u = [1.0]
    x = [1.0]
    @test_nowarn f[1](x, u, p, 0.0)
end
