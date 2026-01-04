using ModelingToolkitBase, Symbolics, Test
using ModelingToolkitBase: get_namespace, has_var, inputs, outputs, is_bound, bound_inputs,
    unbound_inputs, bound_outputs, unbound_outputs, isinput, isoutput,
    ExtraVariablesSystemException
using NonlinearSolve
using SymbolicIndexingInterface: is_parameter
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@variables xx(t) some_input(t) [input = true]
eqs = [D(xx) ~ some_input]
@named model = System(eqs, t)
err = @isdefined(ModelingToolkit) ? ModelingToolkit.StateSelection.ExtraVariablesSystemException : ExtraVariablesSystemException
@test_throws err mtkcompile(model)
if @isdefined(ModelingToolkit)
    err = "In particular, the unset input(s) are:\n some_input(t)"
    @test_throws err mtkcompile(model)
end

# Test input handling
@variables x(t) u(t) [input = true] v(t)[1:2] [input = true]
@test isinput(u)

@named sys = System([D(x) ~ -x + u], t) # both u and x are unbound
@named sys1 = System([D(x) ~ -x + v[1] + v[2]], t) # both v and x are unbound
@named sys2 = System([D(x) ~ -sys.x], t, systems = [sys]) # this binds sys.x in the context of sys2, sys2.x is still unbound
@named sys21 = System([D(x) ~ -sys1.x], t, systems = [sys1]) # this binds sys.x in the context of sys2, sys2.x is still unbound
@named sys3 = System([D(x) ~ -sys.x + sys.u], t, systems = [sys]) # This binds both sys.x and sys.u
@named sys31 = System([D(x) ~ -sys1.x + sys1.v[1]], t, systems = [sys1]) # This binds both sys.x and sys1.v[1]

@named sys4 = System([D(x) ~ -sys.x, u ~ sys.u], t, systems = [sys]) # This binds both sys.x and sys3.u, this system is one layer deeper than the previous. u is directly forwarded to sys.u, and in this case sys.u is bound while u is not

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
ssys = mtkcompile(sys, inputs = [u], outputs = [])
@test is_parameter(ssys, unbound_inputs(ssys)[])
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
@named sys = System([D(x) ~ -x, y ~ x], t) # both y and x are unbound
syss = mtkcompile(sys, outputs = [y]) # This makes y an observed variable

@named sys2 = System([D(x) ~ -sys.x, y ~ sys.y], t, systems = [sys])

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

syss = mtkcompile(sys2, outputs = [sys.y])

@test !is_bound(syss, y)
@test !is_bound(syss, x)
@test is_bound(syss, sys.y)

#@test isequal(unbound_outputs(syss), [y])
@test isequal(bound_outputs(syss), [sys.y])

if @isdefined(ModelingToolkit)
    using ModelingToolkitStandardLibrary
    using ModelingToolkitStandardLibrary.Mechanical.Rotational
    @named inertia1 = Inertia(; J = 1)
    @named inertia2 = Inertia(; J = 1)
    @named spring = Rotational.Spring(; c = 10)
    @named damper = Rotational.Damper(; d = 3)
    @named torque = Torque(; use_support = false)
    @variables y(t) = 0
    eqs = [
        connect(torque.flange, inertia1.flange_a)
        connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
        connect(inertia2.flange_a, spring.flange_b, damper.flange_b)
        y ~ inertia2.w + torque.tau.u
    ]
    model = System(
        eqs, t; systems = [torque, inertia1, inertia2, spring, damper],
        name = :name, guesses = [spring.flange_a.phi => 0.0]
    )
    model_outputs = [inertia1.w, inertia2.w, inertia1.phi, inertia2.phi]
    model_inputs = [torque.tau.u]
    op = Dict(torque.tau.u => 0.0)
    matrices, ssys = linearize(
        model, model_inputs, model_outputs; op
    )
    @test length(ModelingToolkit.outputs(ssys)) == 4

    let # Just to have a local scope for D
        matrices, ssys = linearize(model, model_inputs, [y]; op)
        A, B, C, D = matrices
        obsf = ModelingToolkit.build_explicit_observed_function(
            ssys,
            [y],
            inputs = [torque.tau.u]
        )
        x = randn(size(A, 1))
        u = randn(size(B, 2))
        p = getindex.(
            Ref(ModelingToolkit.initial_conditions_and_guesses(ssys)),
            parameters(ssys)
        )
        y1 = obsf(x, u, p, 0)
        y2 = C * x + D * u
        @test y1[] ≈ y2[]
    end
end

## Code generation with unbound inputs
@testset "generate_control_function with disturbance inputs" begin
    for split in [true, false]
        simplify = true

        @variables x(t) = 0 u(t) = 0 [input = true]
        eqs = [
            D(x) ~ -x + u,
        ]

        @named sys = System(eqs, t)
        f, dvs,
            ps, io_sys = ModelingToolkitBase.generate_control_function(
            sys, [u]; simplify, split
        )

        @test isequal(dvs[], x)
        @test isempty(ps)

        p = [rand()]
        x = [rand()]
        u = [rand()]
        @test f[1](x, u, p, 1) ≈ -x + u

        # With disturbance inputs
        @variables x(t) = 0 u(t) = 0 [input = true] d(t) = 0
        eqs = [
            D(x) ~ -x + u + d^2,
        ]

        @named sys = System(eqs, t)
        f, dvs,
            ps,
            io_sys = ModelingToolkitBase.generate_control_function(
            sys, [u], [d]; simplify, split
        )

        @test isequal(dvs[], x)
        @test isempty(ps)

        p = [rand()]
        x = [rand()]
        u = [rand()]
        @test f[1](x, u, p, 1) ≈ -x + u

        ## With added d argument
        @variables x(t) = 0 u(t) = 0 [input = true] d(t) = 0
        eqs = [
            D(x) ~ -x + u + d^2,
        ]

        @named sys = System(eqs, t)
        f, dvs,
            ps,
            io_sys = ModelingToolkitBase.generate_control_function(
            sys, [u], [d];
            simplify, split, disturbance_argument = true
        )

        @test isequal(dvs[], x)
        @test isempty(ps)

        p = [rand()]
        x = [rand()]
        u = [rand()]
        d = [rand()]
        @test f[1](x, u, p, t, d) ≈ -x + u + [d[]^2]

        ## Test new known_disturbance_inputs parameter (equivalent to disturbance_argument=true)
        @variables x(t) = 0 u(t) = 0 [input = true] d(t) = 0
        eqs = [
            D(x) ~ -x + u + d^2,
        ]

        @named sys = System(eqs, t)
        f, dvs,
            ps,
            io_sys = ModelingToolkitBase.generate_control_function(
            sys, [u];
            known_disturbance_inputs = [d],
            simplify, split
        )

        @test isequal(dvs[], x)
        @test isempty(ps)

        p = [rand()]
        x = [rand()]
        u = [rand()]
        d = [rand()]
        @test f[1](x, u, p, t, d) ≈ -x + u + [d[]^2]

        ## Test mixed known and unknown disturbances
        @variables x(t) = 0 u(t) = 0 [input = true] d1(t) = 0 d2(t) = 0
        eqs = [
            D(x) ~ -x + u + d1^2 + d2^3,
        ]

        @named sys = System(eqs, t)
        f, dvs,
            ps,
            io_sys = ModelingToolkitBase.generate_control_function(
            sys, [u], [d1];           # d1 is unknown (set to zero)
            known_disturbance_inputs = [d2],     # d2 is known (function argument)
            simplify, split
        )

        @test isequal(dvs[], x)
        @test isempty(ps)

        p = [rand()]
        x = [rand()]
        u = [rand()]
        d2_val = [rand()]
        # d1 is set to zero, so only u and d2 affect the dynamics
        @test f[1](x, u, p, t, d2_val) ≈ -x + u + [d2_val[]^3]
    end
end

## more complicated system

@variables u(t) [input = true]

function Mass(; name, m = 1.0, p = 0, v = 0)
    @variables y(t) = 0 [output = true]
    ps = @parameters m = m
    sts = @variables pos(t) = p vel(t) = v
    eqs = [
        D(pos) ~ vel
        y ~ pos
    ]
    return System(eqs, t, [pos, vel, y], ps; name)
end

function MySpring(; name, k = 1.0e4)
    ps = @parameters k = k
    @variables x(t) = 0 # Spring deflection
    return System(Equation[], t, [x], ps; name)
end

function MyDamper(; name, c = 10)
    ps = @parameters c = c
    @variables vel(t) = 0
    return System(Equation[], t, [vel], ps; name)
end

function SpringDamper(; name, k = false, c = false)
    spring = MySpring(; name = :spring, k)
    damper = MyDamper(; name = :damper, c)
    return compose(
        System(Equation[], t; name),
        spring, damper
    )
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

eqs = [
    connect_sd(sd, mass1, mass2)
    D(mass1.vel) ~ (sd_force(sd) + u) / mass1.m
    D(mass2.vel) ~ (-sd_force(sd)) / mass2.m
]
@named _model = System(eqs, t)
@named model = compose(_model, mass1, mass2, sd);

f, dvs, ps, io_sys = ModelingToolkitBase.generate_control_function(
    model, [u]; simplify = true
)
@test length(dvs) == 4
p = MTKParameters(io_sys, [io_sys.u => NaN])
x = ModelingToolkitBase.varmap_to_vars(
    merge(
        ModelingToolkitBase.initial_conditions(model),
        Dict(D.(unknowns(model)) .=> 0.0)
    ), dvs
)
u = [rand()]
out = f[1](x, u, p, 1)
i = findfirst(isequal(u[1]), out)
@test i isa Int
@test iszero(out[[1:(i - 1); (i + 1):end]])

@variables x(t) u(t) [input = true]
eqs = [D(x) ~ u]
@named sys = System(eqs, t)
@test_nowarn mtkcompile(sys, inputs = [u], outputs = [])

@testset "IO canonicalization" begin
    @variables x(t)[1:3] = zeros(3)
    @variables u(t)[1:2]
    y₁, y₂, y₃ = x
    u1, u2 = u
    k₁, k₂, k₃ = 1, 1, 1
    eqs = [
        D(y₁) ~ -k₁ * y₁ + k₃ * y₂ * y₃ + u1
        D(y₂) ~ k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2 + u2
        y₁ + y₂ + y₃ ~ 1
    ]

    @named sys = System(eqs, t)

    @test_throws "entire array must be an input" mtkcompile(sys; inputs = [u1])
    @test_throws "entire array must be an disturbance input" mtkcompile(sys; disturbance_inputs = [x[1]])
    @test_throws "sorted order" mtkcompile(sys; inputs = [u2, u1])
    @test_throws "sorted order" mtkcompile(sys; disturbance_inputs = [u2, u1])

    ss1 = mtkcompile(sys; inputs = [u], outputs = [x])
    @test isequal(ModelingToolkitBase.inputs(ss1), [u1, u2])
    @test isequal(ModelingToolkitBase.outputs(ss1), [x[1], x[2], x[3]])
end

using ModelingToolkitStandardLibrary.Blocks

if @isdefined(ModelingToolkit)
    @testset "Issue#1577" begin
        # https://github.com/SciML/ModelingToolkit.jl/issues/1577
        @named c = Constant(; k = 2)
        @named gain = Gain(1)
        @named int = Integrator(; k = 1)
        @named fb = Feedback()
        @named model = System(
            [
                connect(c.output, fb.input1),
                connect(fb.input2, int.output),
                connect(fb.output, gain.input),
                connect(gain.output, int.input),
            ],
            t,
            systems = [int, gain, c, fb]
        )
        sys = mtkcompile(model)
        @test length(unknowns(sys)) == length(equations(sys)) == 1
    end
end

# Verify using ControlSystemsBase
# P = ss(A,B,C,0)
# G = ss(matrices...)
# @test sminreal(G[1, 3]) ≈ sminreal(P[1,1])*dist

@testset "Observed functions with inputs" begin
    @variables x(t) = 0 u(t) = 0 [input = true]
    eqs = [
        D(x) ~ -x + u,
    ]

    @named sys = System(eqs, t)
    (; io_sys) = ModelingToolkitBase.generate_control_function(
        sys, [u]; simplify = true
    )
    obsfn = ModelingToolkitBase.build_explicit_observed_function(
        io_sys, [x + u * t]; inputs = [u]
    )
    @test obsfn([1.0], [2.0], MTKParameters(io_sys, []), 3.0) ≈ [7.0]
end

# https://github.com/SciML/ModelingToolkit.jl/issues/2896
@testset "Constants substitution" begin
    @constants c = 2.0
    @variables x(t)
    eqs = [D(x) ~ c * x]
    @mtkcompile sys = System(eqs, t, [x], [c])

    f, dvs, ps, io_sys = ModelingToolkitBase.generate_control_function(sys)
    @test f[1]([0.5], nothing, MTKParameters(io_sys, []), 0.0) ≈ [1.0]
end

@testset "With callable symbolic" begin
    @variables x(t) = 0 u(t) = 0 [input = true]
    @parameters p(::Real) = (x -> 2x)
    eqs = [D(x) ~ -x + p(u)]
    @named sys = System(eqs, t)
    f, dvs, ps, io_sys = ModelingToolkitBase.generate_control_function(sys, [u])
    p = MTKParameters(io_sys, [])
    u = [1.0]
    x = [1.0]
    @test_nowarn f[1](x, u, p, 0.0)
end

@testset "Observed inputs and outputs" begin
    @variables x(t) y(t) [input = true] z(t) [output = true]
    eqs = [
        D(x) ~ x + y + z
        y ~ z
    ]
    @named sys = System(eqs, t)
    @test issetequal(ModelingToolkitBase.inputs(sys), [y])
    @test issetequal(ModelingToolkitBase.outputs(sys), [z])

    ss1 = mtkcompile(sys, inputs = [y], outputs = [z])
    @test issetequal(ModelingToolkitBase.inputs(ss1), [y])
    @test issetequal(ModelingToolkitBase.outputs(ss1), [z])

    ss2 = mtkcompile(sys, inputs = [z], outputs = [y])
    @test issetequal(ModelingToolkitBase.inputs(ss2), [z])
    @test issetequal(ModelingToolkitBase.outputs(ss2), [y])
end

@testset "Retain inputs when composing systems" begin
    @variables x(t) y(t) [input = true]
    @named sys = System([D(x) ~ y * x], t)
    csys = compose(System(Equation[], t; name = :outer), sys)
    @test issetequal(ModelingToolkitBase.inputs(csys), [sys.y])

    # More complex hierarchy
    @variables z(t) [input = true] w(t)
    @named sys2 = System([D(w) ~ z - w], t)
    cosys = compose(System(Equation[], t; name = :outermost), [csys, sys2])
    @test issetequal(ModelingToolkitBase.inputs(cosys), [csys.sys.y, sys2.z])
end
