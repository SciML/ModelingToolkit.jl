using ModelingToolkit, Test
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq
using DataInterpolations
using BlockArrays: BlockedArray
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: MTKParameters, ParameterIndex, NONNUMERIC_PORTION
using SciMLStructures: Tunable, Discrete, Constants, Initials
using StaticArrays: SizedVector
using SymbolicIndexingInterface: is_parameter, getp

x = [1, 2.0, false, [1, 2, 3], Parameter(1.0)]

y = ModelingToolkit.promote_to_concrete(x)
@test eltype(y) == Union{Float64, Parameter{Float64}, Vector{Int64}}

y = ModelingToolkit.promote_to_concrete(x; tofloat = false)
@test eltype(y) == Union{Bool, Float64, Int64, Parameter{Float64}, Vector{Int64}}

x = [1, 2.0, false, [1, 2, 3]]
y = ModelingToolkit.promote_to_concrete(x)
@test eltype(y) == Union{Float64, Vector{Int64}}

x = Any[1, 2.0, false]
y = ModelingToolkit.promote_to_concrete(x; tofloat = false)
@test eltype(y) == Union{Bool, Float64, Int64}

y = ModelingToolkit.promote_to_concrete(x; use_union = false)
@test eltype(y) == Float64

x = Float16[1.0, 2.0, 3.0]
y = ModelingToolkit.promote_to_concrete(x)
@test eltype(y) == Float16

# ------------------------ Mixed Single Values and Vector

dt = 4e-4
t_end = 10.0
time = 0:dt:t_end
x = @. time^2 + 1.0

struct Interpolator
    data::Vector{Float64}
    dt::Float64
end

function (i::Interpolator)(t)
    return i.data[round(Int, t / i.dt + 1)]
end
@register_symbolic (i::Interpolator)(t)

get_value(interp::Interpolator, t) = interp(t)
@register_symbolic get_value(interp::Interpolator, t)

Symbolics.derivative(::typeof(get_value), args::NTuple{2, Any}, ::Val{2}) = 0

function Sampled(; name, interp = Interpolator(Float64[], 0.0))
    pars = @parameters begin
        interpolator::Interpolator = interp
    end

    vars = []
    systems = @named begin
        output = RealOutput()
    end

    eqs = [
        output.u ~ get_value(interpolator, t)
    ]

    return System(eqs, t, vars, [interpolator]; name, systems)
end

vars = @variables y(t) dy(t) ddy(t)
@named src = Sampled(; interp = Interpolator(x, dt))
@named int = Integrator()

eqs = [y ~ src.output.u
       D(y) ~ dy
       D(dy) ~ ddy
       connect(src.output, int.input)]

@named sys = System(eqs, t, vars, []; systems = [int, src])
s = complete(sys)
sys = mtkcompile(sys)
prob = ODEProblem(
    sys, [s.src.interpolator => Interpolator(x, dt)], (0.0, t_end);
    tofloat = false)
sol = solve(prob, ImplicitEuler());
@test sol.retcode == ReturnCode.Success
@test sol[y][end] == x[end]

#TODO: remake becomes more complicated now, how to improve?
defs = ModelingToolkit.defaults(sys)
defs[s.src.interpolator] = Interpolator(2x, dt)
p′ = ModelingToolkit.MTKParameters(sys, defs)
prob′ = remake(prob; p = p′)
sol = solve(prob′, ImplicitEuler());
@test sol.retcode == ReturnCode.Success
@test sol[y][end] == 2x[end]

# ------------------------ Mixed Type Converted to float (default behavior)

vars = @variables y(t)=1 dy(t)=0 ddy(t)=0
pars = @parameters a=1.0 b=2.0 c=3
eqs = [D(y) ~ dy * a
       D(dy) ~ ddy * b
       ddy ~ sin(t) * c]

@named model = System(eqs, t, vars, pars)
sys = mtkcompile(model; split = false)

tspan = (0.0, t_end)
prob = ODEProblem(sys, [], tspan; build_initializeprob = false)

@test prob.p isa Vector{Float64}
sol = solve(prob, ImplicitEuler());
@test sol.retcode == ReturnCode.Success

# ------------------------ Mixed Type Conserved

prob = ODEProblem(
    sys, [], tspan; tofloat = false, build_initializeprob = false)

sol = solve(prob, ImplicitEuler());
@test sol.retcode == ReturnCode.Success

# ------------------------- Bug
using ModelingToolkit, LinearAlgebra
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Blocks: t
using ModelingToolkit: connect

"A wrapper function to make symbolic indexing easier"
function wr(sys)
    System(Equation[], ModelingToolkit.get_iv(sys), systems = [sys], name = :a_wrapper)
end
indexof(sym, syms) = findfirst(isequal(sym), syms)

# Parameters
m1 = 1.0
m2 = 1.0
k = 10.0 # Spring stiffness
c = 3.0  # Damping coefficient

@named inertia1 = Inertia(; J = m1)
@named inertia2 = Inertia(; J = m2)
@named spring = Spring(; c = k)
@named damper = Damper(; d = c)
@named torque = Torque(use_support = false)

function SystemModel(u = nothing; name = :model)
    eqs = [connect(torque.flange, inertia1.flange_a)
           connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
           connect(inertia2.flange_a, spring.flange_b, damper.flange_b)]
    if u !== nothing
        push!(eqs, connect(torque.tau, u.output))
        return @named model = System(eqs,
            t;
            systems = [torque, inertia1, inertia2, spring, damper, u])
    end
    System(eqs, t; systems = [torque, inertia1, inertia2, spring, damper],
        name, guesses = [spring.flange_a.phi => 0.0])
end

model = SystemModel() # Model with load disturbance
@named d = Step(start_time = 1.0, duration = 10.0, offset = 0.0, height = 1.0) # Disturbance
model_outputs = [model.inertia1.w, model.inertia2.w, model.inertia1.phi, model.inertia2.phi] # This is the state realization we want to control
inputs = [model.torque.tau.u]
op = [model.torque.tau.u => 0.0]
matrices, ssys = ModelingToolkit.linearize(
    wr(model), inputs, model_outputs; op)

# Design state-feedback gain using LQR
# Define cost matrices
x_costs = [model.inertia1.w => 1.0
           model.inertia2.w => 1.0
           model.inertia1.phi => 1.0
           model.inertia2.phi => 1.0]
L = randn(1, 4) # Post-multiply by `C` to get the correct input to the controller

# This old definition of MatrixGain will work because the parameter space does not include K (an Array term)
# @component function MatrixGainAlt(K::AbstractArray; name)
#     nout, nin = size(K, 1), size(K, 2)
#     @named input = RealInput(; nin = nin)
#     @named output = RealOutput(; nout = nout)
#     eqs = [output.u[i] ~ sum(K[i, j] * input.u[j] for j in 1:nin) for i in 1:nout]
#     compose(System(eqs, t, [], []; name = name), [input, output])
# end

@named state_feedback = MatrixGain(K = -L) # Build negative feedback into the feedback matrix
@named add = Add(; k1 = 1.0, k2 = 1.0) # To add the control signal and the disturbance

connections = [[state_feedback.input.u[i] ~ model_outputs[i] for i in 1:4]
               connect(d.output, :d, add.input1)
               connect(add.input2, state_feedback.output)
               connect(add.output, :u, model.torque.tau)]
@named closed_loop = System(connections, t, systems = [model, state_feedback, add, d])
S = get_sensitivity(closed_loop, :u)

@testset "Indexing MTKParameters with ParameterIndex" begin
    ps = MTKParameters(collect(1.0:10.0), collect(11.0:20.0),
        (BlockedArray([true, false, false, true], [2, 2]),
            BlockedArray([[1 2; 3 4], [2 4; 6 8]], [1, 1])),
        # (BlockedArray([[true, false], [false, true]]), BlockedArray([[[1 2; 3 4]], [[2 4; 6 8]]])),
        ([5, 6],),
        (["hi", "bye"], [:lie, :die]), ())
    @test ps[ParameterIndex(Tunable(), 1)] == 1.0
    @test ps[ParameterIndex(Tunable(), 2:4)] == collect(2.0:4.0)
    @test ps[ParameterIndex(Tunable(), reshape(4:7, 2, 2))] == reshape(4.0:7.0, 2, 2)
    @test ps[ParameterIndex(Initials(), 1)] == 11.0
    @test ps[ParameterIndex(Initials(), 2:4)] == collect(12.0:14.0)
    @test ps[ParameterIndex(Initials(), reshape(4:7, 2, 2))] == reshape(14.0:17.0, 2, 2)
    @test ps[ParameterIndex(Discrete(), (2, 1, 2, 2))] == 4
    @test ps[ParameterIndex(Discrete(), (2, 2))] == [2 4; 6 8]
    @test ps[ParameterIndex(Constants(), (1, 1))] == 5
    @test ps[ParameterIndex(NONNUMERIC_PORTION, (2, 2))] == :die

    ps[ParameterIndex(Tunable(), 1)] = 1.5
    ps[ParameterIndex(Tunable(), 2:4)] = [2.5, 3.5, 4.5]
    ps[ParameterIndex(Tunable(), reshape(5:8, 2, 2))] = [5.5 7.5; 6.5 8.5]
    ps[ParameterIndex(Initials(), 1)] = 11.5
    ps[ParameterIndex(Initials(), 2:4)] = [12.5, 13.5, 14.5]
    ps[ParameterIndex(Initials(), reshape(5:8, 2, 2))] = [15.5 17.5; 16.5 18.5]
    ps[ParameterIndex(Discrete(), (2, 1, 2, 2))] = 5
    @test ps[ParameterIndex(Tunable(), 1:8)] == collect(1.0:8.0) .+ 0.5
    @test ps[ParameterIndex(Initials(), 1:8)] == collect(11.0:18.0) .+ 0.5
    @test ps[ParameterIndex(Discrete(), (2, 1, 2, 2))] == 5
end

@testset "Callable parameters" begin
    @testset "As FunctionWrapper" begin
        _f1(x) = 2x
        struct Foo end
        (::Foo)(x) = 3x
        @variables x(t)
        @parameters fn(::Real) = _f1
        @mtkcompile sys = System(D(x) ~ fn(t), t)
        @test is_parameter(sys, fn)
        @test ModelingToolkit.defaults(sys)[fn] == _f1

        getter = getp(sys, fn)
        prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0))
        @inferred getter(prob)
        # cannot be inferred better since `FunctionWrapper` is only known to return `Real`
        @inferred Vector{<:Real} prob.f(prob.u0, prob.p, prob.tspan[1])
        sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10)
        @test sol.u[end][] ≈ 2.0

        prob = ODEProblem(sys, [x => 1.0, fn => Foo()], (0.0, 1.0))
        @inferred getter(prob)
        @inferred Vector{<:Real} prob.f(prob.u0, prob.p, prob.tspan[1])
        sol = solve(prob; abstol = 1e-10, reltol = 1e-10)
        @test sol.u[end][] ≈ 2.5
    end

    @testset "Concrete function type" begin
        ts = 0.0:0.1:1.0
        interp = LinearInterpolation(
            ts .^ 2, ts; extrapolation = ExtrapolationType.Extension)
        @variables x(t)
        @parameters (fn::typeof(interp))(..)
        @mtkcompile sys = System(D(x) ~ fn(x), t)
        @test is_parameter(sys, fn)
        getter = getp(sys, fn)
        prob = ODEProblem(sys, [x => 1.0, fn => interp], (0.0, 1.0))
        @inferred getter(prob)
        @inferred prob.f(prob.u0, prob.p, prob.tspan[1])
        @test_nowarn sol = solve(prob, Tsit5())
        @test_nowarn prob.ps[fn] = LinearInterpolation(
            ts .^ 3, ts; extrapolation = ExtrapolationType.Extension)
        @test_nowarn sol = solve(prob)
    end
end

@testset "" begin
    @mtkmodel SubSystem begin
        @parameters begin
            c = 1
        end
        @variables begin
            x(t)
        end
        @equations begin
            D(x) ~ c * x
        end
    end

    @mtkmodel ApexSystem begin
        @components begin
            subsys = SubSystem()
        end
        @parameters begin
            k = 1
        end
        @variables begin
            y(t)
        end
        @equations begin
            D(y) ~ k * y + subsys.x
        end
    end

    @named sys = ApexSystem()
    sysref = complete(sys)
    sys2 = complete(sys; split = true, flatten = false)
    ps = Set(full_parameters(sys2))
    @test sysref.k in ps
    @test sysref.subsys.c in ps
    @test length(ps) == 2
end
