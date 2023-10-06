using ModelingToolkit, Test
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq

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

@parameters t
D = Differential(t)

get_value(data, t, dt) = data[round(Int, t / dt + 1)]
@register_symbolic get_value(data, t, dt)

function Sampled(; name, data = Float64[], dt = 0.0)
    pars = @parameters begin
        data = data
        dt = dt
    end

    vars = []
    systems = @named begin
        output = RealOutput()
    end

    eqs = [
        output.u ~ get_value(data, t, dt),
    ]

    return ODESystem(eqs, t, vars, pars; name, systems,
        defaults = [output.u => data[1]])
end

vars = @variables y(t)=1 dy(t)=0 ddy(t)=0
@named src = Sampled(; data = Float64[], dt)
@named int = Integrator()

eqs = [y ~ src.output.u
    D(y) ~ dy
    D(dy) ~ ddy
    connect(src.output, int.input)]

@named sys = ODESystem(eqs, t, vars, []; systems = [int, src])
s = complete(sys)
sys = structural_simplify(sys)
prob = ODEProblem(sys, [], (0.0, t_end), [s.src.data => x]; tofloat = false)
@test prob.p isa Tuple{Vector{Float64}, Vector{Int}, Vector{Vector{Float64}}}
sol = solve(prob, ImplicitEuler());
@test sol.retcode == ReturnCode.Success
@test sol[y][end] == x[end]

#TODO: remake becomes more complicated now, how to improve?
defs = ModelingToolkit.defaults(sys)
defs[s.src.data] = 2x
p′ = ModelingToolkit.varmap_to_vars(defs, parameters(sys); tofloat = false)
p′, = ModelingToolkit.split_parameters_by_type(p′) #NOTE: we need to ensure this is called now before calling remake()
prob′ = remake(prob; p = p′)
sol = solve(prob′, ImplicitEuler());
@test sol.retcode == ReturnCode.Success
@test sol[y][end] == 2x[end]

prob′′ = remake(prob; p = [s.src.data => x])
@test_broken prob′′.p isa Tuple

# ------------------------ Mixed Type Converted to float (default behavior)

vars = @variables y(t)=1 dy(t)=0 ddy(t)=0
pars = @parameters a=1.0 b=2.0 c=3
eqs = [D(y) ~ dy * a
    D(dy) ~ ddy * b
    ddy ~ sin(t) * c]

@named model = ODESystem(eqs, t, vars, pars)
sys = structural_simplify(model)

tspan = (0.0, t_end)
prob = ODEProblem(sys, [], tspan, [])

@test prob.p isa Tuple{Vector{Float64}, Vector{Int}}
sol = solve(prob, ImplicitEuler());
@test sol.retcode == ReturnCode.Success

# ------------------------ Mixed Type Conserved

prob = ODEProblem(sys, [], tspan, []; tofloat = false)

@test prob.p isa Tuple{Vector{Float64}, Vector{Int64}}
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
    ODESystem(Equation[], ModelingToolkit.get_iv(sys), systems = [sys], name = :a_wrapper)
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
        return @named model = ODESystem(eqs,
            t;
            systems = [torque, inertia1, inertia2, spring, damper, u])
    end
    ODESystem(eqs, t; systems = [torque, inertia1, inertia2, spring, damper], name)
end

model = SystemModel() # Model with load disturbance
@named d = Step(start_time = 1.0, duration = 10.0, offset = 0.0, height = 1.0) # Disturbance
model_outputs = [model.inertia1.w, model.inertia2.w, model.inertia1.phi, model.inertia2.phi] # This is the state realization we want to control
inputs = [model.torque.tau.u]
matrices, ssys = ModelingToolkit.linearize(wr(model), inputs, model_outputs)

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
#     compose(ODESystem(eqs, t, [], []; name = name), [input, output])
# end

@named state_feedback = MatrixGain(K = -L) # Build negative feedback into the feedback matrix
@named add = Add(; k1 = 1.0, k2 = 1.0) # To add the control signal and the disturbance

connections = [[state_feedback.input.u[i] ~ model_outputs[i] for i in 1:4]
    connect(d.output, :d, add.input1)
    connect(add.input2, state_feedback.output)
    connect(add.output, :u, model.torque.tau)]
@named closed_loop = ODESystem(connections, t, systems = [model, state_feedback, add, d])
S = get_sensitivity(closed_loop, :u)

@variables t
@parameters a b c
@named s = ODESystem(Equation[], t, [], [a, b, c])
prob = ODEProblem(s, nothing, (0.0, 1.0), Pair[a => 1, b => 2, c => 3.0])
@test prob.p == ([3.0], [1, 2])
@test prob.p isa Tuple{Vector{Float64}, Vector{Int}}
