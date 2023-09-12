using ModelingToolkit, Test
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq

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
prob = ODEProblem(sys, [], (0.0, t_end), [s.src.data => x])
@test prob.p isa Tuple{Vector{Float64}, Vector{Int}, Vector{Vector{Float64}}}
sol = solve(prob, ImplicitEuler());
@test sol.retcode == ReturnCode.Success
@test sol[y][end] == x[end]

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

@test prob.p isa Vector{Float64}
sol = solve(prob, ImplicitEuler());
@test sol.retcode == ReturnCode.Success

# ------------------------ Mixed Type Conserved

prob = ODEProblem(sys, [], tspan, []; tofloat = false)

@test prob.p isa Tuple{Vector{Float64}, Vector{Int64}}
sol = solve(prob, ImplicitEuler());
@test sol.retcode == ReturnCode.Success
