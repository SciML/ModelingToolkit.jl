using ModelingToolkit, Test
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq

dt = 4e-4
t_end = 10.0
time = 0:dt:t_end
x = @. time^2 + 1.0

@parameters t
D = Differential(t)

vars = @variables y(t)=1 dy(t)=0 ddy(t)=0
@named src = SampledData(; data = Float64[], dt)
@named int = Integrator()

eqs = [y ~ src.output.u
    D(y) ~ dy
    D(dy) ~ ddy
    connect(src.output, int.input)]

@named sys = ODESystem(eqs, t, vars, []; systems = [int, src])
s = complete(sys)
sys = structural_simplify(sys)

prob = ODEProblem(sys, [], (0.0, t_end), [s.src.data => x]; split_parameters = true)
@test prob.p isa Tuple{Vector{Float64}, Vector{Int}, Vector{Vector{Float64}}}
@time sol = solve(prob, ImplicitEuler());
prob2 = ODEProblem(sys, [], (0.0, t_end), [s.src.data => x]; to_float = false)
@test prob2.p isa Vector{Union{Float64, Int64, Vector{Float64}}}
@time sol2 = solve(prob2, ImplicitEuler());
