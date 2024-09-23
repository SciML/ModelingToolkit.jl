# Example: Compartmental models in epidemiology
#=
- https://github.com/epirecipes/sir-julia/blob/master/markdown/function_map/function_map.md
- https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#Deterministic_versus_stochastic_epidemic_models
=#
using ModelingToolkit, Test
using ModelingToolkit: t_nounits as t
using ModelingToolkit: get_metadata, MTKParameters

# Make sure positive shifts error
@variables x(t)
k = ShiftIndex(t)
@test_throws ErrorException @mtkbuild sys = DiscreteSystem([x(k + 1) ~ x + x(k - 1)], t)

@inline function rate_to_proportion(r, t)
    1 - exp(-r * t)
end;

# Independent and dependent variables and parameters
@parameters c nsteps δt β γ
@constants h = 1
@variables S(t) I(t) R(t)
infection = rate_to_proportion(
    β * c * I(k - 1) / (S(k - 1) * h + I(k - 1) + R(k - 1)), δt * h) * S(k - 1)
recovery = rate_to_proportion(γ * h, δt) * I(k - 1)

# Equations
eqs = [S ~ S(k - 1) - infection * h,
    I ~ I(k - 1) + infection - recovery,
    R ~ R(k - 1) + recovery]

# System
@named sys = DiscreteSystem(eqs, t, [S, I, R], [c, nsteps, δt, β, γ])
syss = structural_simplify(sys)
@test syss == syss

for df in [
    DiscreteFunction(syss),
    eval(DiscreteFunctionExpr(syss))
]

    # iip
    du = zeros(3)
    u = collect(1:3)
    p = MTKParameters(syss, [c, nsteps, δt, β, γ] .=> collect(1:5))
    df.f(du, u, p, 0)
    @test du ≈ [0.01831563888873422, 0.9816849729159067, 4.999999388195359]

    # oop
    @test df.f(u, p, 0) ≈ [0.01831563888873422, 0.9816849729159067, 4.999999388195359]
end

# Problem
u0 = [S(k - 1) => 990.0, I(k - 1) => 10.0, R(k - 1) => 0.0]
p = [β => 0.05, c => 10.0, γ => 0.25, δt => 0.1, nsteps => 400]
tspan = (0.0, ModelingToolkit.value(substitute(nsteps, p))) # value function (from Symbolics) is used to convert a Num to Float64
prob_map = DiscreteProblem(syss, u0, tspan, p)
@test prob_map.f.sys === syss

# Solution
using OrdinaryDiffEq
sol_map = solve(prob_map, FunctionMap());
@test sol_map[S] isa Vector

# Using defaults constructor
@parameters c=10.0 nsteps=400 δt=0.1 β=0.05 γ=0.25
@variables S(t)=990.0 I(t)=10.0 R(t)=0.0 R2(t)

infection2 = rate_to_proportion(β * c * I(k - 1) / (S(k - 1) + I(k - 1) + R(k - 1)), δt) *
             S(k - 1)
recovery2 = rate_to_proportion(γ, δt) * I(k - 1)

eqs2 = [S ~ S(k - 1) - infection2,
    I ~ I(k - 1) + infection2 - recovery2,
    R ~ R(k - 1) + recovery2,
    R2 ~ R]

@mtkbuild sys = DiscreteSystem(
    eqs2, t, [S, I, R, R2], [c, nsteps, δt, β, γ]; controls = [β, γ], tspan)
@test ModelingToolkit.defaults(sys) != Dict()

prob_map2 = DiscreteProblem(sys)
sol_map2 = solve(prob_map2, FunctionMap());

@test sol_map.u ≈ sol_map2.u
@test sol_map.prob.p == sol_map2.prob.p
@test_throws Any sol_map2[R2]
@test sol_map2[R2(k + 1)][begin:(end - 1)] == sol_map2[R][(begin + 1):end]
# Direct Implementation

function sir_map!(u_diff, u, p, t)
    (S, I, R) = u
    (β, c, γ, δt) = p
    N = S + I + R
    infection = rate_to_proportion(β * c * I / N, δt) * S
    recovery = rate_to_proportion(γ, δt) * I
    @inbounds begin
        u_diff[1] = S - infection
        u_diff[2] = I + infection - recovery
        u_diff[3] = R + recovery
    end
    nothing
end;
u0 = prob_map2.u0;
p = [0.05, 10.0, 0.25, 0.1];
prob_map = DiscreteProblem(sir_map!, u0, tspan, p);
sol_map2 = solve(prob_map, FunctionMap());

@test Array(sol_map) ≈ Array(sol_map2)

# Delayed difference equation
# @variables x(..) y(..) z(t)
# D1 = Difference(t; dt = 1.5)
# D2 = Difference(t; dt = 2)

# @test ModelingToolkit.is_delay_var(Symbolics.value(t), Symbolics.value(x(t - 2)))
# @test ModelingToolkit.is_delay_var(Symbolics.value(t), Symbolics.value(y(t - 1)))
# @test !ModelingToolkit.is_delay_var(Symbolics.value(t), Symbolics.value(z))
# @test_throws ErrorException ModelingToolkit.get_delay_val(Symbolics.value(t),
#     Symbolics.arguments(Symbolics.value(x(t +
#                                           2)))[1])
# @test_throws ErrorException z(t)

# # Equations
# eqs = [
#     D1(x(t)) ~ 0.4x(t) + 0.3x(t - 1.5) + 0.1x(t - 3),
#     D2(y(t)) ~ 0.3y(t) + 0.7y(t - 2) + 0.1z * h,
# ]

# # System
# @named sys = DiscreteSystem(eqs, t, [x(t), x(t - 1.5), x(t - 3), y(t), y(t - 2), z], [])

# eqs2, max_delay = ModelingToolkit.linearize_eqs(sys; return_max_delay = true)

# @test max_delay[Symbolics.operation(Symbolics.value(x(t)))] ≈ 3
# @test max_delay[Symbolics.operation(Symbolics.value(y(t)))] ≈ 2

# linearized_eqs = [eqs
#     x(t - 3.0) ~ x(t - 1.5)
#     x(t - 1.5) ~ x(t)
#     y(t - 2.0) ~ y(t)]
# @test all(eqs2 .== linearized_eqs)

# observed variable handling
@variables x(t) RHS(t)
@parameters τ
@named fol = DiscreteSystem(
    [x ~ (1 - x(k - 1)) / τ], t, [x, RHS], [τ]; observed = [RHS ~ (1 - x) / τ * h])
@test isequal(RHS, @nonamespace fol.RHS)
RHS2 = RHS
@unpack RHS = fol
@test isequal(RHS, RHS2)

# @testset "Preface tests" begin
#     using OrdinaryDiffEq
#     using Symbolics
#     using DiffEqBase: isinplace
#     using ModelingToolkit
#     using SymbolicUtils.Code
#     using SymbolicUtils: Sym

#     c = [0]
#     f = function f(c, d::Vector{Float64}, u::Vector{Float64}, p, t::Float64, dt::Float64)
#         c .= [c[1] + 1]
#         d .= randn(length(u))
#         nothing
#     end

#     dummy_identity(x, _) = x
#     @register_symbolic dummy_identity(x, y)

#     u0 = ones(5)
#     p0 = Float64[]
#     syms = [Symbol(:a, i) for i in 1:5]
#     syms_p = Symbol[]
#     dt = 0.1
#     @assert isinplace(f, 6)
#     wf = let c = c, buffer = similar(u0), u = similar(u0), p = similar(p0), dt = dt
#         t -> (f(c, buffer, u, p, t, dt); buffer)
#     end

#     num = hash(f) ⊻ length(u0) ⊻ length(p0)
#     buffername = Symbol(:fmi_buffer_, num)

#     Δ = DiscreteUpdate(t; dt = dt)
#     us = map(s -> (@variables $s(t))[1], syms)
#     ps = map(s -> (@variables $s(t))[1], syms_p)
#     buffer, = @variables $buffername[1:length(u0)]
#     dummy_var = Sym{Any}(:_) # this is safe because _ cannot be a rvalue in Julia

#     ss = Iterators.flatten((us, ps))
#     vv = Iterators.flatten((u0, p0))
#     defs = Dict{Any, Any}(s => v for (s, v) in zip(ss, vv))

#     preface = [Assignment(dummy_var, SetArray(true, term(getfield, wf, Meta.quot(:u)), us))
#         Assignment(dummy_var, SetArray(true, term(getfield, wf, Meta.quot(:p)), ps))
#         Assignment(buffer, term(wf, t))]
#     eqs = map(1:length(us)) do i
#         Δ(us[i]) ~ dummy_identity(buffer[i], us[i])
#     end

#     @mtkbuild sys = DiscreteSystem(eqs, t, us, ps; defaults = defs, preface = preface)
#     prob = DiscreteProblem(sys, [], (0.0, 1.0))
#     sol = solve(prob, FunctionMap(); dt = dt)
#     @test c[1] + 1 == length(sol)
# end

@variables x(t) y(t)
testdict = Dict([:test => 1])
@named sys = DiscreteSystem([x(k + 1) ~ 1.0], t, [x], []; metadata = testdict)
@test get_metadata(sys) == testdict

@variables x(t) y(t) u(t)
eqs = [u ~ 1
       x ~ x(k - 1) + u
       y ~ x + u]
@mtkbuild de = DiscreteSystem(eqs, t)
prob = DiscreteProblem(de, [x(k - 1) => 0.0], (0, 10))
sol = solve(prob, FunctionMap())

@test reduce(vcat, sol.u) == 1:11

# test that default values apply to the entire history
@variables x(t) = 1.0
@mtkbuild de = DiscreteSystem([x ~ x(k - 1) + x(k - 2)], t)
prob = DiscreteProblem(de, [], (0, 10))
@test prob[x] == 2.0
@test prob[x(k - 1)] == 1.0

# must provide initial conditions for history
@test_throws ErrorException DiscreteProblem(de, [x => 2.0], (0, 10))

# initial values only affect _that timestep_, not the entire history
prob = DiscreteProblem(de, [x(k - 1) => 2.0], (0, 10))
@test prob[x] == 3.0
@test prob[x(k - 1)] == 2.0

# Issue#2585
getdata(buffer, t) = buffer[mod1(Int(t), length(buffer))]
@register_symbolic getdata(buffer::Vector, t)
k = ShiftIndex(t)
function SampledData(; name, buffer)
    L = length(buffer)
    pars = @parameters begin
        buffer[1:L] = buffer
    end
    @variables output(t) time(t)
    eqs = [time ~ time(k - 1) + 1
           output ~ getdata(buffer, time)]
    return DiscreteSystem(eqs, t; name)
end
function System(; name, buffer)
    @named y_sys = SampledData(; buffer = buffer)
    pars = @parameters begin
        α = 0.5, [description = "alpha"]
        β = 0.5, [description = "beta"]
    end
    vars = @variables y(t)=0.0 y_shk(t)=0.0

    eqs = [y_shk ~ y_sys.output
           # y[t] = 0.5 * y[t - 1] + 0.5 * y[t + 1] + y_shk[t]
           y(k - 1) ~ α * y(k - 2) + (β * y(k) + y_shk(k - 1))]

    DiscreteSystem(eqs, t, vars, pars; systems = [y_sys], name = name)
end

@test_nowarn @mtkbuild sys = System(; buffer = ones(10))

# Ensure discrete systems with algebraic equations throw
@variables x(t) y(t)
k = ShiftIndex(t)
@named sys = DiscreteSystem([x ~ x^2 + y^2, y ~ x(k - 1) + y(k - 1)], t)
@test_throws ["algebraic equations", "not yet supported"] structural_simplify(sys)

@testset "Passing `nothing` to `u0`" begin
    @variables x(t) = 1
    k = ShiftIndex()
    @mtkbuild sys = DiscreteSystem([x(k) ~ x(k - 1) + 1], t)
    prob = @test_nowarn DiscreteProblem(sys, nothing, (0.0, 1.0))
    @test_nowarn solve(prob, FunctionMap())
end
