# Example: Compartmental models in epidemiology
#=
- https://github.com/epirecipes/sir-julia/blob/master/markdown/function_map/function_map.md
- https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#Deterministic_versus_stochastic_epidemic_models
=#
using ModelingToolkit, SymbolicIndexingInterface, Test
using ModelingToolkit: t_nounits as t

# Make sure positive shifts error
@variables x(t)
k = ShiftIndex(t)
@test_throws ErrorException @mtkcompile sys = System([x(k + 1) ~ x + x(k - 1)], t)

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
@named sys = System(eqs, t, [S, I, R], [c, nsteps, δt, β, γ, h])
syss = mtkcompile(sys)
@test syss == syss

df = DiscreteFunction(syss)
# iip
du = zeros(3)
u = ModelingToolkit.varmap_to_vars(
    Dict([S(k - 1) => 1, I(k - 1) => 2, R(k - 1) => 3]), unknowns(syss))
p = MTKParameters(syss, [c, nsteps, δt, β, γ] .=> collect(1:5))
df.f(du, u, p, 0)
reorderer = getu(syss, [S(k - 1), I(k - 1), R(k - 1)])
@test reorderer(du) ≈ [0.01831563888873422, 0.9816849729159067, 4.999999388195359]

# oop
@test reorderer(df.f(u, p, 0)) ≈
      [0.01831563888873422, 0.9816849729159067, 4.999999388195359]

# Problem
u0 = [S => 990.0, I => 10.0, R => 0.0]
p = [β => 0.05, c => 10.0, γ => 0.25, δt => 0.1, nsteps => 400]
tspan = (0.0, ModelingToolkit.value(substitute(nsteps, p))) # value function (from Symbolics) is used to convert a Num to Float64
prob_map = DiscreteProblem(
    syss, [u0; p], tspan; guesses = [S(k - 1) => 1.0, I(k - 1) => 1.0, R(k - 1) => 1.0])
@test prob_map.f.sys === syss

# Solution
using OrdinaryDiffEq
sol_map = solve(prob_map, FunctionMap());
@test sol_map[S] isa Vector
@test sol_map[S(k - 1)] isa Vector

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

@mtkcompile sys = System(
    eqs2, t, [S, I, R, R2], [c, nsteps, δt, β, γ])
@test ModelingToolkit.defaults(sys) != Dict()

prob_map2 = DiscreteProblem(sys, [], tspan)
# prob_map2 = DiscreteProblem(sys, [S(k - 1) => S, I(k - 1) => I, R(k - 1) => R], tspan)
sol_map2 = solve(prob_map2, FunctionMap());

@test sol_map.u ≈ sol_map2.u
for p in parameters(sys)
    @test sol_map.prob.ps[p] ≈ sol_map2.prob.ps[p]
end
@test sol_map2[R2][begin:(end - 1)] == sol_map2[R(k - 1)][(begin + 1):end] ==
      sol_map2[R][begin:(end - 1)]
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
u0 = sol_map2[[S, I, R], 1];
p = [0.05, 10.0, 0.25, 0.1];
prob_map = DiscreteProblem(sir_map!, u0, tspan, p);
sol_map2 = solve(prob_map, FunctionMap());

@test reduce(hcat, sol_map[[S, I, R]]) ≈ Array(sol_map2)

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
# @named sys = System(eqs, t, [x(t), x(t - 1.5), x(t - 3), y(t), y(t - 2), z], [])

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
@named fol = System(
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

#     @mtkcompile sys = System(eqs, t, us, ps; defaults = defs, preface = preface)
#     prob = DiscreteProblem(sys, [], (0.0, 1.0))
#     sol = solve(prob, FunctionMap(); dt = dt)
#     @test c[1] + 1 == length(sol)
# end

@variables x(t) y(t) u(t)
eqs = [u ~ 1
       x ~ x(k - 1) + u
       y ~ x + u]
@mtkcompile de = System(eqs, t)
prob = DiscreteProblem(de, [x(k - 1) => 0.0], (0, 10))
sol = solve(prob, FunctionMap())

@test sol[x] == 1:11

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
    return System(eqs, t; name)
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

    System(eqs, t, vars, pars; systems = [y_sys], name = name)
end

@test_nowarn @mtkcompile sys = System(; buffer = ones(10))

@testset "Passing `nothing` to `u0`" begin
    @variables x(t) = 1
    k = ShiftIndex()
    @mtkcompile sys = System([x(k) ~ x(k - 1) + 1], t)
    prob = @test_nowarn DiscreteProblem(sys, nothing, (0.0, 1.0))
    sol = solve(prob, FunctionMap())
    @test SciMLBase.successful_retcode(sol)
end

@testset "Shifted array variables" begin
    @variables x(t)[1:2] y(t)[1:2]
    k = ShiftIndex(t)
    eqs = [
        x(k) ~ x(k - 1) + x(k - 2),
        y[1](k) ~ y[1](k - 1) + y[1](k - 2),
        y[2](k) ~ y[2](k - 1) + y[2](k - 2)
    ]
    @mtkcompile sys = System(eqs, t)
    prob = DiscreteProblem(sys,
        [x(k - 1) => ones(2), x(k - 2) => zeros(2), y[1](k - 1) => 1.0,
            y[1](k - 2) => 0.0, y[2](k - 1) => 1.0, y[2](k - 2) => 0.0],
        (0, 4))
    @test all(isone, prob.u0)
    sol = solve(prob, FunctionMap())
    @test sol[[x..., y...], end] == 8ones(4)
end

@testset "Defaults are totermed appropriately" begin
    @parameters σ ρ β q
    @variables x(t) y(t) z(t)
    k = ShiftIndex(t)
    p = [σ => 28.0,
        ρ => 10.0,
        β => 8 / 3]

    @mtkcompile discsys = System(
        [x ~ x(k - 1) * ρ + y(k - 2), y ~ y(k - 1) * σ - z(k - 2),
            z ~ z(k - 1) * β + x(k - 2)],
        t; defaults = [x => 1.0, y => 1.0, z => 1.0, x(k - 1) => 1.0,
            y(k - 1) => 1.0, z(k - 1) => 1.0])
    discprob = DiscreteProblem(discsys, p, (0, 10))
    sol = solve(discprob, FunctionMap())
    @test SciMLBase.successful_retcode(sol)
end
