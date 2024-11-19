using OrdinaryDiffEq, ModelingToolkit, DataStructures, Test
using Optimization, RecursiveArrayTools, OptimizationOptimJL
using SymbolicIndexingInterface
using ModelingToolkit: t_nounits as t, D_nounits as D
using SciMLBase: parameterless_type

N = 32
const xyd_brusselator = range(0, stop = 1, length = N)
brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0
limit(a, N) = ModelingToolkit.ifelse(a == N + 1, 1, ModelingToolkit.ifelse(a == 0, N, a))
function brusselator_2d_loop(du, u, p, t)
    A, B, alpha, dx = p
    alpha = alpha / dx^2
    @inbounds for I in CartesianIndices((N, N))
        i, j = Tuple(I)
        x, y = xyd_brusselator[I[1]], xyd_brusselator[I[2]]
        ip1, im1, jp1, jm1 = limit(i + 1, N), limit(i - 1, N), limit(j + 1, N),
        limit(j - 1, N)
        du[i, j, 1] = alpha * (u[im1, j, 1] + u[ip1, j, 1] + u[i, jp1, 1] + u[i, jm1, 1] -
                       4u[i, j, 1]) +
                      B + u[i, j, 1]^2 * u[i, j, 2] - (A + 1) * u[i, j, 1] +
                      brusselator_f(x, y, t)
        du[i, j, 2] = alpha * (u[im1, j, 2] + u[ip1, j, 2] + u[i, jp1, 2] + u[i, jm1, 2] -
                       4u[i, j, 2]) +
                      A * u[i, j, 1] - u[i, j, 1]^2 * u[i, j, 2]
    end
end

# Test with tuple parameters
p = (3.4, 1.0, 10.0, step(xyd_brusselator))

function init_brusselator_2d(xyd)
    N = length(xyd)
    u = zeros(N, N, 2)
    for I in CartesianIndices((N, N))
        x = xyd[I[1]]
        y = xyd[I[2]]
        u[I, 1] = 22 * (y * (1 - y))^(3 / 2)
        u[I, 2] = 27 * (x * (1 - x))^(3 / 2)
    end
    u
end
u0 = init_brusselator_2d(xyd_brusselator)

# Test with 3-tensor inputs
prob_ode_brusselator_2d = ODEProblem(brusselator_2d_loop,
    u0, (0.0, 11.5), p)

modelingtoolkitize(prob_ode_brusselator_2d)

## Optimization

rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
p = [1.0, 100.0]

prob = OptimizationProblem(rosenbrock, x0, p)
sys = complete(modelingtoolkitize(prob)) # symbolicitize me captain!

prob = OptimizationProblem(sys, x0, p, grad = true, hess = true)
sol = solve(prob, NelderMead())
@test sol.objective < 1e-8

sol = solve(prob, BFGS())
@test sol.objective < 1e-8

sol = solve(prob, Newton())
@test sol.objective < 1e-8

prob = OptimizationProblem(ones(3); lb = [-Inf, 0.0, 1.0], ub = [Inf, 0.0, 2.0]) do u, p
    sum(abs2, u)
end

sys = complete(modelingtoolkitize(prob))
@test !ModelingToolkit.hasbounds(unknowns(sys)[1])
@test !ModelingToolkit.hasbounds(unknowns(sys)[2])
@test ModelingToolkit.hasbounds(unknowns(sys)[3])
@test ModelingToolkit.getbounds(unknowns(sys)[3]) == (1.0, 2.0)

## SIR System Regression Test

β = 0.01# infection rate
λ_R = 0.05 # inverse of transition time from  infected to recovered
λ_D = 0.83 # inverse of transition time from  infected to dead
i₀ = 0.075 # fraction of initial infected people in every age class
𝒫 = vcat([β, λ_R, λ_D]...)

# regional contact matrix and regional population

## regional contact matrix
regional_all_contact_matrix = [3.45536 0.485314 0.506389 0.123002;
                               0.597721 2.11738 0.911374 0.323385;
                               0.906231 1.35041 1.60756 0.67411;
                               0.237902 0.432631 0.726488 0.979258] # 4x4 contact matrix

## regional population stratified by age
N = [723208, 874150, 1330993, 1411928] # array of 4 elements, each of which representing the absolute amount of population in the corresponding age class.

# Initial conditions
I₀ = repeat([i₀], 4)
S₀ = N .- I₀
R₀ = [0.0 for n in 1:length(N)]
D₀ = [0.0 for n in 1:length(N)]
D_tot₀ = [0.0 for n in 1:length(N)]
ℬ = vcat([S₀, I₀, R₀, D₀, D_tot₀]...)

# Time
final_time = 20
𝒯 = (1.0, final_time);

function SIRD_ac!(du, u, p, t)
    # Parameters to be calibrated
    β, λ_R, λ_D = p

    # initialize this parameter (death probability stratified by age, taken from literature)

    δ₁, δ₂, δ₃, δ₄ = [
        0.003 / 100,
        0.004 / 100,
        (0.015 + 0.030 + 0.064 + 0.213 + 0.718) / (5 * 100),
        (2.384 + 8.466 + 12.497 + 1.117) / (4 * 100)
    ]
    δ = vcat(repeat([δ₁], 1), repeat([δ₂], 1), repeat([δ₃], 1), repeat([δ₄], 4 - 1 - 1 - 1))

    C = regional_all_contact_matrix

    # Unknown variables
    S = @view u[(4 * 0 + 1):(4 * 1)]
    I = @view u[(4 * 1 + 1):(4 * 2)]
    R = @view u[(4 * 2 + 1):(4 * 3)]
    D = @view u[(4 * 3 + 1):(4 * 4)]
    D_tot = @view u[(4 * 4 + 1):(4 * 5)]

    # Differentials
    dS = @view du[(4 * 0 + 1):(4 * 1)]
    dI = @view du[(4 * 1 + 1):(4 * 2)]
    dR = @view du[(4 * 2 + 1):(4 * 3)]
    dD = @view du[(4 * 3 + 1):(4 * 4)]
    dD_tot = @view du[(4 * 4 + 1):(4 * 5)]

    # Force of infection
    Λ = β * [sum([C[i, j] * I[j] / N[j] for j in 1:size(C)[1]]) for i in 1:size(C)[2]]

    # System of equations
    @. dS = -Λ * S
    @. dI = Λ * S - ((1 - δ) * λ_R + δ * λ_D) * I
    @. dR = λ_R * (1 - δ) * I
    @. dD = λ_D * δ * I
    @. dD_tot = dD[1] + dD[2] + dD[3] + dD[4]
end;

# create problem and check it works
problem = ODEProblem(SIRD_ac!, ℬ, 𝒯, 𝒫)
@time solution = solve(problem, Tsit5(), saveat = 1:final_time);

problem = ODEProblem(SIRD_ac!, ℬ, 𝒯, 𝒫)
sys = complete(modelingtoolkitize(problem))
fast_problem = ODEProblem(sys, ℬ, 𝒯, parameters(sys) .=> 𝒫)
@time solution = solve(fast_problem, Tsit5(), saveat = 1:final_time)

## Issue #778

r0 = [1131.340, -2282.343, 6672.423]
v0 = [-5.64305, 4.30333, 2.42879]
Δt = 86400.0 * 365
μ = 398600.4418
rv0 = ArrayPartition(r0, v0)

f = function (dy, y, μ, t)
    r = sqrt(sum(y[1, :] .^ 2))
    dy[1, :] = y[2, :]
    dy[2, :] = -μ .* y[1, :] / r^3
end

prob = ODEProblem(f, rv0, (0.0, Δt), μ)
modelingtoolkitize(prob)

# Index reduction and mass matrix handling
using LinearAlgebra
function pendulum!(du, u, p, t)
    x, dx, y, dy, T = u
    g, L = p
    du[1] = dx
    du[2] = T * x
    du[3] = dy
    du[4] = T * y - g
    du[5] = x^2 + y^2 - L^2
    return nothing
end
pendulum_fun! = ODEFunction(pendulum!, mass_matrix = Diagonal([1, 1, 1, 1, 0]))
u0 = [1.0, 0, 0, 0, 0]
p = [9.8, 1]
tspan = (0, 10.0)
pendulum_prob = ODEProblem(pendulum_fun!, u0, tspan, p)
pendulum_sys_org = complete(modelingtoolkitize(pendulum_prob))
sts = unknowns(pendulum_sys_org)
pendulum_sys = dae_index_lowering(pendulum_sys_org)
prob = ODEProblem(pendulum_sys, Pair[], tspan)
sol = solve(prob, Rodas4())
l2 = sol[sts[1]] .^ 2 + sol[sts[3]] .^ 2
@test all(l -> abs(sqrt(l) - 1) < 0.05, l2)

ff911 = (du, u, p, t) -> begin
    du[1] = u[2] + 1.0
    du[2] = u[1] - 1.0
end
prob = ODEProblem(ff911, zeros(2), (0, 1.0))
@test_nowarn modelingtoolkitize(prob)

k(x, p, t) = p * x
x0 = 1.0
p = 0.98
tspan = (0.0, 1.0)
prob = ODEProblem(k, x0, tspan, p)
sys = modelingtoolkitize(prob)

k(x, p, t) = 0.98 * x
x0 = 1.0
tspan = (0.0, 1.0)
prob = ODEProblem(k, x0, tspan)
sys = modelingtoolkitize(prob)

# https://github.com/SciML/ModelingToolkit.jl/issues/1158

function ode_prob(du, u, p::NamedTuple, t)
    du[1] = u[1] + p.α * u[2]
    du[2] = u[2] + p.β * u[1]
end
params = (α = 1, β = 1)
prob = ODEProblem(ode_prob, [1 1], (0, 1), params)
sys = modelingtoolkitize(prob)
@test nameof.(parameters(sys)) == [:α, :β]

function ode_prob(du, u, p::Tuple, t)
    α, β = p
    du[1] = u[1] + α * u[2]
    du[2] = u[2] + β * u[1]
end

params = (1, 1)
prob = ODEProblem(ode_prob, [1 1], (0, 1), params)
sys = modelingtoolkitize(prob)
@test nameof.(parameters(sys)) == [:α₁, :α₂]

function ode_prob_dict(du, u, p, t)
    du[1] = u[1] + p[:a]
    du[2] = u[2] + p[:b]
    nothing
end
params = OrderedDict(:a => 10, :b => 20)
u0 = [1, 2.0]
prob = ODEProblem(ode_prob_dict, u0, (0.0, 1.0), params)
sys = modelingtoolkitize(prob)
@test [ModelingToolkit.defaults(sys)[s] for s in unknowns(sys)] == u0
@test [ModelingToolkit.defaults(sys)[s] for s in parameters(sys)] == [10, 20]
@test ModelingToolkit.has_tspan(sys)

@parameters sig=10 rho=28.0 beta=8 / 3
@variables x(t)=100 y(t)=1.0 z(t)=1

eqs = [D(x) ~ sig * (y - x),
    D(y) ~ x * (rho - z) - y,
    D(z) ~ x * y - beta * z]

noiseeqs = [0.1 * x,
    0.1 * y,
    0.1 * z]

@named sys = SDESystem(eqs, noiseeqs, t, [x, y, z], [sig, rho, beta]; tspan = (0, 1000.0))
prob = SDEProblem(complete(sys))
sys = modelingtoolkitize(prob)
@test ModelingToolkit.has_tspan(sys)

@testset "Explicit variable names" begin
    function fn(du, u, p::NamedTuple, t)
        du[1] = u[1] + p.a * u[2]
        du[2] = u[2] + p.b * u[1]
    end
    function fn(du, u, p::AbstractDict, t)
        du[1] = u[1] + p[:a] * u[2]
        du[2] = u[2] + p[:b] * u[1]
    end
    function fn(du, u, p, t)
        du[1] = u[1] + p[1] * u[2]
        du[2] = u[2] + p[2] * u[1]
    end
    function fn(du, u, p::Real, t)
        du[1] = u[1] + p * u[2]
        du[2] = u[2] + p * u[1]
    end
    function nl_fn(u, p::NamedTuple)
        [u[1] + p.a * u[2], u[2] + p.b * u[1]]
    end
    function nl_fn(u, p::AbstractDict)
        [u[1] + p[:a] * u[2], u[2] + p[:b] * u[1]]
    end
    function nl_fn(u, p)
        [u[1] + p[1] * u[2], u[2] + p[2] * u[1]]
    end
    function nl_fn(u, p::Real)
        [u[1] + p * u[2], u[2] + p * u[1]]
    end
    params = (a = 1, b = 1)
    odeprob = ODEProblem(fn, [1 1], (0, 1), params)
    nlprob = NonlinearProblem(nl_fn, [1, 1], params)
    optprob = OptimizationProblem(nl_fn, [1, 1], params)

    @testset "$(parameterless_type(prob))" for prob in [optprob]
        sys = modelingtoolkitize(prob, u_names = [:a, :b])
        @test is_variable(sys, sys.a)
        @test is_variable(sys, sys.b)
        @test is_variable(sys, :a)
        @test is_variable(sys, :b)

        @test_throws ["unknowns", "2", "does not match", "names", "3"] modelingtoolkitize(
            prob, u_names = [:a, :b, :c])
        for (pvals, pnames) in [
            ([1, 2], [:p, :q]),
            ((1, 2), [:p, :q]),
            ([1, 2], Dict(1 => :p, 2 => :q)),
            ((1, 2), Dict(1 => :p, 2 => :q)),
            (1.0, :p),
            (1.0, [:p]),
            (1.0, Dict(1 => :p)),
            (Dict(:a => 2, :b => 4), Dict(:a => :p, :b => :q)),
            ((a = 1, b = 2), (a = :p, b = :q)),
            ((a = 1, b = 2), Dict(:a => :p, :b => :q))
        ]
            if pvals isa NamedTuple && prob isa OptimizationProblem
                continue
            end
            sys = modelingtoolkitize(
                remake(prob, p = pvals, interpret_symbolicmap = false), p_names = pnames)
            if pnames isa Symbol
                @test is_parameter(sys, pnames)
                continue
            end
            for p in values(pnames)
                @test is_parameter(sys, p)
            end
        end

        for (pvals, pnames) in [
            ([1, 2], [:p, :q, :r]),
            ((1, 2), [:p, :q, :r]),
            ([1, 2], Dict(1 => :p, 2 => :q, 3 => :r)),
            ((1, 2), Dict(1 => :p, 2 => :q, 3 => :r)),
            (1.0, [:p, :q]),
            (1.0, Dict(1 => :p, 2 => :q)),
            (Dict(:a => 2, :b => 4), Dict(:a => :p, :b => :q, :c => :r)),
            ((a = 1, b = 2), (a = :p, b = :q, c = :r)),
            ((a = 1, b = 2), Dict(:a => :p, :b => :q, :c => :r))
        ]
            newprob = remake(prob, p = pvals, interpret_symbolicmap = false)
            @test_throws [
                "parameters", "$(length(pvals))", "does not match", "$(length(pnames))"] modelingtoolkitize(
                newprob, p_names = pnames)
        end

        sc = SymbolCache([:a, :b], [:p, :q])
        sci_f = parameterless_type(prob.f)(prob.f.f, sys = sc)
        newprob = remake(prob, f = sci_f, p = [1, 2])
        sys = modelingtoolkitize(newprob)
        @test is_variable(sys, sys.a)
        @test is_variable(sys, sys.b)
        @test is_variable(sys, :a)
        @test is_variable(sys, :b)
        @test is_parameter(sys, sys.p)
        @test is_parameter(sys, sys.q)
        @test is_parameter(sys, :p)
        @test is_parameter(sys, :q)
    end

    @testset "From MTK model" begin
        @testset "ODE" begin
            @variables x(t)=1.0 y(t)=2.0
            @parameters p=3.0 q=4.0
            @mtkbuild sys = ODESystem([D(x) ~ p * y, D(y) ~ q * x], t)
            prob1 = ODEProblem(sys, [], (0.0, 5.0))
            newsys = complete(modelingtoolkitize(prob1))
            @test is_variable(newsys, newsys.x)
            @test is_variable(newsys, newsys.y)
            @test is_parameter(newsys, newsys.p)
            @test is_parameter(newsys, newsys.q)
            prob2 = ODEProblem(newsys, [], (0.0, 5.0))

            sol1 = solve(prob1, Tsit5())
            sol2 = solve(prob2, Tsit5())
            @test sol1 ≈ sol2
        end
        @testset "Nonlinear" begin
            @variables x=1.0 y=2.0
            @parameters p=3.0 q=4.0
            @mtkbuild nlsys = NonlinearSystem([0 ~ p * y^2 + x, 0 ~ x + exp(x) * q])
            prob1 = NonlinearProblem(nlsys, [])
            newsys = complete(modelingtoolkitize(prob1))
            @test is_variable(newsys, newsys.x)
            @test is_variable(newsys, newsys.y)
            @test is_parameter(newsys, newsys.p)
            @test is_parameter(newsys, newsys.q)
            prob2 = NonlinearProblem(newsys, [])

            sol1 = solve(prob1)
            sol2 = solve(prob2)
            @test sol1 ≈ sol2
        end
        @testset "Optimization" begin
            @variables begin
                x = 1.0, [bounds = (-2.0, 10.0)]
                y = 2.0, [bounds = (-1.0, 10.0)]
            end
            @parameters p=3.0 q=4.0
            loss = (p - x)^2 + q * (y - x^2)^2
            @mtkbuild optsys = OptimizationSystem(loss, [x, y], [p, q])
            prob1 = OptimizationProblem(optsys, [], grad = true, hess = true)
            newsys = complete(modelingtoolkitize(prob1))
            @test is_variable(newsys, newsys.x)
            @test is_variable(newsys, newsys.y)
            @test is_parameter(newsys, newsys.p)
            @test is_parameter(newsys, newsys.q)
            prob2 = OptimizationProblem(newsys, [], grad = true, hess = true)

            sol1 = solve(prob1, GradientDescent())
            sol2 = solve(prob2, GradientDescent())

            @test sol1 ≈ sol2
        end
    end
end

## NonlinearLeastSquaresProblem

function nlls!(du, u, p)
    du[1] = 2u[1] - 2
    du[2] = u[1] - 4u[2]
    du[3] = 0
end
u0 = [0.0, 0.0]
prob = NonlinearLeastSquaresProblem(
    NonlinearFunction(nlls!, resid_prototype = zeros(3)), u0)
sys = modelingtoolkitize(prob)
@test length(equations(sys)) == 3
@test length(equations(structural_simplify(sys; fully_determined = false))) == 0

@testset "`modelingtoolkitize(::SDEProblem)` sets defaults" begin
    function sdeg!(du, u, p, t)
        du[1] = 0.3 * u[1]
        du[2] = 0.3 * u[2]
        du[3] = 0.3 * u[3]
    end
    function sdef!(du, u, p, t)
        x, y, z = u
        sigma, rho, beta = p
        du[1] = sigma * (y - x)
        du[2] = x * (rho - z) - y
        du[3] = x * y - beta * z
    end
    u0 = [1.0, 0.0, 0.0]
    tspan = (0.0, 100.0)
    p = [10.0, 28.0, 2.66]
    sprob = SDEProblem(sdef!, sdeg!, u0, tspan, p)
    sys = complete(modelingtoolkitize(sprob))
    @test length(ModelingToolkit.defaults(sys)) == 6
    sprob2 = SDEProblem(sys, [], tspan)

    truevals = similar(u0)
    sprob.f(truevals, u0, p, tspan[1])
    mtkvals = similar(u0)
    sprob2.f(mtkvals, sprob2.u0, sprob2.p, tspan[1])
    @test mtkvals ≈ truevals
end
