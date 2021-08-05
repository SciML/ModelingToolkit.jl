using OrdinaryDiffEq, ModelingToolkit, Test
using GalacticOptim, Optim, RecursiveArrayTools

N = 32
const xyd_brusselator = range(0,stop=1,length=N)
brusselator_f(x, y, t) = (((x-0.3)^2 + (y-0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.
limit(a, N) = ModelingToolkit.ifelse(a == N+1, 1, ModelingToolkit.ifelse(a == 0, N, a))
function brusselator_2d_loop(du, u, p, t)
  A, B, alpha, dx = p
  alpha = alpha/dx^2
  @inbounds for I in CartesianIndices((N, N))
    i, j = Tuple(I)
    x, y = xyd_brusselator[I[1]], xyd_brusselator[I[2]]
    ip1, im1, jp1, jm1 = limit(i+1, N), limit(i-1, N), limit(j+1, N), limit(j-1, N)
    du[i,j,1] = alpha*(u[im1,j,1] + u[ip1,j,1] + u[i,jp1,1] + u[i,jm1,1] - 4u[i,j,1]) +
                B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1] + brusselator_f(x, y, t)
    du[i,j,2] = alpha*(u[im1,j,2] + u[ip1,j,2] + u[i,jp1,2] + u[i,jm1,2] - 4u[i,j,2]) +
                A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]
    end
end

# Test with tuple parameters
p = (3.4, 1., 10., step(xyd_brusselator))

function init_brusselator_2d(xyd)
  N = length(xyd)
  u = zeros(N, N, 2)
  for I in CartesianIndices((N, N))
    x = xyd[I[1]]
    y = xyd[I[2]]
    u[I,1] = 22*(y*(1-y))^(3/2)
    u[I,2] = 27*(x*(1-x))^(3/2)
  end
  u
end
u0 = init_brusselator_2d(xyd_brusselator)

# Test with 3-tensor inputs
prob_ode_brusselator_2d = ODEProblem(brusselator_2d_loop,
                                     u0,(0.,11.5),p)

modelingtoolkitize(prob_ode_brusselator_2d)

## Optimization

rosenbrock(x,p) =  (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
p  = [1.0,100.0]

prob = OptimizationProblem(rosenbrock,x0,p)
sys = modelingtoolkitize(prob) # symbolicitize me captain!

prob = OptimizationProblem(sys,x0,p,grad=true,hess=true)
sol = solve(prob,NelderMead())
@test sol.minimum < 1e-8

sol = solve(prob,BFGS())
@test sol.minimum < 1e-8

sol = solve(prob,Newton())
@test sol.minimum < 1e-8

## SIR System Regression Test

β = 0.01# infection rate
λ_R = 0.05 # inverse of transition time from  infected to recovered
λ_D = 0.83 # inverse of transition time from  infected to dead
i₀ = 0.075 # fraction of initial infected people in every age class
𝒫 = vcat([β, λ_R, λ_D]...)

# regional contact matrix and regional population

## regional contact matrix
regional_all_contact_matrix = [3.45536   0.485314  0.506389  0.123002 ; 0.597721  2.11738   0.911374  0.323385 ; 0.906231  1.35041   1.60756   0.67411 ; 0.237902  0.432631  0.726488  0.979258] # 4x4 contact matrix

## regional population stratified by age
N = [723208 , 874150, 1330993, 1411928] # array of 4 elements, each of which representing the absolute amount of population in the corresponding age class.


# Initial conditions
I₀ = repeat([i₀],4)
S₀ = N.-I₀
R₀ = [0.0 for n in 1:length(N)]
D₀ = [0.0 for n in 1:length(N)]
D_tot₀ = [0.0 for n in 1:length(N)]
ℬ = vcat([S₀, I₀, R₀, D₀, D_tot₀]...)

# Time
final_time = 20
𝒯 = (1.0,final_time);




function SIRD_ac!(du,u,p,t)
    # Parameters to be calibrated
    β, λ_R, λ_D = p

    # initialize this parameter (death probability stratified by age, taken from literature)

    δ₁, δ₂, δ₃, δ₄ = [0.003/100, 0.004/100, (0.015+0.030+0.064+0.213+0.718)/(5*100), (2.384+8.466+12.497+1.117)/(4*100)]
    δ = vcat(repeat([δ₁],1),repeat([δ₂],1),repeat([δ₃],1),repeat([δ₄],4-1-1-1))


    C = regional_all_contact_matrix


    # State variables
    S = @view u[4*0+1:4*1]
    I = @view u[4*1+1:4*2]
    R = @view u[4*2+1:4*3]
    D = @view u[4*3+1:4*4]
    D_tot = @view u[4*4+1:4*5]

    # Differentials
    dS = @view du[4*0+1:4*1]
    dI = @view du[4*1+1:4*2]
    dR = @view du[4*2+1:4*3]
    dD = @view du[4*3+1:4*4]
    dD_tot = @view du[4*4+1:4*5]

    # Force of infection
    Λ = β*[sum([C[i,j]*I[j]/N[j] for j in 1:size(C)[1]]) for i in 1:size(C)[2]]

    # System of equations
    @. dS = -Λ*S
    @. dI = Λ*S - ((1-δ)*λ_R + δ*λ_D)*I
    @. dR = λ_R*(1-δ)*I
    @. dD = λ_D*δ*I
    @. dD_tot = dD[1]+dD[2]+dD[3]+dD[4]


end;


# create problem and check it works
problem = ODEProblem(SIRD_ac!, ℬ, 𝒯, 𝒫)
@time solution = solve(problem, Tsit5(), saveat = 1:final_time);

problem = ODEProblem(SIRD_ac!, ℬ, 𝒯, 𝒫)
sys = modelingtoolkitize(problem)
fast_problem = ODEProblem(sys,ℬ, 𝒯, 𝒫 )
@time solution = solve(fast_problem, Tsit5(), saveat = 1:final_time)

## Issue #778

r0 = [1131.340, -2282.343, 6672.423]
v0 = [-5.64305, 4.30333, 2.42879]
Δt = 86400.0*365
μ = 398600.4418
rv0 = ArrayPartition(r0,v0)

function f(dy, y, μ, t)
    r = sqrt(sum(y[1,:].^2))
    dy[1,:] = y[2,:]
    dy[2,:] = -μ .* y[1,:] / r^3
end

prob = ODEProblem(f, rv0, (0.0, Δt), μ)
sol = solve(prob, Vern8())

modelingtoolkitize(prob)

# Index reduction and mass matrix handling
using LinearAlgebra
function pendulum!(du, u, p, t)
    x, dx, y, dy, T = u
    g, L = p
    du[1] = dx
    du[2] = T*x
    du[3] = dy
    du[4] = T*y - g
    du[5] = x^2 + y^2 - L^2
    return nothing
end
pendulum_fun! = ODEFunction(pendulum!, mass_matrix=Diagonal([1,1,1,1,0]))
u0 = [1.0, 0, 0, 0, 0]
p = [9.8, 1]
tspan = (0, 10.0)
pendulum_prob = ODEProblem(pendulum_fun!, u0, tspan, p)
pendulum_sys_org = modelingtoolkitize(pendulum_prob)
sts = states(pendulum_sys_org)
pendulum_sys = dae_index_lowering(pendulum_sys_org)
prob = ODEProblem(pendulum_sys, Pair[], tspan)
sol = solve(prob, Rodas4())
l2 = sol[sts[1]].^2 + sol[sts[3]].^2
@test all(l->abs(sqrt(l) - 1) < 0.05, l2)

ff911 = (du,u,p,t) -> begin
    du[1] = u[2] + 1.0
    du[2] = u[1] - 1.0
end
prob = ODEProblem(ff911, zeros(2), (0, 1.0))
@test_nowarn modelingtoolkitize(prob)

k(x,p,t) = p*x
x0 = 1.0
p = 0.98
tspan = (0.0,1.0)
prob = ODEProblem(k,x0,tspan,p)
sys = modelingtoolkitize(prob)

k(x,p,t) = 0.98*x
x0 = 1.0
tspan = (0.0,1.0)
prob = ODEProblem(k,x0,tspan)
sys = modelingtoolkitize(prob)


## https://github.com/SciML/ModelingToolkit.jl/issues/1054
using LabelledArrays
using ModelingToolkit

# ODE model: simple SIR model with seasonally forced contact rate
function SIR!(du,u,p,t)

    # states
    (S, I, R) = u[1:3]
    N = S + I + R

    # params
    β = p.β
    η = p.η
    φ = p.φ
    ω = 1.0/p.ω
    μ = p.μ
    σ = p.σ

    # FOI
    βeff = β * (1.0+η*cos(2.0*π*(t-φ)/365.0))
    λ = βeff*I/N

    # change in states
    du[1] = (μ*N - λ*S - μ*S + ω*R)
    du[2] = (λ*S - σ*I - μ*I)
    du[3] = (σ*I - μ*R - ω*R)
    du[4] = (σ*I) # cumulative incidence

end

# Solver settings
tmin = 0.0
tmax = 10.0*365.0
tspan = (tmin, tmax)

# Initiate ODE problem
theta_fix =  [1.0/(80*365)]
theta_est =  [0.28, 0.07, 1.0/365.0, 1.0 ,1.0/5.0]
p = @LArray [theta_est; theta_fix] (:β, :η, :ω, :φ, :σ, :μ)
u0 = @LArray [9998.0,1.0,1.0,1.0] (:S,:I,:R,:C)

# Initiate ODE problem
problem = ODEProblem(SIR!,u0,tspan,p)
sys = modelingtoolkitize(problem)

@parameters t
@test all(isequal.(parameters(sys),getproperty.(@variables(β, η, ω, φ, σ, μ),:val)))
@test all(isequal.(Symbol.(states(sys)),Symbol.(@variables(S(t),I(t),R(t),C(t)))))

# https://github.com/SciML/ModelingToolkit.jl/issues/1158

function ode_prob(du, u, p::NamedTuple, t)
    du[1] = u[1]+p.α*u[2]
    du[2] = u[2]+p.β*u[1]
end
params = (α = 1, β = 1)
prob = ODEProblem(ode_prob, [1 1], (0, 1), params)
sys = modelingtoolkitize(prob)
@test nameof.(parameters(sys)) == [:α,:β]

function ode_prob(du, u, p::Tuple, t)
    α, β = p
    du[1] = u[1]+α*u[2]
    du[2] = u[2]+β*u[1]
end

params = (1, 1)
prob = ODEProblem(ode_prob, [1 1], (0, 1), params)
sys = modelingtoolkitize(prob)
@test nameof.(parameters(sys)) == [:α₁,:α₂]
