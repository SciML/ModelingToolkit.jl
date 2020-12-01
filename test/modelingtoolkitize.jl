using OrdinaryDiffEq, ModelingToolkit, Test
using GalacticOptim, Optim

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
