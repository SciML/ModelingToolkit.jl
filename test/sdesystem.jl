using ModelingToolkit, StaticArrays, LinearAlgebra
using StochasticDiffEq, SparseArrays
using Random,Test

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

noiseeqs = [0.1*x,
            0.1*y,
            0.1*z]

de = SDESystem(eqs,noiseeqs,t,[x,y,z],[σ,ρ,β])
f = eval(generate_diffusion_function(de)[1])
@test f(ones(3),rand(3),nothing) == 0.1ones(3)

f = SDEFunction(de)
prob = SDEProblem(SDEFunction(de),f.g,[1.0,0.0,0.0],(0.0,100.0),(10.0,26.0,2.33))
sol = solve(prob,SRIW1(),seed=1)

probexpr = SDEProblem(SDEFunction(de),f.g,[1.0,0.0,0.0],(0.0,100.0),(10.0,26.0,2.33))
solexpr = solve(eval(probexpr),SRIW1(),seed=1)

@test all(x->x==0,Array(sol-solexpr))

# Test no error
@test_nowarn SDEProblem(de,nothing,(0, 10.0))

noiseeqs_nd = [0.01*x 0.01*x*y 0.02*x*z
               σ      0.01*y   0.02*x*z
               ρ      β        0.01*z  ]
de = SDESystem(eqs,noiseeqs_nd,t,[x,y,z],[σ,ρ,β])
f = eval(generate_diffusion_function(de)[1])
@test f([1,2,3.0],[0.1,0.2,0.3],nothing) == [0.01*1   0.01*1*2   0.02*1*3
                                             0.1      0.01*2     0.02*1*3
                                             0.2      0.3        0.01*3  ]

f = eval(generate_diffusion_function(de)[2])
du = ones(3,3)
f(du,[1,2,3.0],[0.1,0.2,0.3],nothing)
@test du == [0.01*1   0.01*1*2   0.02*1*3
             0.1      0.01*2     0.02*1*3
             0.2      0.3        0.01*3  ]

f = SDEFunction(de)
prob = SDEProblem(SDEFunction(de),f.g,[1.0,0.0,0.0],(0.0,100.0),(10.0,26.0,2.33),
                  noise_rate_prototype = zeros(3,3))
sol = solve(prob,EM(),dt=0.001)

u0map = [
    x => 1.0,
    y => 0.0,
    z => 0.0
]

parammap = [
    σ => 10.0,
    β => 26.0,
    ρ => 2.33
]

prob = SDEProblem(de,u0map,(0.0,100.0),parammap)
@test size(prob.noise_rate_prototype) == (3,3)
@test prob.noise_rate_prototype isa Matrix
sol = solve(prob,EM(),dt=0.001)

prob = SDEProblem(de,u0map,(0.0,100.0),parammap,sparsenoise=true)
@test size(prob.noise_rate_prototype) == (3,3)
@test prob.noise_rate_prototype isa SparseMatrixCSC
sol = solve(prob,EM(),dt=0.001)

# Test eval_expression=false
function test_SDEFunction_no_eval()
    # Need to test within a function scope to trigger world age issues
    f = SDEFunction(de, eval_expression=false)
    @test f([1.0,0.0,0.0], (10.0,26.0,2.33), (0.0,100.0)) ≈ [-10.0, 26.0, 0.0]
end
test_SDEFunction_no_eval()


# modelingtoolkitize and Ito <-> Stratonovich sense
seed = 10
Random.seed!(seed)


# simple 2D diagonal noise
u0 = rand(2)
t =  randn()
trange = (0.0,100.0)
p = [1.01,0.87]
f1!(du,u,p,t) = (du .= p[1]*u)
σ1!(du,u,p,t) = (du .= p[2]*u)

prob = SDEProblem(f1!,σ1!,u0,trange,p)
# no correction
sys = modelingtoolkitize(prob)
fdrift = eval(generate_function(sys)[1])
fdif = eval(generate_diffusion_function(sys)[1])
@test fdrift(u0,p,t) == p[1]*u0
@test fdif(u0,p,t) == p[2]*u0
fdrift! = eval(generate_function(sys)[2])
fdif! = eval(generate_diffusion_function(sys)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == p[1]*u0
fdif!(du,u0,p,t)
@test du == p[2]*u0

# Ito -> Strat
sys2 = stochastic_integral_transform(sys,-1//2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0,p,t) == p[1]*u0 - 1//2*p[2]^2*u0
@test fdif(u0,p,t) == p[2]*u0
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == p[1]*u0 - 1//2*p[2]^2*u0
fdif!(du,u0,p,t)
@test du == p[2]*u0

# Strat -> Ito
sys2 = stochastic_integral_transform(sys,1//2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0,p,t) == p[1]*u0 + 1//2*p[2]^2*u0
@test fdif(u0,p,t) == p[2]*u0
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == p[1]*u0 + 1//2*p[2]^2*u0
fdif!(du,u0,p,t)
@test du == p[2]*u0

# somewhat complicated 1D without explicit parameters but with explicit time-dependence
f2!(du,u,p,t) = (du[1] = sin(t) + cos(u[1]))
σ2!(du,u,p,t) = (du[1] = pi + atan(u[1]))

u0 = rand(1)
prob = SDEProblem(f2!,σ2!,u0,trange)
# no correction
sys = modelingtoolkitize(prob)
fdrift = eval(generate_function(sys)[1])
fdif = eval(generate_diffusion_function(sys)[1])
@test fdrift(u0,p,t) ==  @. sin(t) + cos(u0)
@test fdif(u0,p,t) ==  pi .+ atan.(u0)
fdrift! = eval(generate_function(sys)[2])
fdif! = eval(generate_diffusion_function(sys)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == @. sin(t) + cos(u0)
fdif!(du,u0,p,t)
@test du == pi .+ atan.(u0)

# Ito -> Strat
sys2 = stochastic_integral_transform(sys,-1//2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0,p,t) ==  @. sin(t) + cos(u0) - 1//2*1/(1 + u0^2)*(pi + atan(u0))
@test fdif(u0,p,t) == pi .+ atan.(u0)
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du ==  @. sin(t) + cos(u0) - 1//2*1/(1 + u0^2)*(pi + atan(u0))
fdif!(du,u0,p,t)
@test du == pi .+ atan.(u0)

# Strat -> Ito
sys2 = stochastic_integral_transform(sys,1//2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0,p,t) ≈  @. sin(t) + cos(u0) + 1//2*1/(1 + u0^2)*(pi + atan(u0))
@test fdif(u0,p,t) == pi .+ atan.(u0)
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du ≈  @. sin(t) + cos(u0) + 1//2*1/(1 + u0^2)*(pi + atan(u0))
fdif!(du,u0,p,t)
@test du == pi .+ atan.(u0)


# 2D diagonal noise with mixing terms (no parameters, no time-dependence)
u0 = rand(2)
t =  randn()
function f3!(du,u,p,t)
  du[1] = u[1]/2
  du[2] = u[2]/2
  return nothing
end
function σ3!(du,u,p,t)
  du[1] = u[2]
  du[2] = u[1]
  return nothing
end

prob = SDEProblem(f3!,σ3!,u0,trange,p)
# no correction
sys = modelingtoolkitize(prob)
fdrift = eval(generate_function(sys)[1])
fdif = eval(generate_diffusion_function(sys)[1])
@test fdrift(u0,p,t) == u0/2
@test fdif(u0,p,t) == reverse(u0)
fdrift! = eval(generate_function(sys)[2])
fdif! = eval(generate_diffusion_function(sys)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == u0/2
fdif!(du,u0,p,t)
@test du ==  reverse(u0)

# Ito -> Strat
sys2 = stochastic_integral_transform(sys,-1//2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0,p,t) == u0*0
@test fdif(u0,p,t) == reverse(u0)
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == u0*0
fdif!(du,u0,p,t)
@test du == reverse(u0)

# Strat -> Ito
sys2 = stochastic_integral_transform(sys,1//2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0,p,t) == u0
@test fdif(u0,p,t) == reverse(u0)
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == u0
fdif!(du,u0,p,t)
@test du == reverse(u0)



# simple 2D diagonal noise oop
u0 = rand(2)
t =  randn()
p = [1.01,0.87]
f1(u,p,t) = p[1]*u
σ1(u,p,t) = p[2]*u

prob = SDEProblem(f1,σ1,u0,trange,p)
# no correction
sys = modelingtoolkitize(prob)
fdrift = eval(generate_function(sys)[1])
fdif = eval(generate_diffusion_function(sys)[1])
@test fdrift(u0,p,t) == p[1]*u0
@test fdif(u0,p,t) == p[2]*u0
fdrift! = eval(generate_function(sys)[2])
fdif! = eval(generate_diffusion_function(sys)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == p[1]*u0
fdif!(du,u0,p,t)
@test du == p[2]*u0

# Ito -> Strat
sys2 = stochastic_integral_transform(sys,-1//2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0,p,t) == p[1]*u0 - 1//2*p[2]^2*u0
@test fdif(u0,p,t) == p[2]*u0
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == p[1]*u0 - 1//2*p[2]^2*u0
fdif!(du,u0,p,t)
@test du == p[2]*u0

# Strat -> Ito
sys2 = stochastic_integral_transform(sys,1//2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0,p,t) == p[1]*u0 + 1//2*p[2]^2*u0
@test fdif(u0,p,t) == p[2]*u0
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == p[1]*u0 + 1//2*p[2]^2*u0
fdif!(du,u0,p,t)
@test du == p[2]*u0


# non-diagonal noise
u0 = rand(2)
t =  randn()
p = [1.01,0.3,0.6,1.2,0.2]
f4!(du,u,p,t) = du .= p[1]*u
function g4!(du,u,p,t)
  du[1,1] = p[2]*u[1]
  du[1,2] = p[3]*u[1]
  du[2,1] = p[4]*u[1]
  du[2,2] = p[5]*u[2]
  return nothing
end

prob = SDEProblem(f4!,g4!,u0,trange,noise_rate_prototype=zeros(2,2),p)
# no correction
sys = modelingtoolkitize(prob)
fdrift = eval(generate_function(sys)[1])
fdif = eval(generate_diffusion_function(sys)[1])
@test fdrift(u0,p,t) == p[1]*u0
@test fdif(u0,p,t) == [p[2]*u0[1]   p[3]*u0[1]
                      p[4]*u0[1]     p[5]*u0[2] ]
fdrift! = eval(generate_function(sys)[2])
fdif! = eval(generate_diffusion_function(sys)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == p[1]*u0
du = similar(u0, size(prob.noise_rate_prototype))
fdif!(du,u0,p,t)
@test du == [p[2]*u0[1]   p[3]*u0[1]
             p[4]*u0[1]   p[5]*u0[2]  ]

# Ito -> Strat
sys2 = stochastic_integral_transform(sys,-1//2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0,p,t) == [p[1]*u0[1] - 1//2*(p[2]^2*u0[1]+p[3]^2*u0[1]), p[1]*u0[2] - 1//2*(p[2]*p[4]*u0[1]+p[5]^2*u0[2])]
@test fdif(u0,p,t) == [p[2]*u0[1]   p[3]*u0[1]
                      p[4]*u0[1]     p[5]*u0[2] ]
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == [p[1]*u0[1] - 1//2*(p[2]^2*u0[1]+p[3]^2*u0[1]), p[1]*u0[2] - 1//2*(p[2]*p[4]*u0[1]+p[5]^2*u0[2])]
du = similar(u0, size(prob.noise_rate_prototype))
fdif!(du,u0,p,t)
@test du == [p[2]*u0[1]   p[3]*u0[1]
            p[4]*u0[1]     p[5]*u0[2] ]

# Strat -> Ito
sys2 = stochastic_integral_transform(sys,1//2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0,p,t) == [p[1]*u0[1] + 1//2*(p[2]^2*u0[1]+p[3]^2*u0[1]), p[1]*u0[2] + 1//2*(p[2]*p[4]*u0[1]+p[5]^2*u0[2])]
@test fdif(u0,p,t) == [p[2]*u0[1]   p[3]*u0[1]
                      p[4]*u0[1]     p[5]*u0[2] ]
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == [p[1]*u0[1] + 1//2*(p[2]^2*u0[1]+p[3]^2*u0[1]), p[1]*u0[2] + 1//2*(p[2]*p[4]*u0[1]+p[5]^2*u0[2])]
du = similar(u0, size(prob.noise_rate_prototype))
fdif!(du,u0,p,t)
@test du == [p[2]*u0[1]   p[3]*u0[1]
            p[4]*u0[1]     p[5]*u0[2] ]


# non-diagonal noise: Torus -- Strat and Ito are identical
u0 = rand(2)
t =  randn()
p =  rand(1)
f5!(du,u,p,t) = du .= false
function g5!(du,u,p,t)
  du[1,1] = cos(p[1])*sin(u[1])
  du[1,2] = cos(p[1])*cos(u[1])
  du[1,3] = -sin(p[1])*sin(u[2])
  du[1,4] = -sin(p[1])*cos(u[2])
  du[2,1] = sin(p[1])*sin(u[1])
  du[2,2] = sin(p[1])*cos(u[1])
  du[2,3] = cos(p[1])*sin(u[2])
  du[2,4] = cos(p[1])*cos(u[2])
  return nothing
end

prob = SDEProblem(f5!,g5!,u0,trange,noise_rate_prototype=zeros(2,4),p)
# no correction
sys = modelingtoolkitize(prob)
fdrift = eval(generate_function(sys)[1])
fdif = eval(generate_diffusion_function(sys)[1])
@test fdrift(u0,p,t) == 0*u0
@test fdif(u0,p,t) == [ cos(p[1])*sin(u0[1])   cos(p[1])*cos(u0[1])   -sin(p[1])*sin(u0[2])   -sin(p[1])*cos(u0[2])
                        sin(p[1])*sin(u0[1])   sin(p[1])*cos(u0[1])    cos(p[1])*sin(u0[2])    cos(p[1])*cos(u0[2])]
fdrift! = eval(generate_function(sys)[2])
fdif! = eval(generate_diffusion_function(sys)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du == 0*u0
du = similar(u0, size(prob.noise_rate_prototype))
fdif!(du,u0,p,t)
@test du == [ cos(p[1])*sin(u0[1])   cos(p[1])*cos(u0[1])   -sin(p[1])*sin(u0[2])   -sin(p[1])*cos(u0[2])
                        sin(p[1])*sin(u0[1])   sin(p[1])*cos(u0[1])    cos(p[1])*sin(u0[2])    cos(p[1])*cos(u0[2])]

# Ito -> Strat
sys2 = stochastic_integral_transform(sys,-1//2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0,p,t) == 0*u0
@test fdif(u0,p,t) == [ cos(p[1])*sin(u0[1])   cos(p[1])*cos(u0[1])   -sin(p[1])*sin(u0[2])   -sin(p[1])*cos(u0[2])
                        sin(p[1])*sin(u0[1])   sin(p[1])*cos(u0[1])    cos(p[1])*sin(u0[2])    cos(p[1])*cos(u0[2])]
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du ==  0*u0
du = similar(u0, size(prob.noise_rate_prototype))
fdif!(du,u0,p,t)
@test du == [ cos(p[1])*sin(u0[1])   cos(p[1])*cos(u0[1])   -sin(p[1])*sin(u0[2])   -sin(p[1])*cos(u0[2])
              sin(p[1])*sin(u0[1])   sin(p[1])*cos(u0[1])    cos(p[1])*sin(u0[2])    cos(p[1])*cos(u0[2])]

# Strat -> Ito
sys2 = stochastic_integral_transform(sys,1//2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0,p,t) ==  0*u0
@test fdif(u0,p,t) ==  [ cos(p[1])*sin(u0[1])   cos(p[1])*cos(u0[1])   -sin(p[1])*sin(u0[2])   -sin(p[1])*cos(u0[2])
                         sin(p[1])*sin(u0[1])   sin(p[1])*cos(u0[1])    cos(p[1])*sin(u0[2])    cos(p[1])*cos(u0[2])]
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du,u0,p,t)
@test  du ==  0*u0
du = similar(u0, size(prob.noise_rate_prototype))
fdif!(du,u0,p,t)
@test du == [ cos(p[1])*sin(u0[1])   cos(p[1])*cos(u0[1])   -sin(p[1])*sin(u0[2])   -sin(p[1])*cos(u0[2])
              sin(p[1])*sin(u0[1])   sin(p[1])*cos(u0[1])    cos(p[1])*sin(u0[2])    cos(p[1])*cos(u0[2])]
