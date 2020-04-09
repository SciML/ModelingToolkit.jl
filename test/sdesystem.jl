using ModelingToolkit, StaticArrays, LinearAlgebra
using DiffEqBase, StochasticDiffEq
using Test

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

noiseeqs = [0.01*x,
            0.01*y,
            0.01*z]

de = SDESystem(eqs,noiseeqs,t,[x,y,z],[σ,ρ,β])
f = eval(generate_diffusion_function(de)[1])
@test f(ones(3),rand(3),nothing) == 0.01ones(3)

f = SDEFunction(de)
prob = SDEProblem(SDEFunction(de),f.g,[1.0,0.0,0.0],(0.0,100.0),(10.0,26.0,2.33))
sol = solve(prob,SOSRI())

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
sol = solve(prob,EM(),dt=0.1)
