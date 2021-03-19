using ModelingToolkit, LinearAlgebra, OrdinaryDiffEq, Test
MT = ModelingToolkit

# Repressilator model
@parameters t α₀ α K n δ β μ 
@variables m(t) P(t) R(t)
rxs = [Reaction(α₀, nothing, [m]),
       Reaction(α / (1 + (R/K)^n), nothing, [m]),
       Reaction(δ, [m], nothing),
       Reaction(β, [m], [m,P]),
       Reaction(μ, [P], nothing)    
]

specs = [m,P,R]
pars  = [α₀,α,K,n,δ,β,μ]
@named rs = ReactionSystem(rxs, t, specs, pars)

# using ODESystem components
@named os₁ = convert(ODESystem, rs)
@named os₂ = convert(ODESystem, rs)
@named os₃ = convert(ODESystem, rs)
connections = [os₁.R ~ os₃.P, 
               os₂.R ~ os₁.P,
               os₃.R ~ os₂.P]
@named connected = ODESystem(connections, t, [], [], systems=[os₁,os₂,os₃])
oderepressilator = structural_simplify(connected)

pvals = [os₁.α₀ => 5e-4, 
         os₁.α => .5, 
         os₁.K => 40.0, 
         os₁.n => 2, 
         os₁.δ => (log(2)/120), 
         os₁.β => (20*log(2)/120),
         os₁.μ => (log(2)/600),
         os₂.α₀ => 5e-4, 
         os₂.α => .5, 
         os₂.K => 40.0, 
         os₂.n => 2, 
         os₂.δ => (log(2)/120), 
         os₂.β => (20*log(2)/120),
         os₂.μ => (log(2)/600),
         os₃.α₀ => 5e-4, 
         os₃.α => .5, 
         os₃.K => 40.0, 
         os₃.n => 2, 
         os₃.δ => (log(2)/120), 
         os₃.β => (20*log(2)/120),
         os₃.μ => (log(2)/600)]
u₀    = [os₁.m => 0.0, os₁.P => 20.0, os₂.m => 0.0, os₂.P => 0.0, os₃.m => 0.0, os₃.P => 0.0]
tspan = (0.0, 100000.0)
oprob = ODEProblem(oderepressilator, u₀, tspan, pvals)
sol = solve(oprob, Tsit5())

# hardcoded network
function repress!(f, y, p, t)
    α = p.α; α₀ = p.α₀; β = p.β; δ = p.δ; μ = p.μ; K = p.K; n = p.n
    f[1] = α / (1 + (y[6] / K)^n) - δ * y[1] + α₀
    f[2] = α / (1 + (y[4] / K)^n) - δ * y[2] + α₀
    f[3] = α / (1 + (y[5] / K)^n) - δ * y[3] + α₀
    f[4] = β * y[1] - μ * y[4]
    f[5] = β * y[2] - μ * y[5]    
    f[6] = β * y[3] - μ * y[6]    
    nothing
end
ps = (α₀=5e-4, α=.5, K=40.0, n=2, δ=(log(2)/120), β=(20*log(2)/120), μ=(log(2)/600))
u0 = [0.0,0.0,0.0,20.0,0.0,0.0]
oprob2 = ODEProblem(repress!, u0, tspan, ps)
sol2 = solve(oprob2, Tsit5())
tvs = 0:1:tspan[end]

indexof(sym,syms) = findfirst(isequal(sym),syms)
i = indexof(os₁.P, states(oderepressilator))
@test all(isapprox(u[1],u[2],atol=1e-4) for u in zip(sol(tvs, idxs=2), sol2(tvs, idxs=4)))

# using ReactionSystem components

# @named rs₁ = ReactionSystem(rxs, t, specs, pars)
# @named rs₂ = ReactionSystem(rxs, t, specs, pars)
# @named rs₃ = ReactionSystem(rxs, t, specs, pars)
# connections = [rs₁.R ~ rs₃.P, 
#                rs₂.R ~ rs₁.P,
#                rs₃.R ~ rs₂.P]
# @named csys = ODESystem(connections, t, [], [])
# @named repressilator = ReactionSystem(t; systems=[csys,rs₁,rs₂,rs₃])
# @named oderepressilator2 = convert(ODESystem, repressilator)
# sys2 = structural_simplify(oderepressilator2)  # FAILS currently