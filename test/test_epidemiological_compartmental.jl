# Example: Compartmental models in epidemiology
#=
- https://github.com/epirecipes/sir-julia/blob/master/markdown/function_map/function_map.md
- https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#Deterministic_versus_stochastic_epidemic_models
=#
using ModelingToolkit

@inline function rate_to_proportion(r,t)
    1-exp(-r*t)
end;

# Independent and dependent variables and parameters
@parameters t c nsteps δt β γ
@variables S(t) I(t) R(t) next_S(t) next_I(t) next_R(t)

infection = rate_to_proportion(β*c*I/nsteps,δt)*S
recovery = rate_to_proportion(γ,δt)*I

# Equations
eqs = [next_S ~ S-infection,
       next_I ~ I+infection-recovery,
       next_R ~ R+recovery]

# System
sys = DiscreteSystem(eqs,t,[S,I,R],[c,nsteps,δt,β,γ])

# Problem
u0 = [S => 990.0, I => 10.0, R => 0.0]
p = [β => 0.05, c => 10.0, γ => 0.25, δt => 0.1, nsteps => 400]
tspan = (0.0,ModelingToolkit.value(substitute(nsteps,p))) # value function (from Symbolics) is used to convert a Num to Float64
prob_map = DiscreteProblem(sys,u0,tspan,p)

# Solution
using OrdinaryDiffEq
sol_map = solve(prob_map,solver=FunctionMap);
