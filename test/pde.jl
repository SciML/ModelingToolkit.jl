using ModelingToolkit, DiffEqOperators, DiffEqBase, LinearAlgebra

# Define some variables
@parameters t x
@variables u(..)
@derivatives Dt'~t
@derivatives Dxx''~x
eq  = Dt(u(t,x)) ~ Dxx(u(t,x))
bcs = [u(0,x) ~ - x * (x-1) * sin(x),
           u(t,0) ~ 0, u(t,1) ~ 0]

domains = [t ∈ IntervalDomain(0.0,1.0),
           x ∈ IntervalDomain(0.0,1.0)]

pdesys = PDESystem(eq,bcs,domains,[t,x],[u])
discretization = MOLFiniteDifference(0.1)
prob = discretize(pdesys,discretization) # This gives an ODEProblem since it's time-dependent

using OrdinaryDiffEq
sol = solve(prob,Tsit5(),saveat=0.1)
