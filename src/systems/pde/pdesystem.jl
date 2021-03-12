"""
$(TYPEDEF)

A system of partial differential equations.

# Fields
$(FIELDS)

# Example

```julia
using ModelingToolkit

@parameters t x
@variables u(..)
Dxx = Differential(x)^2
Dtt = Differential(t)^2
Dt = Differential(t)

#2D PDE
C=1
eq  = Dtt(u(t,x)) ~ C^2*Dxx(u(t,x))

# Initial and boundary conditions
bcs = [u(t,0) ~ 0.,# for all t > 0
       u(t,1) ~ 0.,# for all t > 0
       u(0,x) ~ x*(1. - x), #for all 0 < x < 1
       Dt(u(0,x)) ~ 0. ] #for all  0 < x < 1]

# Space and time domains
domains = [t ∈ IntervalDomain(0.0,1.0),
           x ∈ IntervalDomain(0.0,1.0)]

pde_system = PDESystem(eq,bcs,domains,[t,x],[u])
```
"""
struct PDESystem <: ModelingToolkit.AbstractSystem
  "The equations which define the PDE"
  eqs
  "The boundary conditions"
  bcs
  "The domain for the independent variables."
  domain
  "The independent variables"
  indvars
  "The dependent variables"
  depvars
  "The parameters"
  ps
  "The default values of the parameters"
  default_p
  @add_kwonly function PDESystem(eqs, bcs, domain, indvars, depvars, ps = SciMLBase.NullParameters(), default_p = nothing)
      new(eqs, bcs, domain, indvars, depvars, ps, default_p)
  end
end

Base.getproperty(x::PDESystem, sym::Symbol) = getfield(x, sym)

Base.summary(prob::PDESystem) = string(nameof(typeof(prob)))
function Base.show(io::IO, sys::PDESystem)
  println(io,summary(sys))
  println(io,"Equations: ", sys.eqs)
  println(io,"Boundary Conditions: ", sys.bcs)
  println(io,"Domain: ", sys.domain)
  println(io,"Dependent Variables: ", sys.depvars)
  println(io,"Independent Variables: ", sys.indvars)
  println(io,"Parameters: ", sys.ps)
  println(io,"Default Parameter Values", sys.default_p)
  return nothing
end
