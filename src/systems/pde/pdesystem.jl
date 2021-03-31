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
    """
    defaults: The default values to use when initial conditions and/or
    parameters are not supplied in `ODEProblem`.
    """
    defaults::Dict
    @add_kwonly function PDESystem(eqs, bcs, domain, indvars, depvars, ps = SciMLBase.NullParameters(), defaults = Dict())
        new(eqs, bcs, domain, indvars, depvars, ps, defaults)
    end
end

Base.getproperty(x::PDESystem, sym::Symbol) = getfield(x, sym)

Base.summary(prob::PDESystem) = string(nameof(typeof(prob)))
function Base.show(io::IO, ::MIME"text/plain", sys::PDESystem)
    println(io,summary(sys))
    println(io,"Equations: ", get_eqs(sys))
    println(io,"Boundary Conditions: ", get_bcs(sys))
    println(io,"Domain: ", get_domain(sys))
    println(io,"Dependent Variables: ", get_depvars(sys))
    println(io,"Independent Variables: ", get_indvars(sys))
    println(io,"Parameters: ", get_ps(sys))
    print(io,"Default Parameter Values", get_defaults(sys))
    return nothing
end
