"""
$(TYPEDEF)

A system of partial differential equations.

# Fields
$(FIELDS)

# Example

```julia
using ModelingToolkit

@parameters x
@variables t u(..)
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
domains = [t ∈ (0.0,1.0),
           x ∈ (0.0,1.0)]

@named pde_system = PDESystem(eq,bcs,domains,[t,x],[u])
```
"""
struct PDESystem <: ModelingToolkit.AbstractMultivariateSystem
    "The equations which define the PDE"
    eqs::Any
    "The boundary conditions"
    bcs::Any
    "The domain for the independent variables."
    domain::Any
    "The independent variables"
    ivs::Any
    "The dependent variables"
    dvs::Any
    "The parameters"
    ps::Any
    """
    defaults: The default values to use when initial conditions and/or
    parameters are not supplied in `ODEProblem`.
    """
    defaults::Dict
    """
    type: type of the system
    """
    connector_type::Any
    """
    systems: The internal systems. These are required to have unique names.
    """
    systems::Vector
    """
    name: the name of the system
    """
    name::Symbol
    """
    metadata: metadata for the system, to be used by downstream packages.
    """
    metadata::Any
    @add_kwonly function PDESystem(eqs, bcs, domain, ivs, dvs,
                                   ps = SciMLBase.NullParameters();
                                   defaults = Dict(),
                                   systems = [],
                                   connector_type = nothing,
                                   checks::Bool = true,
                                   name, metadata = nothing)
        if checks
            all_dimensionless([dvs; ivs; ps]) || check_units(eqs)
        end
        eqs = eqs isa Vector ? eqs : [eqs]
        new(eqs, bcs, domain, ivs, dvs, ps, defaults, connector_type, systems, name, metadata)
    end
end

function Base.getproperty(x::PDESystem, sym::Symbol)
    if sym == :indvars
        return getfield(x, :ivs)
        Base.depwarn("`sys.indvars` is deprecated, please use `get_ivs(sys)`", :getproperty,
                     force = true)

    elseif sym == :depvars
        return getfield(x, :dvs)
        Base.depwarn("`sys.depvars` is deprecated, please use `get_dvs(sys)`", :getproperty,
                     force = true)

    else
        return getfield(x, sym)
    end
end

Base.summary(prob::PDESystem) = string(nameof(typeof(prob)))
function Base.show(io::IO, ::MIME"text/plain", sys::PDESystem)
    println(io, summary(sys))
    println(io, "Equations: ", get_eqs(sys))
    println(io, "Boundary Conditions: ", get_bcs(sys))
    println(io, "Domain: ", get_domain(sys))
    println(io, "Dependent Variables: ", get_dvs(sys))
    println(io, "Independent Variables: ", get_ivs(sys))
    println(io, "Parameters: ", get_ps(sys))
    print(io, "Default Parameter Values", get_defaults(sys))
    return nothing
end
