"""
$(TYPEDEF)

A system of difference equations.

# Fields
$(FIELDS)

# Example

```
using ModelingToolkit

@parameters t σ ρ β
@variables x(t) y(t) z(t) next_x(t) next_y(t) next_z(t)

eqs = [next_x ~ σ*(y-x),
       next_y ~ x*(ρ-z)-y,
       next_z ~ x*y - β*z]

de = DiscreteSystem(eqs,t,[x,y,z],[σ,ρ,β])
```
"""
struct DiscreteSystem <: AbstractSystem
    """The differential equations defining the discrete system."""
    eqs::Vector{Equation}
    """Independent variable."""
    iv::Sym
    """Dependent (state) variables."""
    states::Vector
    """Parameter variables."""
    ps::Vector
    observed::Vector{Equation}
    """
    Name: the name of the system
    """
    name::Symbol
    """
    systems: The internal systems. These are required to have unique names.
    """
    systems::Vector{DiscreteSystem}
    """
    default_u0: The default initial conditions to use when initial conditions
    are not supplied in `DiscreteSystem`.
    """
    default_u0::Dict
    """
    default_p: The default parameters to use when parameters are not supplied
    in `DiscreteSystem`.
    """
    default_p::Dict
end

"""
    $(TYPEDSIGNATURES)

Constructs a DiscreteSystem.
"""
function DiscreteSystem(
                   discreteEqs::AbstractVector{<:Equation}, iv, dvs, ps;
                   observed = Num[],
                   systems = DiscreteSystem[],
                   name=gensym(:DiscreteSystem),
                   default_u0=Dict(),
                   default_p=Dict(),
                  )
    iv′ = value(iv)
    dvs′ = value.(dvs)
    ps′ = value.(ps)

    default_u0 isa Dict || (default_u0 = Dict(default_u0))
    default_p isa Dict || (default_p = Dict(default_p))
    default_u0 = Dict(value(k) => value(default_u0[k]) for k in keys(default_u0))
    default_p = Dict(value(k) => value(default_p[k]) for k in keys(default_p))

    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end
    DiscreteSystem(discreteEqs, iv′, dvs′, ps′, observed, name, systems, default_u0, default_p)
end

"""
    $(TYPEDSIGNATURES)

Generates an DiscreteProblem from an DiscreteSystem.
"""
function DiffEqBase.DiscreteProblem(sys::DiscreteSystem,u0map,tspan,parammap=DiffEqBase.NullParameters();kwargs...)
    dvs = states(sys)
    ps = parameters(sys)
    eqs = equations(sys)
    # defs = defaults(sys)
    t = get_iv(sys)
    u0 = varmap_to_vars(u0map,dvs)
    rhss = [eq.rhs for eq in eqs]
    u = dvs
    p = varmap_to_vars(parammap,ps)
    eval_module = @__MODULE__
    eval_expression = true
    f_gen = build_function(rhss, dvs, ps, t; expression=Val{eval_expression}, expression_module=eval_module)
    f_oop,f_iip = (@RuntimeGeneratedFunction(eval_module, ex) for ex in f_gen)
    f(u,p,t) = f_oop(u,p,t)
    DiscreteProblem(f,u0,tspan,p)
end
