"""
$(TYPEDEF)

A nonlinear system of equations.

# Fields
$(FIELDS)

# Examples

```
@variables x y z
@parameters σ ρ β

eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs, [x,y,z],[σ,ρ,β])
```
"""
struct NonlinearSystem <: AbstractSystem
    """Vector of equations defining the system."""
    eqs::Vector{Equation}
    """Unknown variables."""
    states::Vector
    """Parameters."""
    ps::Vector
    pins::Vector
    observed::Vector{Equation}
    """
    Name: the name of the system
    """
    name::Symbol
    """
    systems: The internal systems
    """
    systems::Vector{NonlinearSystem}
    """
    default_u0: The default initial conditions to use when initial conditions
    are not supplied in `ODEProblem`.
    """
    default_u0::Dict
    """
    default_p: The default parameters to use when parameters are not supplied
    in `ODEProblem`.
    """
    default_p::Dict
end

function NonlinearSystem(eqs, states, ps;
                         pins = [],
                         observed = [],
                         name = gensym(:NonlinearSystem),
                         default_u0=Dict(),
                         default_p=Dict(),
                         systems = NonlinearSystem[])
    default_u0 isa Dict || (default_u0 = Dict(default_u0))
    default_p isa Dict || (default_p = Dict(default_p))
    NonlinearSystem(eqs, value.(states), value.(ps), value.(pins), observed, name, systems, default_u0, default_p)
end

function calculate_jacobian(sys::NonlinearSystem;sparse=false,simplify=false)
    rhs = [eq.rhs for eq ∈ equations(sys)]
    vals = [dv for dv in states(sys)]
    if sparse
        jac = sparsejacobian(rhs, vals, simplify=simplify)
    else
        jac = jacobian(rhs, vals, simplify=simplify)
    end
    return jac
end

function generate_jacobian(sys::NonlinearSystem, vs = states(sys), ps = parameters(sys);
                           sparse = false, simplify=false, kwargs...)
    jac = calculate_jacobian(sys,sparse=sparse, simplify=simplify)
    return build_function(jac, vs, ps;
                          conv = AbstractSysToExpr(sys), kwargs...)
end

function generate_function(sys::NonlinearSystem, vs = states(sys), ps = parameters(sys); kwargs...)
    rhss = [eq.rhs for eq ∈ sys.eqs]
    return build_function(rhss, vs, ps;
                          conv = AbstractSysToExpr(sys), kwargs...)
end

jacobian_sparsity(sys::NonlinearSystem) =
    jacobian_sparsity([eq.rhs for eq ∈ equations(sys)],
                      states(sys))

"""
```julia
function DiffEqBase.NonlinearProblem{iip}(sys::NonlinearSystem,u0map,
                                          parammap=DiffEqBase.NullParameters();
                                          jac = false, sparse=false,
                                          checkbounds = false,
                                          linenumbers = true, parallel=SerialForm(),
                                          kwargs...) where iip
```

Generates an NonlinearProblem from a NonlinearSystem and allows for automatically
symbolically calculating numerical enhancements.
"""
function DiffEqBase.NonlinearProblem{iip}(sys::NonlinearSystem,u0map,
                                          parammap=DiffEqBase.NullParameters();
                                          jac = false, sparse=false,
                                          checkbounds = false,
                                          linenumbers = true, parallel=SerialForm(),
                                          kwargs...) where iip
    dvs = states(sys)
    ps = parameters(sys)

    f = generate_function(sys;checkbounds=checkbounds,linenumbers=linenumbers,
                              parallel=parallel,sparse=sparse,expression=Val{false})
    u0 = varmap_to_vars(u0map,dvs; defaults=get_default_u0(sys))
    p = varmap_to_vars(parammap,ps; defaults=get_default_p(sys))
    NonlinearProblem(f,u0,p;kwargs...)
end

"""
```julia
function DiffEqBase.NonlinearProblemExpr{iip}(sys::NonlinearSystem,u0map,
                                          parammap=DiffEqBase.NullParameters();
                                          jac = false, sparse=false,
                                          checkbounds = false,
                                          linenumbers = true, parallel=SerialForm(),
                                          kwargs...) where iip
```

Generates a Julia expression for a NonlinearProblem from a
NonlinearSystem and allows for automatically symbolically calculating
numerical enhancements.
"""
struct NonlinearProblemExpr{iip} end

function NonlinearProblemExpr{iip}(sys::NonlinearSystem,u0map,
                                          parammap=DiffEqBase.NullParameters();
                                          jac = false, sparse=false,
                                          checkbounds = false,
                                          linenumbers = false, parallel=SerialForm(),
                                          kwargs...) where iip
    dvs = states(sys)
    ps = parameters(sys)

    f = generate_function(sys;checkbounds=checkbounds,linenumbers=linenumbers,
                              parallel=parallel,sparse=sparse,expression=Val{true})
    u0 = varmap_to_vars(u0map,dvs; defaults=get_default_u0(sys))
    p = varmap_to_vars(parammap,ps; defaults=get_default_p(sys))
    quote
        f = $f
        u0 = $u0
        p = $p
        NonlinearProblem(f,u0,p;kwargs...)
    end
end
