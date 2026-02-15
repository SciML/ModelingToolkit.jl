"""
$(TYPEDEF)

A system of partial differential equations.

# Fields
$(FIELDS)

# Example

```julia
using ModelingToolkit

@parameters x t
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
domains = [t ∈ (0.0,1.0),
           x ∈ (0.0,1.0)]

@named pde_system = PDESystem(eq,bcs,domains,[t,x],[u])
```
"""
struct PDESystem <: AbstractSystem
    "The equations which define the PDE."
    eqs::Any
    "The boundary conditions."
    bcs::Any
    "The domain for the independent variables."
    domain::Any
    "The independent variables."
    ivs::Any
    "The dependent variables."
    dvs::Any
    "The parameters."
    ps::Any
    """
    Initial conditions for variables (unknowns/observables/parameters) which can be
    changed/overridden. When constructing a numerical problem from the system.
    """
    initial_conditions::SymmapT
    """
    Type of the system.
    """
    connector_type::Any
    """
    The internal systems. These are required to have unique names.
    """
    systems::Vector
    """
    A vector of explicit symbolic expressions for the analytic solutions of each
    dependent variable. e.g. `analytic = [u(t, x) ~ a*sin(c*t) * cos(k*x)]`.
    """
    analytic::Any
    """
    A vector of functions for the analytic solutions of each dependent
    variable. Will be generated from `analytic` if not provided. Should have the same
    argument signature as the variable, and a `ps` argument as the last argument,
    which takes an indexable of parameter values in the order you specified them in `ps`.
    e.g. `analytic_func = [u(t, x) => (ps, t, x) -> ps[1]*sin(ps[2]*t) * cos(ps[3]*x)]`.
    """
    analytic_func::Any
    """
    The name of the system.
    """
    name::Symbol
    """
    A description of the system.
    """
    description::String
    """
    Metadata for the system, to be used by downstream packages.
    """
    metadata::Any
    """
    Metadata for MTK GUI.
    """
    gui_metadata::Union{Nothing, GUIMetadata}
    """
    Mapping from variable names to symbolic variables, enabling the `sys.x` accessor pattern.
    """
    var_to_name::Dict{Symbol, SymbolicT}
    @add_kwonly function PDESystem(
            eqs, bcs, domain, ivs, dvs,
            ps = SciMLBase.NullParameters();
            initial_conditions = Dict(),
            systems = [],
            connector_type = nothing,
            metadata = nothing,
            analytic = nothing,
            analytic_func = nothing,
            gui_metadata = nothing,
            eval_module = @__MODULE__,
            checks::Union{Bool, Int} = true,
            description = "",
            name
        )
        if checks == true || (checks & CheckUnits) > 0
            u = __get_unit_type(dvs, ivs, ps)
            check_units(u, eqs)
        end

        eqs = eqs isa Vector ? eqs : [eqs]

        if !isnothing(analytic)
            analytic = analytic isa Vector ? analytic : [analytic]
            if length(analytic) != length(dvs)
                throw(ArgumentError("The number of analytic solutions must match the number of dependent variables"))
            end

            if isnothing(analytic_func)
                analytic_func = map(analytic) do eq
                    args = arguments(eq.lhs)
                    p = ps isa SciMLBase.NullParameters ? [] : ps
                    args = vcat(DestructuredArgs(p), args)
                    ex = Func(args, [], eq.rhs) |> toexpr
                    eq.lhs => drop_expr(
                        RuntimeGeneratedFunction(
                            eval_module, eval_module, ex
                        )
                    )
                end
            end
        end

        if !isnothing(analytic_func)
            analytic_func = analytic_func isa Dict ? analytic_func : analytic_func |> Dict
        end

        var_to_name = Dict{Symbol, SymbolicT}()
        dvs_vec = SymbolicT[unwrap(v) for v in dvs]
        collect_var_to_name!(var_to_name, dvs_vec)
        if !(ps isa SciMLBase.NullParameters)
            ps_vec = SymbolicT[unwrap(p) for p in ps]
            collect_var_to_name!(var_to_name, ps_vec)
        end

        new(
            eqs, bcs, domain, ivs, dvs, ps, initial_conditions, connector_type, systems, analytic,
            analytic_func, name, description, metadata, gui_metadata, var_to_name
        )
    end
end

function Base.getproperty(
        x::PDESystem, sym::Symbol; namespace = does_namespacing(x)
    )
    if sym == :indvars
        Base.depwarn(
            "`sys.indvars` is deprecated, please use `get_ivs(sys)`", :getproperty,
            force = true
        )
        return getfield(x, :ivs)
    elseif sym == :depvars
        Base.depwarn(
            "`sys.depvars` is deprecated, please use `get_dvs(sys)`", :getproperty,
            force = true
        )
        return getfield(x, :dvs)
    elseif hasfield(PDESystem, sym)
        return getfield(x, sym)
    end
    # Check independent variables first (IVs are not namespaced, matching System behavior)
    for iv in getfield(x, :ivs)
        if getname(unwrap(iv)) == sym
            return iv
        end
    end
    # Fall through to AbstractSystem's getvar for variable lookup (sys.x pattern)
    return wrap(getvar(x, sym; namespace = namespace))
end

# Accessor bridges: PDESystem stores dependent variables as `dvs` and independent
# variables as `ivs` (plural), whereas the AbstractSystem API expects `unknowns` and `iv`.
get_unknowns(sys::PDESystem) = getfield(sys, :dvs)
has_unknowns(::PDESystem) = true

function Base.propertynames(sys::PDESystem; private = false)
    if private
        return fieldnames(PDESystem)
    else
        names = Symbol[]
        for s in get_systems(sys)
            push!(names, getname(s))
        end
        for s in get_dvs(sys)
            push!(names, getname(s))
        end
        ps = get_ps(sys)
        if !(ps isa SciMLBase.NullParameters)
            for s in ps
                hasname(s) || continue
                push!(names, getname(s))
            end
        end
        for s in get_ivs(sys)
            push!(names, getname(s))
        end
        return names
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
    print(io, "Default Parameter Values", get_initial_conditions(sys))
    return nothing
end
