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
    """Dependent (state) variables. Must not contain the independent variable."""
    states::Vector
    """Parameter variables. Must not contain the independent variable."""
    ps::Vector
    """Control parameters (some subset of `ps`)."""
    ctrls::Vector
    """Observed states."""
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
    function DiscreteSystem(discreteEqs, iv, dvs, ps, ctrls, observed, name, systems, default_u0, default_p)
        check_variables(dvs,iv)
        check_parameters(ps,iv)
        new(discreteEqs, iv, dvs, ps, ctrls, observed, name, systems, default_u0, default_p)
    end
end

"""
    $(TYPEDSIGNATURES)

Constructs a DiscreteSystem.
"""
function DiscreteSystem(
                   discreteEqs::AbstractVector{<:Equation}, iv, dvs, ps;
                   controls = Num[],
                   observed = Num[],
                   systems = DiscreteSystem[],
                   name=gensym(:DiscreteSystem),
                   default_u0=Dict(),
                   default_p=Dict(),
                   defaults=_merge(Dict(default_u0), Dict(default_p)),
                  )
    iv′ = value(iv)
    dvs′ = value.(dvs)
    ps′ = value.(ps)
    ctrl′ = value.(controls)

    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn("`default_u0` and `default_p` are deprecated. Use `defaults` instead.", :ODESystem, force=true)
    end
    defaults = todict(defaults)
    defaults = Dict(value(k) => value(v) for (k, v) in pairs(defaults))

    collect_defaults!(defaults, dvs′)
    collect_defaults!(defaults, ps′)

    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end
    DiscreteSystem(discreteEqs, iv′, dvs′, ps′, ctrl′, observed, name, systems, default_u0, default_p)
end

"""
    $(TYPEDSIGNATURES)

Generates an DiscreteProblem from an DiscreteSystem.
"""
function DiffEqBase.DiscreteProblem(sys::DiscreteSystem,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    eval_module = @__MODULE__,
                                    eval_expression = true,
                                    kwargs...)
    dvs = states(sys)
    ps = parameters(sys)
    eqs = equations(sys)
    # defs = defaults(sys)
    t = get_iv(sys)
    u0 = varmap_to_vars(u0map,dvs)
    rhss = [eq.rhs for eq in eqs]
    u = dvs
    p = varmap_to_vars(parammap,ps)

    f_gen = generate_function(sys; expression=Val{eval_expression}, expression_module=eval_module)
    f_oop, _ = (@RuntimeGeneratedFunction(eval_module, ex) for ex in f_gen)
    f(u,p,t) = f_oop(u,p,t)
    DiscreteProblem(f,u0,tspan,p;kwargs...)
end

isdifference(expr) = istree(expr) && operation(expr) isa Difference
isdifferenceeq(eq) = isdifference(eq.lhs)

difference_vars(x::Sym) = Set([x])
difference_vars(exprs::Symbolic) = difference_vars([exprs])
difference_vars(exprs) = foldl(difference_vars!, exprs; init = Set())
difference_vars!(difference_vars, eq::Equation) = (difference_vars!(difference_vars, eq.lhs); difference_vars!(difference_vars, eq.rhs); difference_vars)
function difference_vars!(difference_vars, O)
    if isa(O, Sym)
        return push!(difference_vars, O)
    end
    !istree(O) && return difference_vars

    operation(O) isa Difference && return push!(difference_vars, O)

    if operation(O) === (getindex) &&
        first(arguments(O)) isa Symbolic

        return push!(difference_vars, O)
    end

    symtype(operation(O)) <: FnType && push!(difference_vars, O)
    for arg in arguments(O)
        difference_vars!(difference_vars, arg)
    end

    return difference_vars
end

function collect_difference_variables(sys::DiscreteSystem)
    eqs = equations(sys)
    vars = Set()
    difference_vars = Set()
    for eq in eqs
        difference_vars!(vars, eq)
        for v in vars
            isdifference(v) || continue
            push!(difference_vars, arguments(v)[1])
        end
        empty!(vars)
    end
    return difference_vars
end

@noinline function throw_invalid_difference(difvar, eq)
    msg = "The difference variable must be isolated to the left-hand " *
    "side of the equation like `$difvar ~ ...`.\n Got $eq."
    throw(InvalidSystemException(msg))
end

function check_difference_variables(eq, expr=eq.rhs)
    istree(expr) || return nothing
    if operation(expr) isa Difference
        throw_invalid_difference(expr, eq)
    end
    foreach(Base.Fix1(check_difference_variables, eq), arguments(expr))
end

function generate_function(
        sys::DiscreteSystem, dvs = states(sys), ps = parameters(sys);
        kwargs...
    )
    eqs = equations(sys)
    foreach(check_difference_variables, eqs)
    # substitute x(t) by just x
    rhss = [isdifference(eq.lhs) ? arguments(eq.lhs)[1] + eq.rhs : eq.rhs for eq in eqs]

    u = map(x->time_varying_as_func(value(x), sys), dvs)
    p = map(x->time_varying_as_func(value(x), sys), ps)
    t = get_iv(sys)
    
    build_function(rhss, u, p, t; kwargs...)
end
