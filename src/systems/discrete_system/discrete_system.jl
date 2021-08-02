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
    """Array variables."""
    var_to_name
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
    function DiscreteSystem(discreteEqs, iv, dvs, ps, var_to_name, ctrls, observed, name, systems, default_u0, default_p)
        check_variables(dvs,iv)
        check_parameters(ps,iv)
        new(discreteEqs, iv, dvs, ps, var_to_name, ctrls, observed, name, systems, default_u0, default_p)
    end
end

"""
    $(TYPEDSIGNATURES)

Constructs a DiscreteSystem.
"""
function DiscreteSystem(
                   eqs::AbstractVector{<:Equation}, iv, dvs, ps;
                   controls = Num[],
                   observed = Num[],
                   systems = DiscreteSystem[],
                   name=gensym(:DiscreteSystem),
                   default_u0=Dict(),
                   default_p=Dict(),
                   defaults=_merge(Dict(default_u0), Dict(default_p)),
                  )
    eqs = collect(eqs)
    iv′ = value(iv)
    dvs′ = value.(dvs)
    ps′ = value.(ps)
    ctrl′ = value.(controls)

    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn("`default_u0` and `default_p` are deprecated. Use `defaults` instead.", :ODESystem, force=true)
    end
    defaults = todict(defaults)
    defaults = Dict(value(k) => value(v) for (k, v) in pairs(defaults))

    var_to_name = Dict()
    process_variables!(var_to_name, defaults, dvs′)
    process_variables!(var_to_name, defaults, ps′)

    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end
    DiscreteSystem(eqs, iv′, dvs′, ps′, var_to_name, ctrl′, observed, name, systems, default_u0, default_p)
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
    eqs, max_delay = linearize_eqs(sys, eqs)
    # TODO: Add equations to update delayed variables
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

function linearize_eqs(sys, eqs=sys.eqs)
    max_delay = 0
    r = @rule ~t => begin
        if ~t isa Symbolics.Term && any(isequal((~t).f), Symbolics.operation.(sys.states)) && is_delay(sys, Symbolics.arguments(~t))
            delay = get_delay_val(sys, first(Symbolics.arguments(~t)))
            if delay > max_delay
                max_delay = delay
            end
            var_sym = Symbol((~t).f.name, Symbol("("), join(Symbolics.tosymbol.(Symbolics.arguments(~t)), Symbol(",")), Symbol(")"))
            Symbolics.variable(var_sym)
        else
            ~t
        end
    end
    rw = SymbolicUtils.Postwalk(r)
    [
        Symbolics.value(eq.lhs) ~ Symbolics.value(rw(eq.rhs))
        for eq in eqs
    ], max_delay
end

function get_delay_val(sys, x)
    delay = x - sys.iv
    delay > 0 && error("Forward delay not permitted")
    return -delay
end

function is_delay(sys, args::AbstractVector)
    length(args) > 1 && return false
    isequal(first(args), sys.iv) && return false
    delay = sys.iv - first(args)
    delay isa Integer || 
    (delay isa AbstractFloat && isinteger(delay)) ||
    (delay isa Num && isinteger(delay.val))
    
end

check_difference_variables(eq) = check_operator_variables(eq, Difference)

function generate_function(
        sys::DiscreteSystem, dvs = states(sys), ps = parameters(sys);
        kwargs...
    )
    eqs = equations(sys)
    foreach(check_difference_variables, eqs)
    # substitute x(t) by just x
    rhss = [eq.rhs for eq in eqs]

    u = map(x->time_varying_as_func(value(x), sys), dvs)
    p = map(x->time_varying_as_func(value(x), sys), ps)
    t = get_iv(sys)
    
    build_function(rhss, u, p, t; kwargs...)
end
