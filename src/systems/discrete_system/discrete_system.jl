"""
$(TYPEDEF)

A system of difference equations.

# Fields
$(FIELDS)

# Example

```
using ModelingToolkit

@parameters σ ρ β
@variables t x(t) y(t) z(t) next_x(t) next_y(t) next_z(t)

eqs = [next_x ~ σ*(y-x),
       next_y ~ x*(ρ-z)-y,
       next_z ~ x*y - β*z]

@named de = DiscreteSystem(eqs,t,[x,y,z],[σ,ρ,β])
```
"""
struct DiscreteSystem <: AbstractTimeDependentSystem
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
    function DiscreteSystem(discreteEqs, iv, dvs, ps, var_to_name, ctrls, observed, name, systems, default_u0, default_p; checks::Bool = true)
        if checks
            check_variables(dvs, iv)
            check_parameters(ps, iv)
            all_dimensionless([dvs;ps;iv;ctrls]) ||check_units(discreteEqs)
        end
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
                   name=nothing,
                   default_u0=Dict(),
                   default_p=Dict(),
                   defaults=_merge(Dict(default_u0), Dict(default_p)),
                   kwargs...,
                  )
    name === nothing && throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
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
    DiscreteSystem(eqs, iv′, dvs′, ps′, var_to_name, ctrl′, observed, name, systems, default_u0, default_p, kwargs...)
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
    eqs = linearize_eqs(sys, eqs)
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

function linearize_eqs(sys, eqs=get_eqs(sys); return_max_delay=false)
    unique_states = unique(operation.(states(sys)))
    max_delay = Dict(v=>0.0 for v in unique_states)

    r = @rule ~t::(t -> istree(t) && any(isequal(operation(t)), operation.(states(sys))) && is_delay_var(get_iv(sys), t)) => begin
        delay = get_delay_val(get_iv(sys), first(arguments(~t)))
        if delay > max_delay[operation(~t)]
            max_delay[operation(~t)] = delay
        end
        nothing
    end
    SymbolicUtils.Postwalk(r).(rhss(eqs))

    if any(values(max_delay) .> 0)

        dts = Dict(v=>Any[] for v in unique_states)
        state_ops = Dict(v=>Any[] for v in unique_states)
        for v in unique_states
            for eq in eqs
                if isdifferenceeq(eq) && istree(arguments(eq.lhs)[1]) && isequal(v, operation(arguments(eq.lhs)[1]))
                    append!(dts[v], [operation(eq.lhs).dt])
                    append!(state_ops[v], [operation(eq.lhs)])
                end
            end
        end

        all(length.(unique.(values(state_ops))) .<= 1) || error("Each state should be used with single difference operator.")
        
        dts_gcd = Dict()
        for v in keys(dts)
            dts_gcd[v] = (length(dts[v]) > 0) ? first(dts[v]) : nothing
        end

        lin_eqs = [
            v(get_iv(sys) - (t)) ~ v(get_iv(sys) - (t-dts_gcd[v]))
            for v in unique_states if max_delay[v] > 0 && dts_gcd[v]!==nothing for t in collect(max_delay[v]:(-dts_gcd[v]):0)[1:end-1] 
        ]
        eqs = vcat(eqs, lin_eqs)
    end
    if return_max_delay return eqs, max_delay end
    eqs
end

function get_delay_val(iv, x)
    delay = x - iv
    isequal(delay > 0, true) && error("Forward delay not permitted")
    return -delay
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
