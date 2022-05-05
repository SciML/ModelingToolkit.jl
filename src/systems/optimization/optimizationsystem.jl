"""
$(TYPEDEF)

A scalar equation for optimization.

# Fields
$(FIELDS)

# Examples

```julia
@variables x y z
@parameters σ ρ β

op = σ*(y-x) + x*(ρ-z)-y + x*y - β*z
@named os = OptimizationSystem(op, [x,y,z],[σ,ρ,β])
```
"""
struct OptimizationSystem <: AbstractTimeIndependentSystem
    """Vector of equations defining the system."""
    op::Any
    """Unknown variables."""
    states::Vector
    """Parameters."""
    ps::Vector
    """Array variables."""
    var_to_name
    observed::Vector{Equation}
    equality_constraints::Vector{Equation}
    inequality_constraints::Vector
    """
    Name: the name of the system.  These are required to have unique names.
    """
    name::Symbol
    """
    systems: The internal systems
    """
    systems::Vector{OptimizationSystem}
    """
    defaults: The default values to use when initial conditions and/or
    parameters are not supplied in `ODEProblem`.
    """
    defaults::Dict
    function OptimizationSystem(op, states, ps, var_to_name, observed, equality_constraints, inequality_constraints, name, systems, defaults; checks::Bool = true)
        if checks
            check_units(op)
            check_units(observed)
            check_units(equality_constraints)
            all_dimensionless([states;ps]) || check_units(inequality_constraints)
        end
        new(op, states, ps, var_to_name, observed, equality_constraints, inequality_constraints, name, systems, defaults)
    end
end

function OptimizationSystem(op, states, ps;
                            observed = [],
                            equality_constraints = Equation[],
                            inequality_constraints = [],
                            default_u0=Dict(),
                            default_p=Dict(),
                            defaults=_merge(Dict(default_u0), Dict(default_p)),
                            name=nothing,
                            systems = OptimizationSystem[],
                            checks = true)
    name === nothing && throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn("`default_u0` and `default_p` are deprecated. Use `defaults` instead.", :OptimizationSystem, force=true)
    end
    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end
    defaults = todict(defaults)
    defaults = Dict(value(k) => value(v) for (k, v) in pairs(defaults))

    states, ps = value.(states), value.(ps)
    var_to_name = Dict()
    process_variables!(var_to_name, defaults, states)
    process_variables!(var_to_name, defaults, ps)
    isempty(observed) || collect_var_to_name!(var_to_name, (eq.lhs for eq in observed))
    OptimizationSystem(
                       value(op), states, ps, var_to_name,
                       observed,
                       equality_constraints, inequality_constraints,
                       name, systems, defaults; checks = checks
                      )
end

function calculate_jacobian(sys::OptimizationSystem)
    expand_derivatives.(jacobian(equations(sys), states(sys)))
end

function generate_jacobian(sys::OptimizationSystem, vs = states(sys), ps = parameters(sys); kwargs...)
    jac = calculate_jacobian(sys)
    return build_function(jac, vs, ps; conv =  AbstractSysToExpr(sys), kwargs...)
end

function calculate_gradient(sys::OptimizationSystem)
    expand_derivatives.(gradient(equations(sys), states(sys)))
end

function generate_gradient(sys::OptimizationSystem, vs = states(sys), ps = parameters(sys); kwargs...)
    grad = calculate_gradient(sys)
    return build_function(grad, vs, ps;
                          conv = AbstractSysToExpr(sys),kwargs...)
end

function calculate_hessian(sys::OptimizationSystem)
    expand_derivatives.(hessian(equations(sys), states(sys)))
end

function generate_hessian(sys::OptimizationSystem, vs = states(sys), ps = parameters(sys);
                          sparse = false, kwargs...)
    if sparse
        hess = sparsehessian(equations(sys),states(sys))
    else
        hess = calculate_hessian(sys)
    end
    return build_function(hess, vs, ps;
                          conv = AbstractSysToExpr(sys),kwargs...)
end

function generate_function(sys::OptimizationSystem, vs = states(sys), ps = parameters(sys); kwargs...)
    return build_function(equations(sys), vs, ps;
                          conv = AbstractSysToExpr(sys),kwargs...)
end

equations(sys::OptimizationSystem) = isempty(get_systems(sys)) ? get_op(sys) : get_op(sys) + reduce(+,namespace_expr.(get_systems(sys)))
namespace_expr(sys::OptimizationSystem) = namespace_expr(get_op(sys), sys)

hessian_sparsity(sys::OptimizationSystem) = hessian_sparsity(get_op(sys), states(sys))
