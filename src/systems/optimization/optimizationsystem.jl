"""
$(TYPEDEF)

A scalar equation for optimization.

# Fields
* `op` - The objective function

# Examples

```
@variables x y z
@parameters σ ρ β

op = σ*(y-x) + x*(ρ-z)-y + x*y - β*z
os = OptimizationSystem(eqs, [x,y,z],[σ,ρ,β])
```
"""
struct OptimizationSystem <: AbstractSystem
    """Vector of equations defining the system."""
    op::Operation
    """Unknown variables."""
    states::Vector{Expression}
    """Parameters."""
    ps::Vector{Variable}
    """
    Name: the name of the system
    """
    name::Symbol
    """
    systems: The internal systems
    """
    systems::Vector{NonlinearSystem}
end

function NonlinearSystem(eqs, states, ps;
                         name = gensym(:NonlinearSystem),
                         systems = NonlinearSystem[])
    NonlinearSystem(eqs, states, convert.(Variable,ps), name, systems)
end

function calculate_jacobian(sys::NonlinearSystem)
    rhs = [eq.rhs for eq in sys.eqs]
    jac = expand_derivatives.(calculate_jacobian(rhs, sys.states))
    return jac
end

function generate_jacobian(sys::NonlinearSystem, vs = sys.states, ps = sys.ps, expression = Val{true}; kwargs...)
    jac = calculate_jacobian(sys)
    return build_function(jac, convert.(Variable,vs), convert.(Variable,ps), (), NLSysToExpr(sys))
end

struct NLSysToExpr
    sys::NonlinearSystem
end
function (f::NLSysToExpr)(O::Operation)
    any(isequal(O), f.sys.states) && return O.op.name  # variables
    if isa(O.op, Variable)
        isempty(O.args) && return O.op.name  # 0-ary parameters
        return build_expr(:call, Any[O.op.name; f.(O.args)])
    end
    return build_expr(:call, Any[O.op; f.(O.args)])
end
(f::NLSysToExpr)(x) = convert(Expr, x)

function generate_function(sys::NonlinearSystem, vs = sys.states, ps = sys.ps, expression = Val{true}; kwargs...)
    rhss = [eq.rhs for eq ∈ sys.eqs]
    vs′ = convert.(Variable,vs)
    ps′ = convert.(Variable,ps)
    return build_function(rhss, vs′, ps′, (), NLSysToExpr(sys), expression; kwargs...)
end
