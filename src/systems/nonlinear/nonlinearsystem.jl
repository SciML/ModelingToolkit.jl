"""
$(TYPEDEF)

A nonlinear system of equations.

# Fields
* `eqs` - Vector of equations defining the system.

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
    states::Vector{Variable}
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
    NonlinearSystem(eqs, convert.(Variable,states), convert.(Variable,ps), name, systems)
end

function calculate_jacobian(sys::NonlinearSystem)
    rhs = [eq.rhs for eq in sys.eqs]
    jac = expand_derivatives.(calculate_jacobian(rhs, [dv() for dv in states(sys)]))
    return jac
end

function generate_jacobian(sys::NonlinearSystem, vs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)
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

function generate_function(sys::NonlinearSystem, vs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)
    rhss = [eq.rhs for eq ∈ sys.eqs]
    vs′ = convert.(Variable,vs)
    ps′ = convert.(Variable,ps)
    return build_function(rhss, vs′, ps′, (), NLSysToExpr(sys), expression; kwargs...)
end

function DiffEqBase.NonlinearProblem{iip}(sys::NonlinearSystem,u0map,tspan,
                                          parammap=DiffEqBase.NullParameters();
                                          jac = false,
                                          checkbounds = false,
                                          linenumbers = true, multithread=false,
                                          kwargs...) where iip
    dvs = states(sys)
    ps = parameters(sys)

    f = generate_function(sys;checkbounds=checkbounds,linenumbers=linenumbers,
                              multithread=multithread)
    u0 = varmap_to_vars(u0map,dvs)
    p = varmap_to_vars(parammap,ps)
    NonlinearProblem(f,u0,tspan,p;kwargs...)
end
