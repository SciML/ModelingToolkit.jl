export NonlinearSystem


struct NLEq
    rhs::Expression
end
function Base.convert(::Type{NLEq}, eq::Equation)
    isequal(eq.lhs, Constant(0)) || return NLEq(eq.rhs - eq.lhs)
    return NLEq(eq.rhs)
end
Base.:(==)(a::NLEq, b::NLEq) = a.rhs == b.rhs

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
ns = NonlinearSystem(eqs, [x,y,z])
```
"""
struct NonlinearSystem <: AbstractSystem
    """Vector of equations defining the system."""
    eqs::Vector{NLEq}
    """Unknown variables."""
    vs::Vector{Expression}
    """Parameters."""
    ps::Vector{Variable}
    function NonlinearSystem(eqs, vs)
        rhss = [eq.rhs for eq ∈ eqs]
        ps = reduce(∪, map(_find_params(vs), rhss); init = vnil())
        new(eqs, vs, collect(ps))
    end
end

vnil() = Set{Variable}()
_find_params(vs) = Base.Fix2(_find_params, vs)
function _find_params(O, vs)
    isa(O, Operation) || return vnil()
    any(isequal(O), vs) && return vnil()
    ps = reduce(∪, map(_find_params(vs), O.args); init = vnil())
    isa(O.op, Variable) && push!(ps, O.op)
    return ps
end


independent_variables(sys::NonlinearSystem) = Set{Variable}()
dependent_variables(sys::NonlinearSystem) = Set{Expression}(sys.vs)
parameters(sys::NonlinearSystem) = Set{Variable}(sys.ps)


function calculate_jacobian(sys::NonlinearSystem)
    rhs = [eq.rhs for eq in sys.eqs]
    jac = expand_derivatives.(calculate_jacobian(rhs, sys.vs))
    return jac
end

function generate_jacobian(sys::NonlinearSystem; version::FunctionVersion = ArrayFunction)
    jac = calculate_jacobian(sys)
    return build_function(jac, clean.(sys.vs), sys.ps, (), NLSysToExpr(sys); version = version)
end

struct NLSysToExpr
    sys::NonlinearSystem
end
function (f::NLSysToExpr)(O::Operation)
    any(isequal(O), f.sys.vs) && return O.op.name  # variables
    if isa(O.op, Variable)
        isempty(O.args) && return O.op.name  # 0-ary parameters
        return build_expr(:call, Any[O.op.name; f.(O.args)])
    end
    return build_expr(:call, Any[O.op; f.(O.args)])
end
(f::NLSysToExpr)(x) = convert(Expr, x)


function generate_function(sys::NonlinearSystem, vs, ps; version::FunctionVersion = ArrayFunction)
    rhss = [eq.rhs for eq ∈ sys.eqs]
    vs′ = [clean(v) for v ∈ vs]
    ps′ = [clean(p) for p ∈ ps]
    return build_function(rhss, vs′, ps′, (), NLSysToExpr(sys); version = version)
end
