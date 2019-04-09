export NonlinearSystem


struct NLEq
    rhs::Expression
end
function Base.convert(::Type{NLEq}, eq::Equation)
    isequal(eq.lhs, Constant(0)) || return NLEq(eq.rhs - eq.lhs)
    return NLEq(eq.rhs)
end
Base.:(==)(a::NLEq, b::NLEq) = a.rhs == b.rhs

struct NonlinearSystem <: AbstractSystem
    eqs::Vector{NLEq}
    vs::Vector{Expression}
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


function calculate_jacobian(sys::NonlinearSystem)
    rhs = [eq.rhs for eq in sys.eqs]
    jac = expand_derivatives.(calculate_jacobian(rhs, sys.vs))
    return jac
end

function generate_jacobian(sys::NonlinearSystem; version::FunctionVersion = ArrayFunction)
    jac = calculate_jacobian(sys)
    return build_function(jac, clean.(sys.vs), sys.ps; version = version)
end

function generate_function(sys::NonlinearSystem, vs, ps; version::FunctionVersion = ArrayFunction)
    rhss = [eq.rhs for eq ∈ sys.eqs]
    vs′ = [clean(v) for v ∈ vs]
    ps′ = [clean(p) for p ∈ ps]
    return build_function(rhss, vs′, ps′; version = version)
end
