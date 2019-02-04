export NonlinearSystem


struct NLEq
    rhs::Expression
end
function Base.convert(::Type{NLEq}, eq::Equation)
    isequal(eq.lhs, Constant(0)) || return NLEq(eq.rhs - eq.lhs)
    return NLEq(eq.rhs)
end
Base.:(==)(a::NLEq, b::NLEq) = a.rhs == b.rhs
get_args(eq::NLEq) = Expression[eq.rhs]

struct NonlinearSystem <: AbstractSystem
    eqs::Vector{NLEq}
    vs::Vector{Variable}
    ps::Vector{Variable}
end

function NonlinearSystem(eqs)
    vs, ps = extract_elements(eqs, [_is_unknown, _is_known])
    NonlinearSystem(eqs, vs, ps)
end


function calculate_jacobian(sys::NonlinearSystem)
    rhs = [eq.rhs for eq in sys.eqs]
    jac = expand_derivatives.(calculate_jacobian(rhs, sys.vs))
    return jac
end

function generate_jacobian(sys::NonlinearSystem; version::FunctionVersion = ArrayFunction)
    jac = calculate_jacobian(sys)
    return build_function(jac, sys.vs, sys.ps; version = version)
end

function generate_function(sys::NonlinearSystem; version::FunctionVersion = ArrayFunction)
    rhss = [eq.rhs for eq ∈ sys.eqs]
    return build_function(rhss, sys.vs, sys.ps; version = version)
end
