export NonlinearSystem


struct NLEq
    rhs::Expression
end
function Base.convert(::Type{NLEq}, eq::Equation)
    isequal(eq.lhs, Constant(0)) || return NLEq(eq.rhs - eq.lhs)
    return NLEq(eq.rhs)
end
Base.convert(::Type{Equation}, eq::NLEq) = Equation(0, eq.rhs)
Base.:(==)(a::NLEq, b::NLEq) = a.rhs == b.rhs
get_args(eq::NLEq) = Expression[eq.rhs]

struct NonlinearSystem <: AbstractSystem
    eqs::Vector{NLEq}
    vs::Vector{Variable}
    ps::Vector{Variable}
end

function NonlinearSystem(eqs)
    vs, ps = extract_elements(eqs, [_subtype(:Unknown), _subtype(:Parameter)])
    NonlinearSystem(eqs, vs, ps)
end


function calculate_jacobian(sys::NonlinearSystem)
    rhs = [eq.rhs for eq in sys.eqs]
    jac = expand_derivatives.(calculate_jacobian(rhs, sys.vs))
    return jac
end

system_eqs(sys::NonlinearSystem) = collect(Equation, sys.eqs)
system_vars(sys::NonlinearSystem) = sys.vs
system_params(sys::NonlinearSystem) = sys.ps
