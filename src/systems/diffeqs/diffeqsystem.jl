export DiffEqSystem, ODEFunction


using Base: RefValue


isintermediate(eq::Equation) = !(isa(eq.lhs, Operation) && isa(eq.lhs.op, Differential))

function flatten_differential(O::Operation)
    @assert is_derivative(O) "invalid differential: $O"
    is_derivative(O.args[1]) || return (O.args[1], O.op.x, 1)
    (x, t, order) = flatten_differential(O.args[1])
    isequal(t, O.op.x) || throw(ArgumentError("non-matching differentials on lhs: $t, $(O.op.x)"))
    return (x, t, order + 1)
end


struct DiffEq  # dⁿx/dtⁿ = rhs
    x::Expression
    t::Variable
    n::Int
    rhs::Expression
end
function Base.convert(::Type{DiffEq}, eq::Equation)
    isintermediate(eq) && throw(ArgumentError("intermediate equation received"))
    (x, t, n) = flatten_differential(eq.lhs)
    return DiffEq(x, t, n, eq.rhs)
end
Base.:(==)(a::DiffEq, b::DiffEq) = isequal((a.x, a.t, a.n, a.rhs), (b.x, b.t, b.n, b.rhs))
get_args(eq::DiffEq) = Expression[eq.x, eq.t, eq.rhs]

struct DiffEqSystem <: AbstractSystem
    eqs::Vector{DiffEq}
    iv::Variable
    dvs::Vector{Variable}
    ps::Vector{Variable}
    jac::RefValue{Matrix{Expression}}
    function DiffEqSystem(eqs, iv, dvs, ps)
        jac = RefValue(Matrix{Expression}(undef, 0, 0))
        new(eqs, iv, dvs, ps, jac)
    end
end

function DiffEqSystem(eqs)
    dvs, = extract_elements(eqs, [_is_dependent])
    ivs = unique(vcat((dv.dependents for dv ∈ dvs)...))
    length(ivs) == 1 || throw(ArgumentError("one independent variable currently supported"))
    iv = first(ivs)
    ps, = extract_elements(eqs, [_is_parameter(iv)])
    DiffEqSystem(eqs, iv, dvs, ps)
end

function DiffEqSystem(eqs, iv)
    dvs, ps = extract_elements(eqs, [_is_dependent, _is_parameter(iv)])
    DiffEqSystem(eqs, iv, dvs, ps)
end


function calculate_jacobian(sys::DiffEqSystem)
    isempty(sys.jac[]) || return sys.jac[]  # use cached Jacobian, if possible
    rhs = [eq.rhs for eq in sys.eqs]

    jac = expand_derivatives.(calculate_jacobian(rhs, sys.dvs))
    sys.jac[] = jac  # cache Jacobian
    return jac
end

function generate_jacobian(sys::DiffEqSystem; version::FunctionVersion = ArrayFunction)
    jac = calculate_jacobian(sys)
    return build_function(jac, sys.dvs, sys.ps, (sys.iv.name,); version = version)
end

function generate_function(sys::DiffEqSystem; version::FunctionVersion = ArrayFunction)
    rhss = [eq.rhs for eq ∈ sys.eqs]
    return build_function(rhss, sys.dvs, sys.ps, (sys.iv.name,); version = version)
end


function generate_factorized_W(sys::DiffEqSystem, simplify=true; version::FunctionVersion = ArrayFunction)
    jac = calculate_jacobian(sys)

    gam = Variable(:gam; known = true)

    W = LinearAlgebra.I - gam*jac
    Wfact = lu(W, Val(false), check=false).factors

    if simplify
        Wfact = simplify_constants.(Wfact)
    end

    W_t = LinearAlgebra.I/gam - jac
    Wfact_t = lu(W_t, Val(false), check=false).factors
    if simplify
        Wfact_t = simplify_constants.(Wfact_t)
    end

    vs, ps = sys.dvs, sys.ps
    Wfact_func   = build_function(Wfact  , vs, ps, (:gam,:t); version = version)
    Wfact_t_func = build_function(Wfact_t, vs, ps, (:gam,:t); version = version)

    return (Wfact_func, Wfact_t_func)
end

function DiffEqBase.ODEFunction(sys::DiffEqSystem; version::FunctionVersion = ArrayFunction)
    expr = generate_function(sys; version = version)
    if version === ArrayFunction
        ODEFunction{true}(eval(expr))
    elseif version === SArrayFunction
        ODEFunction{false}(eval(expr))
    end
end
