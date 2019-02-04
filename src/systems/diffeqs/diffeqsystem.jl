export DiffEqSystem, ODEFunction


using Base: RefValue


isintermediate(eq::Equation) = !(isa(eq.lhs, Operation) && isa(eq.lhs.op, Differential))

function _unwrap_differenital(O)
    isa(O, Operation)       || return (O, nothing, 0)
    isa(O.op, Differential) || return (O, nothing, 0)
    (x, t, order) = _unwrap_differenital(O.args[1])
    t === nothing && (t = O.op.x)
    t == O.op.x || throw(ArgumentError("non-matching differentials on lhs"))
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
    (x, t, n) = _unwrap_differenital(eq.lhs)
    return DiffEq(x, t, n, eq.rhs)
end
Base.:(==)(a::DiffEq, b::DiffEq) = (a.x, a.t, a.n, a.rhs) == (b.x, b.t, b.n, b.rhs)
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


function generate_ode_iW(sys::DiffEqSystem, simplify=true; version::FunctionVersion = ArrayFunction)
    jac = calculate_jacobian(sys)

    gam = Parameter(:gam)

    W = LinearAlgebra.I - gam*jac
    W = SMatrix{size(W,1),size(W,2)}(W)
    iW = inv(W)

    if simplify
        iW = simplify_constants.(iW)
    end

    W = inv(LinearAlgebra.I/gam - jac)
    W = SMatrix{size(W,1),size(W,2)}(W)
    iW_t = inv(W)
    if simplify
        iW_t = simplify_constants.(iW_t)
    end

    vs, ps = sys.dvs, sys.ps
    iW_func   = build_function(iW  , vs, ps, (:gam,:t); version = version)
    iW_t_func = build_function(iW_t, vs, ps, (:gam,:t); version = version)

    return (iW_func, iW_t_func)
end

function DiffEqBase.ODEFunction(sys::DiffEqSystem; version::FunctionVersion = ArrayFunction)
    expr = generate_function(sys; version = version)
    if version === ArrayFunction
        ODEFunction{true}(eval(expr))
    elseif version === SArrayFunction
        ODEFunction{false}(eval(expr))
    end
end
