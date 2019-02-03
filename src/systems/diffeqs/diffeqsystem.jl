export DiffEqSystem, ODEFunction


using Base: RefValue


isintermediate(eq::Equation) = !(isa(eq.lhs, Operation) && isa(eq.lhs.op, Differential))

struct DiffEq  # D(x) = t
    D::Differential  # D
    var::Variable    # x
    rhs::Expression  # t
end
function Base.convert(::Type{DiffEq}, eq::Equation)
    isintermediate(eq) && throw(ArgumentError("intermediate equation received"))
    return DiffEq(eq.lhs.op, eq.lhs.args[1], eq.rhs)
end
Base.:(==)(a::DiffEq, b::DiffEq) = (a.D, a.var, a.rhs) == (b.D, b.var, b.rhs)
get_args(eq::DiffEq) = Expression[eq.var, eq.rhs]

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
