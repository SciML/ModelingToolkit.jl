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
    x::Variable
    n::Int
    rhs::Expression
end
function to_diffeq(eq::Equation)
    isintermediate(eq) && throw(ArgumentError("intermediate equation received"))
    (x, t, n) = flatten_differential(eq.lhs)
    (isa(t, Operation) && isa(t.op, Variable) && isempty(t.args)) ||
        throw(ArgumentError("invalid independent variable $t"))
    (isa(x, Operation) && isa(x.op, Variable) && length(x.args) == 1 && isequal(first(x.args), t)) ||
        throw(ArgumentError("invalid dependent variable $x"))
    return t.op, DiffEq(x.op, n, eq.rhs)
end
Base.:(==)(a::DiffEq, b::DiffEq) = isequal((a.x, a.t, a.n, a.rhs), (b.x, b.t, b.n, b.rhs))

struct DiffEqSystem <: AbstractSystem
    eqs::Vector{DiffEq}
    iv::Variable
    dvs::Vector{Variable}
    ps::Vector{Variable}
    jac::RefValue{Matrix{Expression}}
    function DiffEqSystem(eqs)
        reformatted = to_diffeq.(eqs)

        ivs = unique(r[1] for r ∈ reformatted)
        length(ivs) == 1 || throw(ArgumentError("one independent variable currently supported"))
        iv = first(ivs)

        deqs = [r[2] for r ∈ reformatted]

        dvs = [deq.x for deq ∈ deqs]
        ps = filter(vars(deq.rhs for deq ∈ deqs)) do x
            x.known & !isequal(x, iv)
        end |> collect

        jac = RefValue(Matrix{Expression}(undef, 0, 0))

        new(deqs, iv, dvs, ps, jac)
    end
end


function calculate_jacobian(sys::DiffEqSystem)
    isempty(sys.jac[]) || return sys.jac[]  # use cached Jacobian, if possible
    rhs = [eq.rhs for eq ∈ sys.eqs]

    jac = expand_derivatives.(calculate_jacobian(rhs, sys.dvs, sys.iv))
    sys.jac[] = jac  # cache Jacobian
    return jac
end

function generate_jacobian(sys::DiffEqSystem; version::FunctionVersion = ArrayFunction)
    jac = calculate_jacobian(sys)
    return build_function(jac, sys.dvs, sys.ps, (sys.iv.name,); version = version)
end

struct DiffEqToExpr
    sys::DiffEqSystem
end
function (f::DiffEqToExpr)(O::Operation)
    if isa(O.op, Variable)
        isequal(O.op, f.sys.iv) && return O.op.name  # independent variable
        O.op ∈ f.sys.dvs        && return O.op.name  # dependent variables
        isempty(O.args)         && return O.op.name  # 0-ary parameters
        return build_expr(:call, Any[O.op.name; f.(O.args)])
    end
    return build_expr(:call, Any[O.op; f.(O.args)])
end
(f::DiffEqToExpr)(x) = convert(Expr, x)

function generate_function(sys::DiffEqSystem, vs, ps; version::FunctionVersion = ArrayFunction)
    rhss = [deq.rhs for deq ∈ sys.eqs]
    vs′ = [clean(v) for v ∈ vs]
    ps′ = [clean(p) for p ∈ ps]
    return build_function(rhss, vs′, ps′, (sys.iv.name,), DiffEqToExpr(sys); version = version)
end


function generate_ode_iW(sys::DiffEqSystem, simplify=true; version::FunctionVersion = ArrayFunction)
    jac = calculate_jacobian(sys)

    gam = Variable(:gam; known = true)()

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
