using Base: RefValue


isintermediate(eq::Equation) = !(isa(eq.lhs, Operation) && isa(eq.lhs.op, Differential))

mutable struct DiffEq  # D(x) = t
    D::Differential  # D
    var::Variable    # x
    rhs::Expression  # t
end
function Base.convert(::Type{DiffEq}, eq::Equation)
    isintermediate(eq) && throw(ArgumentError("intermediate equation received"))
    return DiffEq(eq.lhs.op, eq.lhs.args[1], eq.rhs)
end
get_args(eq::DiffEq) = Expression[eq.var, eq.rhs]

struct DiffEqSystem <: AbstractSystem
    eqs::Vector{DiffEq}
    ivs::Vector{Variable}
    dvs::Vector{Variable}
    ps::Vector{Variable}
    jac::RefValue{Matrix{Expression}}
    function DiffEqSystem(eqs, ivs, dvs, ps)
        jac = RefValue(Matrix{Expression}(undef, 0, 0))
        new(eqs, ivs, dvs, ps, jac)
    end
end

function DiffEqSystem(eqs)
    dvs, = extract_elements(eqs, [_is_dependent])
    ivs = unique(vcat((dv.dependents for dv ∈ dvs)...))
    ps, = extract_elements(eqs, [_is_parameter(ivs)])
    DiffEqSystem(eqs, ivs, dvs, ps)
end

function DiffEqSystem(eqs, ivs)
    dvs, ps = extract_elements(eqs, [_is_dependent, _is_parameter(ivs)])
    DiffEqSystem(eqs, ivs, dvs, ps)
end


function generate_ode_function(sys::DiffEqSystem; version::FunctionVersion = ArrayFunction)
    var_exprs = [:($(sys.dvs[i].name) = u[$i]) for i in eachindex(sys.dvs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in eachindex(sys.ps)]
    sys_exprs = build_equals_expr.(sys.eqs)
    if version === ArrayFunction
        dvar_exprs = [:(du[$i] = $(Symbol("$(sys.dvs[i].name)_$(sys.ivs[1].name)"))) for i in eachindex(sys.dvs)]
        exprs = vcat(var_exprs,param_exprs,sys_exprs,dvar_exprs)
        block = expr_arr_to_block(exprs)
        :((du,u,p,t)->$(toexpr(block)))
    elseif version === SArrayFunction
        dvar_exprs = [:($(Symbol("$(sys.dvs[i].name)_$(sys.ivs[1].name)"))) for i in eachindex(sys.dvs)]
        svector_expr = quote
            E = eltype(tuple($(dvar_exprs...)))
            T = StaticArrays.similar_type(typeof(u), E)
            T($(dvar_exprs...))
        end
        exprs = vcat(var_exprs,param_exprs,sys_exprs,svector_expr)
        block = expr_arr_to_block(exprs)
        :((u,p,t)->$(toexpr(block)))
    end
end

function build_equals_expr(eq::DiffEq)
    lhs = Symbol(eq.var.name, :_, eq.D.x.name)
    return :($lhs = $(convert(Expr, eq.rhs)))
end

function calculate_jacobian(sys::DiffEqSystem, simplify=true)
    isempty(sys.jac[]) || return sys.jac[]  # use cached Jacobian, if possible
    rhs = [eq.rhs for eq in sys.eqs]

    jac = expand_derivatives.(calculate_jacobian(rhs, sys.dvs))
    sys.jac[] = jac  # cache Jacobian
    return jac
end

function generate_ode_jacobian(sys::DiffEqSystem, simplify=true)
    var_exprs = [:($(sys.dvs[i].name) = u[$i]) for i in eachindex(sys.dvs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in eachindex(sys.ps)]
    jac = calculate_jacobian(sys, simplify)
    jac_exprs = [:(J[$i,$j] = $(convert(Expr, jac[i,j]))) for i in 1:size(jac,1), j in 1:size(jac,2)]
    exprs = vcat(var_exprs,param_exprs,vec(jac_exprs))
    block = expr_arr_to_block(exprs)
    :((J,u,p,t)->$(block))
end

function generate_ode_iW(sys::DiffEqSystem, simplify=true)
    var_exprs = [:($(sys.dvs[i].name) = u[$i]) for i in eachindex(sys.dvs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in eachindex(sys.ps)]
    jac = calculate_jacobian(sys, simplify)

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

    iW_exprs = [:(iW[$i,$j] = $(convert(Expr, iW[i,j]))) for i in 1:size(iW,1), j in 1:size(iW,2)]
    exprs = vcat(var_exprs,param_exprs,vec(iW_exprs))
    block = expr_arr_to_block(exprs)

    iW_t_exprs = [:(iW[$i,$j] = $(convert(Expr, iW_t[i,j]))) for i in 1:size(iW_t,1), j in 1:size(iW_t,2)]
    exprs = vcat(var_exprs,param_exprs,vec(iW_t_exprs))
    block2 = expr_arr_to_block(exprs)
    :((iW,u,p,gam,t)->$(block)),:((iW,u,p,gam,t)->$(block2))
end

function DiffEqBase.ODEFunction(sys::DiffEqSystem; version::FunctionVersion = ArrayFunction)
    expr = generate_ode_function(sys; version = version)
    if version === ArrayFunction
        ODEFunction{true}(eval(expr))
    elseif version === SArrayFunction
        ODEFunction{false}(eval(expr))
    end
end


export DiffEqSystem, ODEFunction
export generate_ode_function
