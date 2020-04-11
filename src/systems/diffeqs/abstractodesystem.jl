function calculate_tgrad(sys::AbstractODESystem)
  isempty(sys.tgrad[]) || return sys.tgrad[]  # use cached tgrad, if possible
  rhs = [detime_dvs(eq.rhs) for eq ∈ equations(sys)]
  iv = sys.iv()
  notime_tgrad = [expand_derivatives(ModelingToolkit.Differential(iv)(r)) for r in rhs]
  tgrad = retime_dvs.(notime_tgrad,(states(sys),),iv)
  sys.tgrad[] = tgrad
  return tgrad
end

function calculate_jacobian(sys::AbstractODESystem)
    isempty(sys.jac[]) || return sys.jac[]  # use cached Jacobian, if possible
    rhs = [eq.rhs for eq ∈ equations(sys)]

    iv = sys.iv()
    dvs = [dv(iv) for dv ∈ states(sys)]

    jac = expand_derivatives.(calculate_jacobian(rhs, dvs))
    sys.jac[] = jac  # cache Jacobian
    return jac
end

struct ODEToExpr
    sys::AbstractODESystem
end
function (f::ODEToExpr)(O::Operation)
    if isa(O.op, Variable)
        isequal(O.op, f.sys.iv) && return O.op.name  # independent variable
        O.op ∈ f.sys.dvs        && return O.op.name  # dependent variables
        isempty(O.args)         && return O.op.name  # 0-ary parameters
        return build_expr(:call, Any[O.op.name; f.(O.args)])
    end
    return build_expr(:call, Any[Symbol(O.op); f.(O.args)])
end
(f::ODEToExpr)(x) = convert(Expr, x)

function generate_tgrad(sys::AbstractODESystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)
    tgrad = calculate_tgrad(sys)
    return build_function(tgrad, dvs, ps, (sys.iv.name,), ODEToExpr(sys), expression; kwargs...)
end

function generate_jacobian(sys::AbstractODESystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)
    jac = calculate_jacobian(sys)
    return build_function(jac, dvs, ps, (sys.iv.name,), ODEToExpr(sys), expression; kwargs...)
end

function generate_function(sys::AbstractODESystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)
    rhss = [deq.rhs for deq ∈ equations(sys)]
    dvs′ = [clean(dv) for dv ∈ dvs]
    ps′ = [clean(p) for p ∈ ps]
    return build_function(rhss, dvs′, ps′, (sys.iv.name,), ODEToExpr(sys), expression; kwargs...)
end

function calculate_factorized_W(sys::AbstractODESystem, simplify=true)
    isempty(sys.Wfact[]) || return (sys.Wfact[],sys.Wfact_t[])

    jac = calculate_jacobian(sys)
    gam = Variable(gensym(:gamma))()

    W = - LinearAlgebra.I + gam*jac
    Wfact = lu(W, Val(false), check=false).factors

    if simplify
        Wfact = simplify_constants.(Wfact)
    end

    W_t = - LinearAlgebra.I/gam + jac
    Wfact_t = lu(W_t, Val(false), check=false).factors
    if simplify
        Wfact_t = simplify_constants.(Wfact_t)
    end
    sys.Wfact[] = Wfact
    sys.Wfact_t[] = Wfact_t

    (Wfact,Wfact_t)
end

function generate_factorized_W(sys::AbstractODESystem, vs = states(sys), ps = parameters(sys), simplify=true, expression = Val{true}; kwargs...)
    (Wfact,Wfact_t) = calculate_factorized_W(sys,simplify)
    siz = size(Wfact)
    constructor = :(x -> begin
                        A = SMatrix{$siz...}(x)
                        StaticArrays.LU(LowerTriangular( SMatrix{$siz...}(UnitLowerTriangular(A)) ), UpperTriangular(A), SVector(ntuple(n->n, max($siz...))))
                    end)

    Wfact_func   = build_function(Wfact  , vs, ps, (:gam,:t), ODEToExpr(sys), expression;constructor=constructor,kwargs...)
    Wfact_t_func = build_function(Wfact_t, vs, ps, (:gam,:t), ODEToExpr(sys), expression;constructor=constructor,kwargs...)

    return (Wfact_func, Wfact_t_func)
end

function calculate_massmatrix(sys::AbstractODESystem, simplify=true)
    eqs = equations(sys)
    dvs = states(sys)
    M = zeros(length(eqs),length(eqs))
    for (i,eq) in enumerate(eqs)
        if eq.lhs isa Constant
            @assert eq.lhs.value == 0
        elseif eq.lhs.op isa Differential
            j = findfirst(x->isequal(x.name,var_from_nested_derivative(eq.lhs)[1].name),dvs)
            M[i,j] = 1
        else
            error("Only semi-explicit mass matrices are currently supported")
        end
    end
    M = simplify ? simplify_constants.(M) : M
    M == I ? I : M
end

"""
$(SIGNATURES)

Create an `ODEFunction` from the [`ODESystem`](@ref). The arguments `dvs` and `ps`
are used to set the order of the dependent variable and parameter vectors,
respectively.
"""
function DiffEqBase.ODEFunction{iip}(sys::AbstractODESystem, dvs = states(sys),
                                     ps = parameters(sys);
                                     version = nothing, tgrad=false,
                                     jac = false, Wfact = false) where {iip}
    f_oop,f_iip = generate_function(sys, dvs, ps, Val{false})

    f(u,p,t) = f_oop(u,p,t)
    f(du,u,p,t) = f_iip(du,u,p,t)

    if tgrad
        tgrad_oop,tgrad_iip = generate_tgrad(sys, dvs, ps, Val{false})
        _tgrad(u,p,t) = tgrad_oop(u,p,t)
        _tgrad(J,u,p,t) = tgrad_iip(J,u,p,t)
    else
        _tgrad = nothing
    end

    if jac
        jac_oop,jac_iip = generate_jacobian(sys, dvs, ps, Val{false})
        _jac(u,p,t) = jac_oop(u,p,t)
        _jac(J,u,p,t) = jac_iip(J,u,p,t)
    else
        _jac = nothing
    end

    if Wfact
        tmp_Wfact,tmp_Wfact_t = generate_factorized_W(sys, dvs, ps, true, Val{false})
        Wfact_oop, Wfact_iip = tmp_Wfact
        Wfact_oop_t, Wfact_iip_t = tmp_Wfact_t
        _Wfact(u,p,dtgamma,t) = Wfact_oop(u,p,dtgamma,t)
        _Wfact(W,u,p,dtgamma,t) = Wfact_iip(W,u,p,dtgamma,t)
        _Wfact_t(u,p,dtgamma,t) = Wfact_oop_t(u,p,dtgamma,t)
        _Wfact_t(W,u,p,dtgamma,t) = Wfact_iip_t(W,u,p,dtgamma,t)
    else
        _Wfact,_Wfact_t = nothing,nothing
    end

    M = calculate_massmatrix(sys)

    ODEFunction{iip}(f,jac=_jac,
                      tgrad = _tgrad,
                      Wfact = _Wfact,
                      Wfact_t = _Wfact_t,
                      mass_matrix = M,
                      syms = Symbol.(sys.dvs))
end

renamespace(namespace,name) = Symbol(string(namespace)*"′"*string(name))

function DiffEqBase.ODEFunction(sys::AbstractODESystem, args...; kwargs...)
    ODEFunction{true}(sys, args...; kwargs...)
end

function namespace_variables(sys::AbstractODESystem)
    [rename(x,renamespace(sys.name,x.name)) for x in states(sys)]
end

function namespace_parameters(sys::AbstractODESystem)
    [rename(x,renamespace(sys.name,x.name)) for x in parameters(sys)]
end

namespace_equations(sys::AbstractODESystem) = namespace_equation.(equations(sys),sys.name,sys.iv.name)

function namespace_equation(eq::Equation,name,ivname)
    _lhs = namespace_operation(eq.lhs,name,ivname)
    _rhs = namespace_operation(eq.rhs,name,ivname)
    _lhs ~ _rhs
end

function namespace_operation(O::Operation,name,ivname)
    if O.op isa Variable && O.op.name != ivname
        Operation(rename(O.op,renamespace(name,O.op.name)),namespace_operation.(O.args,name,ivname))
    else
        Operation(O.op,namespace_operation.(O.args,name,ivname))
    end
end
namespace_operation(O::Constant,name,ivname) = O

independent_variable(sys::AbstractODESystem) = sys.iv
states(sys::AbstractODESystem) = isempty(sys.systems) ? sys.dvs : [sys.dvs;reduce(vcat,namespace_variables.(sys.systems))]
parameters(sys::AbstractODESystem) = isempty(sys.systems) ? sys.ps : [sys.ps;reduce(vcat,namespace_parameters.(sys.systems))]

function equations(sys::AbstractODESystem)
    isempty(sys.systems) ? sys.eqs : [sys.eqs;reduce(vcat,namespace_equations.(sys.systems))]
end

function states(sys::AbstractODESystem,name::Symbol)
    x = sys.dvs[findfirst(x->x.name==name,sys.dvs)]
    Variable(Symbol(string(sys.name)*"′"*string(x.name)),known=x.known)(sys.iv())
end

function parameters(sys::AbstractODESystem,name::Symbol)
    x = sys.ps[findfirst(x->x.name==name,sys.ps)]
    Variable(Symbol(string(sys.name)*"′"*string(x.name)),known=x.known)(sys.iv())
end

function states(sys::AbstractODESystem,args...)
    name = last(args)
    extra_names = reduce(*,["′$(x.name)" for x in args[1:end-1]])
    Variable(Symbol(string(sys.name)*extra_names*"′"*string(name)))(sys.iv())
end

function parameters(sys::AbstractODESystem,args...)
    name = last(args)
    extra_names = reduce(*,["′$(x.name)" for x in args[1:end-1]])
    Variable(Symbol(string(sys.name)*extra_names*"′"*string(name)))(sys.iv())
end

function _eq_unordered(a, b)
    length(a) === length(b) || return false
    n = length(a)
    idxs = Set(1:n)
    for x ∈ a
        idx = findfirst(isequal(x), b)
        idx === nothing && return false
        idx ∈ idxs      || return false
        delete!(idxs, idx)
    end
    return true
end
