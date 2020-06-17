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
    states::Vector{Variable}
end
ODEToExpr(sys) = ODEToExpr(sys,states(sys))
function (f::ODEToExpr)(O::Operation)
    if isa(O.op, Variable)
        isequal(O.op, f.sys.iv) && return O.op.name  # independent variable
        O.op ∈ f.states         && return O.op.name  # dependent variables
        isempty(O.args)         && return O.op.name  # 0-ary parameters
        return build_expr(:call, Any[O.op.name; f.(O.args)])
    end
    return build_expr(:call, Any[Symbol(O.op); f.(O.args)])
end
(f::ODEToExpr)(x) = convert(Expr, x)

function generate_tgrad(sys::AbstractODESystem, dvs = states(sys), ps = parameters(sys); kwargs...)
    tgrad = calculate_tgrad(sys)
    return build_function(tgrad, dvs, ps, sys.iv;
                          conv = ODEToExpr(sys), kwargs...)
end

function generate_jacobian(sys::AbstractODESystem, dvs = states(sys), ps = parameters(sys); sparse = false, kwargs...)
    jac = calculate_jacobian(sys)
    if sparse
        jac = SparseArrays.sparse(jac)
    end
    return build_function(jac, dvs, ps, sys.iv;
                          conv = ODEToExpr(sys), kwargs...)
end

function generate_function(sys::AbstractODESystem, dvs = states(sys), ps = parameters(sys); kwargs...)
    rhss = [deq.rhs for deq ∈ equations(sys)]
    dvs′ = convert.(Variable,dvs)
    ps′ = convert.(Variable,ps)
    return build_function(rhss, dvs′, ps′, sys.iv;
                          conv = ODEToExpr(sys),kwargs...)
end

function calculate_factorized_W(sys::AbstractODESystem, simplify=true)
    isempty(sys.Wfact[]) || return (sys.Wfact[],sys.Wfact_t[])

    jac = calculate_jacobian(sys)
    M = calculate_massmatrix(sys)
    gam = Variable(:__MTKWgamma)()

    W = - M + gam*jac
    Wfact = lu(W, Val(false), check=false).factors

    if simplify
        Wfact = ModelingToolkit.simplify.(Wfact)
    end

    W_t = - M/gam + jac
    Wfact_t = lu(W_t, Val(false), check=false).factors
    if simplify
        Wfact_t = ModelingToolkit.simplify.(Wfact_t)
    end
    sys.Wfact[] = Wfact
    sys.Wfact_t[] = Wfact_t

    (Wfact,Wfact_t)
end

function generate_factorized_W(sys::AbstractODESystem, vs = states(sys), ps = parameters(sys), simplify=true, expression = Val{true}; kwargs...)
    Wfact,Wfact_t = calculate_factorized_W(sys,simplify)
    siz = size(Wfact)
    constructor = :(x -> begin
                        A = SMatrix{$siz...}(x)
                        StaticArrays.LU(LowerTriangular( SMatrix{$siz...}(UnitLowerTriangular(A)) ), UpperTriangular(A), SVector(ntuple(n->n, max($siz...))))
                    end)

    Wfact_func   = build_function(Wfact  , vs, ps, Variable(:__MTKWgamma), sys.iv;
                                  conv = ODEToExpr(sys), expression = expression, constructor=constructor,kwargs...)
    Wfact_t_func = build_function(Wfact_t, vs, ps, Variable(:__MTKWgamma), sys.iv;
                                  conv = ODEToExpr(sys), expression = expression, constructor=constructor,kwargs...)

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
            error("Only semi-explicit constant mass matrices are currently supported")
        end
    end
    M = simplify ? ModelingToolkit.simplify.(M) : M
    # M should only contain concrete numbers
    M = map(x->x isa Constant ? x.value : x, M)
    M == I ? I : M
end

function DiffEqBase.ODEFunction(sys::AbstractODESystem, args...; kwargs...)
    ODEFunction{true}(sys, args...; kwargs...)
end

"""
```julia
function DiffEqBase.ODEFunction{iip}(sys::AbstractODESystem, dvs = states(sys),
                                     ps = parameters(sys);
                                     version = nothing, tgrad=false,
                                     jac = false, Wfact = false,
                                     sparse = false,
                                     kwargs...) where {iip}
```

Create an `ODEFunction` from the [`ODESystem`](@ref). The arguments `dvs` and `ps`
are used to set the order of the dependent variable and parameter vectors,
respectively.
"""
function DiffEqBase.ODEFunction{iip}(sys::AbstractODESystem, dvs = states(sys),
                                     ps = parameters(sys), u0 = nothing;
                                     version = nothing, tgrad=false,
                                     jac = false, Wfact = false,
                                     sparse = false,
                                     kwargs...) where {iip}

    f_oop,f_iip = ModelingToolkit.eval.(generate_function(sys, dvs, ps; expression=Val{true}, kwargs...))
    f(u,p,t) = f_oop(u,p,t)
    f(du,u,p,t) = f_iip(du,u,p,t)

    if tgrad
        tgrad_oop,tgrad_iip = ModelingToolkit.eval.(generate_tgrad(sys, dvs, ps; expression=Val{true}, kwargs...))
        _tgrad(u,p,t) = tgrad_oop(u,p,t)
        _tgrad(J,u,p,t) = tgrad_iip(J,u,p,t)
    else
        _tgrad = nothing
    end

    if jac
        jac_oop,jac_iip = ModelingToolkit.eval.(generate_jacobian(sys, dvs, ps; sparse = sparse, expression=Val{true}, kwargs...))
        _jac(u,p,t) = jac_oop(u,p,t)
        _jac(J,u,p,t) = jac_iip(J,u,p,t)
    else
        _jac = nothing
    end

    if Wfact
        tmp_Wfact,tmp_Wfact_t = generate_factorized_W(sys, dvs, ps; expression=Val{true}, kwargs...)
        Wfact_oop, Wfact_iip =  ModelingToolkit.eval.(tmp_Wfact)
        Wfact_oop_t, Wfact_iip_t =  ModelingToolkit.eval.(tmp_Wfact_t)
        _Wfact(u,p,dtgamma,t) = Wfact_oop(u,p,dtgamma,t)
        _Wfact(W,u,p,dtgamma,t) = Wfact_iip(W,u,p,dtgamma,t)
        _Wfact_t(u,p,dtgamma,t) = Wfact_oop_t(u,p,dtgamma,t)
        _Wfact_t(W,u,p,dtgamma,t) = Wfact_iip_t(W,u,p,dtgamma,t)
    else
        _Wfact,_Wfact_t = nothing,nothing
    end

    M = calculate_massmatrix(sys)

    _M = (u0 === nothing || M == I) ? M : ArrayInterface.restructure(u0 .* u0',M)

    ODEFunction{iip}(DiffEqBase.EvalFunc(f),
                      jac = _jac === nothing ? nothing : DiffEqBase.EvalFunc(_jac),
                      tgrad = _tgrad === nothing ? nothing : DiffEqBase.EvalFunc(_tgrad),
                      Wfact = _Wfact === nothing ? nothing : DiffEqBase.EvalFunc(_Wfact),
                      Wfact_t = _Wfact_t === nothing ? nothing : DiffEqBase.EvalFunc(_Wfact_t),
                      mass_matrix = _M,
                      syms = Symbol.(states(sys)))
end

"""
```julia
function DiffEqBase.ODEFunctionExpr{iip}(sys::AbstractODESystem, dvs = states(sys),
                                     ps = parameters(sys);
                                     version = nothing, tgrad=false,
                                     jac = false, Wfact = false,
                                     sparse = false,
                                     kwargs...) where {iip}
```

Create a Julia expression for an `ODEFunction` from the [`ODESystem`](@ref).
The arguments `dvs` and `ps` are used to set the order of the dependent
variable and parameter vectors, respectively.
"""
struct ODEFunctionExpr{iip} end

function ODEFunctionExpr{iip}(sys::AbstractODESystem, dvs = states(sys),
                                     ps = parameters(sys), u0 = nothing;
                                     version = nothing, tgrad=false,
                                     jac = false, Wfact = false,
                                     sparse = false,linenumbers = false,
                                     kwargs...) where {iip}

    idx = iip ? 2 : 1
    f = generate_function(sys, dvs, ps; expression=Val{true}, kwargs...)[idx]
    if tgrad
        _tgrad = generate_tgrad(sys, dvs, ps; expression=Val{true}, kwargs...)[idx]
    else
        _tgrad = :nothing
    end

    if jac
        _jac = generate_jacobian(sys, dvs, ps; sparse = sparse, expression=Val{true}, kwargs...)[idx]
    else
        _jac = :nothing
    end

    if Wfact
        tmp_Wfact,tmp_Wfact_t = generate_factorized_W(sys, dvs, ps; expression=Val{true}, kwargs...)
        _Wfact =  tmp_Wfact[idx]
        _Wfact_t =  tmp_Wfact_t[idx]
    else
        _Wfact,_Wfact_t = :nothing,:nothing
    end

    M = calculate_massmatrix(sys)

    _M = (u0 === nothing || M == I) ? M : ArrayInterface.restructure(u0 .* u0',M)

    ex = quote
        f = $f
        tgrad = $_tgrad
        jac = $_jac
        Wfact = $_Wfact
        Wfact_t = $_Wfact_t
        M = $_M

        ODEFunction{iip}(f,
                         jac = jac,
                         tgrad = tgrad,
                         Wfact = Wfact,
                         Wfact_t = Wfact_t,
                         mass_matrix = M,
                         syms = $(Symbol.(states(sys))))
    end
    !linenumbers ? striplines(ex) : ex
end


function ODEFunctionExpr(sys::AbstractODESystem, args...; kwargs...)
    ODEFunctionExpr{true}(sys, args...; kwargs...)
end


function DiffEqBase.ODEProblem(sys::AbstractODESystem, args...; kwargs...)
    ODEProblem{true}(sys, args...; kwargs...)
end

"""
```julia
function DiffEqBase.ODEProblem{iip}(sys::AbstractODESystem,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false, Wfact = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = true, parallel=SerialForm(),
                                    kwargs...) where iip
```

Generates an ODEProblem from an ODESystem and allows for automatically
symbolically calculating numerical enhancements.
"""
function DiffEqBase.ODEProblem{iip}(sys::AbstractODESystem,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false, Wfact = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = true, parallel=SerialForm(),
                                    kwargs...) where iip
    dvs = states(sys)
    ps = parameters(sys)
    u0 = varmap_to_vars(u0map,dvs)
    p = varmap_to_vars(parammap,ps)
    f = ODEFunction{iip}(sys,dvs,ps,u0;tgrad=tgrad,jac=jac,Wfact=Wfact,checkbounds=checkbounds,
                        linenumbers=linenumbers,parallel=parallel,
                        sparse=sparse)
    ODEProblem{iip}(f,u0,tspan,p;kwargs...)
end

"""
```julia
function DiffEqBase.ODEProblemExpr{iip}(sys::AbstractODESystem,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false, Wfact = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = true, parallel=SerialForm(),
                                    kwargs...) where iip
```

Generates a Julia expression for constructing an ODEProblem from an
ODESystem and allows for automatically symbolically calculating
numerical enhancements.
"""
struct ODEProblemExpr{iip} end

function ODEProblemExpr{iip}(sys::AbstractODESystem,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false, Wfact = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = false, parallel=SerialForm(),
                                    kwargs...) where iip
    dvs = states(sys)
    ps = parameters(sys)
    u0 = varmap_to_vars(u0map,dvs)
    p = varmap_to_vars(parammap,ps)
    f = ODEFunctionExpr{iip}(sys,dvs,ps,u0;tgrad=tgrad,jac=jac,Wfact=Wfact,checkbounds=checkbounds,
                        linenumbers=linenumbers,parallel=parallel,
                        sparse=sparse)
    ex = quote
        f = $f
        u0 = $u0
        tspan = $tspan
        p = $p
        ODEProblem(f,u0,tspan,p;$(kwargs...))
    end
    !linenumbers ? striplines(ex) : ex
end

function ODEProblemExpr(sys::AbstractODESystem, args...; kwargs...)
    ODEProblemExpr{true}(sys, args...; kwargs...)
end


### Enables Steady State Problems ###
function DiffEqBase.SteadyStateProblem(sys::AbstractODESystem, args...; kwargs...)
    SteadyStateProblem{true}(sys, args...; kwargs...)
end

"""
```julia
function DiffEqBase.SteadyStateProblem(sys::AbstractODESystem,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false, Wfact = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = true, parallel=SerialForm(),
                                    kwargs...) where iip
```
Generates an SteadyStateProblem from an ODESystem and allows for automatically
symbolically calculating numerical enhancements.
"""
function DiffEqBase.SteadyStateProblem{iip}(sys::AbstractODESystem,u0map,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false, Wfact = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = true, parallel=SerialForm(),
                                    kwargs...) where iip
    dvs = states(sys)
    ps = parameters(sys)
    u0 = varmap_to_vars(u0map,dvs)
    p = varmap_to_vars(parammap,ps)
    f = ODEFunction(sys,dvs,ps,u0;tgrad=tgrad,jac=jac,Wfact=Wfact,checkbounds=checkbounds,
                        linenumbers=linenumbers,parallel=parallel,
                        sparse=sparse)
    SteadyStateProblem(f,u0,p;kwargs...)
end

struct SteadyStateProblemExpr{iip} end

"""
```julia
function DiffEqBase.SteadyStateProblem(sys::AbstractODESystem,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false, Wfact = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = true, parallel=SerialForm(),
                                    kwargs...) where iip
```
Generates an SteadyStateProblem from an ODESystem and allows for automatically
symbolically calculating numerical enhancements.
"""
function DiffEqBase.SteadyStateProblemExpr{iip}(sys::AbstractODESystem,u0map,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false, Wfact = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = true, parallel=SerialForm(),
                                    kwargs...) where iip
    dvs = states(sys)
    ps = parameters(sys)
    u0 = varmap_to_vars(u0map,dvs)
    p = varmap_to_vars(parammap,ps)
    f = ODEFunctionExpr(sys,dvs,ps,u0;tgrad=tgrad,jac=jac,Wfact=Wfact,checkbounds=checkbounds,
                        linenumbers=linenumbers,parallel=parallel,
                        sparse=sparse)
    ex = quote
        f = $f
        u0 = $u0
        tspan = $tspan
        p = $p
        SteadyStateProblem(f,u0,tspan,p;$(kwargs...))
    end
    !linenumbers ? striplines(ex) : ex
end

function DiffEqBase.SteadyStateProblemExpr(sys::AbstractODESystem, args...; kwargs...)
    SteadyStateProblemExpr{true}(sys, args...; kwargs...)
end
