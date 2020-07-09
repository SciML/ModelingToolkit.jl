function calculate_tgrad(sys::AbstractODESystem;
                         simplify=true)
  isempty(sys.tgrad[]) || return sys.tgrad[]  # use cached tgrad, if possible
  rhs = [detime_dvs(eq.rhs) for eq ∈ equations(sys)]
  iv = sys.iv()
  notime_tgrad = [expand_derivatives(ModelingToolkit.Differential(iv)(r)) for r in rhs]
  tgrad = retime_dvs.(notime_tgrad,(states(sys),),iv)
  if simplify
      tgrad = ModelingToolkit.simplify.(tgrad)
  end
  sys.tgrad[] = tgrad
  return tgrad
end

function calculate_jacobian(sys::AbstractODESystem;
                            sparse=false, simplify=true)
    isempty(sys.jac[]) || return sys.jac[]  # use cached Jacobian, if possible
    rhs = [eq.rhs for eq ∈ equations(sys)]

    iv = sys.iv()
    dvs = [dv(iv) for dv ∈ states(sys)]

    if sparse
        jac = sparsejacobian(rhs, dvs, simplify=simplify)
    else
        jac = jacobian(rhs, dvs, simplify=simplify)
    end

    sys.jac[] = jac  # cache Jacobian
    return jac
end

struct ODEToExpr
    sys::AbstractODESystem
    states::Vector{Variable}
end
ODEToExpr(@nospecialize(sys)) = ODEToExpr(sys,states(sys))
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

function generate_tgrad(sys::AbstractODESystem, dvs = states(sys), ps = parameters(sys);
                        simplify = true, kwargs...)
    tgrad = calculate_tgrad(sys,simplify=simplify)
    return build_function(tgrad, dvs, ps, sys.iv;
                          conv = ODEToExpr(sys), kwargs...)
end

function generate_jacobian(sys::AbstractODESystem, dvs = states(sys), ps = parameters(sys);
                           simplify = true, sparse = false, kwargs...)
    jac = calculate_jacobian(sys;simplify=simplify,sparse=sparse)
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

function calculate_massmatrix(sys::AbstractODESystem; simplify=true)
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
                                     jac = false,
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
                                     jac = false,
                                     eval_expression = true,
                                     sparse = false, simplify = true,
                                     kwargs...) where {iip}

    f_gen = generate_function(sys, dvs, ps; expression=Val{eval_expression}, kwargs...)
    f_oop,f_iip = eval_expression ? ModelingToolkit.eval.(f_gen) : f_gen
    f(u,p,t) = f_oop(u,p,t)
    f(du,u,p,t) = f_iip(du,u,p,t)

    if tgrad
        tgrad_gen = generate_tgrad(sys, dvs, ps;
                                   simplify=simplify,
                                   expression=Val{eval_expression}, kwargs...)
        tgrad_oop,tgrad_iip = eval_expression ? ModelingToolkit.eval.(tgrad_gen) : tgrad_gen
        _tgrad(u,p,t) = tgrad_oop(u,p,t)
        _tgrad(J,u,p,t) = tgrad_iip(J,u,p,t)
    else
        _tgrad = nothing
    end

    if jac
        jac_gen = generate_jacobian(sys, dvs, ps;
                                    simplify=simplify, sparse = sparse,
                                    expression=Val{eval_expression}, kwargs...)
        jac_oop,jac_iip = eval_expression ? ModelingToolkit.eval.(jac_gen) : jac_gen
        _jac(u,p,t) = jac_oop(u,p,t)
        _jac(J,u,p,t) = jac_iip(J,u,p,t)
    else
        _jac = nothing
    end

    M = calculate_massmatrix(sys)

    _M = (u0 === nothing || M == I) ? M : ArrayInterface.restructure(u0 .* u0',M)

    ODEFunction{iip}(DiffEqBase.EvalFunc(f),
                      jac = _jac === nothing ? nothing : DiffEqBase.EvalFunc(_jac),
                      tgrad = _tgrad === nothing ? nothing : DiffEqBase.EvalFunc(_tgrad),
                      mass_matrix = _M,
                      jac_prototype = sparse ? similar(sys.jac[],Float64) : nothing,
                      syms = Symbol.(states(sys)))
end

"""
```julia
function DiffEqBase.ODEFunctionExpr{iip}(sys::AbstractODESystem, dvs = states(sys),
                                     ps = parameters(sys);
                                     version = nothing, tgrad=false,
                                     jac = false,
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
                                     jac = false,
                                     linenumbers = false,
                                     sparse = false, simplify = true,
                                     kwargs...) where {iip}

    idx = iip ? 2 : 1
    f = generate_function(sys, dvs, ps; expression=Val{true}, kwargs...)[idx]
    if tgrad
        _tgrad = generate_tgrad(sys, dvs, ps;
                                simplify=simplify,
                                expression=Val{true}, kwargs...)[idx]
    else
        _tgrad = :nothing
    end

    if jac
        _jac = generate_jacobian(sys, dvs, ps;
                                 sparse=sparse, simplify=simplify,
                                 expression=Val{true}, kwargs...)[idx]
    else
        _jac = :nothing
    end

    M = calculate_massmatrix(sys)

    _M = (u0 === nothing || M == I) ? M : ArrayInterface.restructure(u0 .* u0',M)

    jp_expr = sparse ? :(similar($(sys.jac[]),Float64)) : :nothing

    ex = quote
        f = $f
        tgrad = $_tgrad
        jac = $_jac
        M = $_M

        ODEFunction{$iip}(f,
                         jac = jac,
                         tgrad = tgrad,
                         mass_matrix = M,
                         jac_prototype = $jp_expr,
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
                                    jac = false,
                                    checkbounds = false, sparse = false,
                                    simplify = true,
                                    linenumbers = true, parallel=SerialForm(),
                                    kwargs...) where iip
```

Generates an ODEProblem from an ODESystem and allows for automatically
symbolically calculating numerical enhancements.
"""
function DiffEqBase.ODEProblem{iip}(sys::AbstractODESystem,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false,
                                    checkbounds = false, sparse = false,
                                    simplify = true,
                                    linenumbers = true, parallel=SerialForm(),
                                    eval_expression = true,
                                    kwargs...) where iip
    dvs = states(sys)
    ps = parameters(sys)
    u0 = varmap_to_vars(u0map,dvs)
    p = varmap_to_vars(parammap,ps)
    f = ODEFunction{iip}(sys,dvs,ps,u0;tgrad=tgrad,jac=jac,checkbounds=checkbounds,
                        linenumbers=linenumbers,parallel=parallel,simplify=simplify,
                        sparse=sparse,eval_expression=eval_expression,kwargs...)
    ODEProblem{iip}(f,u0,tspan,p;kwargs...)
end

"""
```julia
function DiffEqBase.ODEProblemExpr{iip}(sys::AbstractODESystem,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = true, parallel=SerialForm(),
                                    skipzeros=true, fillzeros=true,
                                    simplify = true,
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
                                    jac = false,
                                    checkbounds = false, sparse = false,
                                    simplify = true,
                                    linenumbers = false, parallel=SerialForm(),
                                    kwargs...) where iip
    dvs = states(sys)
    ps = parameters(sys)
    u0 = varmap_to_vars(u0map,dvs)
    p = varmap_to_vars(parammap,ps)
    f = ODEFunctionExpr{iip}(sys,dvs,ps,u0;tgrad=tgrad,jac=jac,checkbounds=checkbounds,
                        linenumbers=linenumbers,parallel=parallel,
                        simplify=simplify,
                        sparse=sparse,kwargs...)
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
                                    jac = false,
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
                                    jac = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = true, parallel=SerialForm(),
                                    kwargs...) where iip
    dvs = states(sys)
    ps = parameters(sys)
    u0 = varmap_to_vars(u0map,dvs)
    p = varmap_to_vars(parammap,ps)
    f = ODEFunction(sys,dvs,ps,u0;tgrad=tgrad,jac=jac,checkbounds=checkbounds,
                        linenumbers=linenumbers,parallel=parallel,
                        sparse=sparse,kwargs...)
    SteadyStateProblem(f,u0,p;kwargs...)
end

"""
```julia
function DiffEqBase.SteadyStateProblemExpr(sys::AbstractODESystem,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false,
                                    checkbounds = false, sparse = false,
                                    skipzeros=true, fillzeros=true,
                                    linenumbers = true, parallel=SerialForm(),
                                    kwargs...) where iip
```
Generates a Julia expression for building a SteadyStateProblem from
an ODESystem and allows for automatically symbolically calculating
numerical enhancements.
"""
struct SteadyStateProblemExpr{iip} end

function SteadyStateProblemExpr{iip}(sys::AbstractODESystem,u0map,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = true, parallel=SerialForm(),
                                    kwargs...) where iip
    dvs = states(sys)
    ps = parameters(sys)
    u0 = varmap_to_vars(u0map,dvs)
    p = varmap_to_vars(parammap,ps)
    f = ODEFunctionExpr(sys,dvs,ps,u0;tgrad=tgrad,jac=jac,checkbounds=checkbounds,
                        linenumbers=linenumbers,parallel=parallel,
                        sparse=sparse,kwargs...)
    ex = quote
        f = $f
        u0 = $u0
        p = $p
        SteadyStateProblem(f,u0,p;$(kwargs...))
    end
    !linenumbers ? striplines(ex) : ex
end

function SteadyStateProblemExpr(sys::AbstractODESystem, args...; kwargs...)
    SteadyStateProblemExpr{true}(sys, args...; kwargs...)
end
