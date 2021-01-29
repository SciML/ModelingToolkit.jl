"""
```julia
calculate_tgrad(sys::AbstractSystem)
```

Calculate the time gradient of a system.

Returns a vector of [`Num`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_tgrad end

"""
```julia
calculate_gradient(sys::AbstractSystem)
```

Calculate the gradient of a scalar system.

Returns a vector of [`Num`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_gradient end

"""
```julia
calculate_jacobian(sys::AbstractSystem)
```

Calculate the jacobian matrix of a system.

Returns a matrix of [`Num`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_jacobian end

"""
```julia
calculate_factorized_W(sys::AbstractSystem)
```

Calculate the factorized W-matrix of a system.

Returns a matrix of [`Num`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_factorized_W end

"""
```julia
calculate_hessian(sys::AbstractSystem)
```

Calculate the hessian matrix of a scalar system.

Returns a matrix of [`Num`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_hessian end

"""
```julia
generate_tgrad(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)
```

Generates a function for the time gradient of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_tgrad end

"""
```julia
generate_gradient(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)
```

Generates a function for the gradient of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_gradient end

"""
```julia
generate_jacobian(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; sparse = false, kwargs...)
```

Generates a function for the jacobian matrix matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_jacobian end

"""
```julia
generate_factorized_W(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; sparse = false, kwargs...)
```

Generates a function for the factorized W-matrix matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_factorized_W end

"""
```julia
generate_hessian(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; sparse = false, kwargs...)
```

Generates a function for the hessian matrix matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_hessian end

"""
```julia
generate_function(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)
```

Generate a function to evaluate the system's equations.
"""
function generate_function end

function getname(t)
    if istree(t)
        operation(t) isa Sym ? getname(operation(t)) : error("Cannot get name of $t")
    else
        nameof(t)
    end
end

function Base.getproperty(sys::AbstractSystem, name::Symbol)

    if name ∈ fieldnames(typeof(sys))
        return getfield(sys,name)
    elseif !isempty(sys.systems)
        i = findfirst(x->x.name==name,sys.systems)
        if i !== nothing
            return rename(sys.systems[i],renamespace(sys.name,name))
        end
    end

    i = findfirst(x->getname(x) == name, sys.states)

    if i !== nothing
        return rename(sys.states[i],renamespace(sys.name,name))
    end

    if :ps ∈ fieldnames(typeof(sys))
        i = findfirst(x->getname(x) == name,sys.ps)
        if i !== nothing
            return rename(sys.ps[i],renamespace(sys.name,name))
        end
    end

    if :observed ∈ fieldnames(typeof(sys))
        i = findfirst(x->getname(x.lhs)==name,sys.observed)
        if i !== nothing
            return rename(sys.observed[i].lhs,renamespace(sys.name,name))
        end
    end

    throw(error("Variable $name does not exist"))
end

function renamespace(namespace, x)
    if x isa Num
        renamespace(namespace, value(x))
    elseif istree(x)
        renamespace(namespace, operation(x))(arguments(x)...)
    elseif x isa Sym
        Sym{symtype(x)}(renamespace(namespace,nameof(x)))
    else
        Symbol(namespace,:₊,x)
    end
end

function namespace_variables(sys::AbstractSystem)
    [renamespace(sys.name,x) for x in states(sys)]
end

function namespace_parameters(sys::AbstractSystem)
    [toparam(renamespace(sys.name,x)) for x in parameters(sys)]
end

function namespace_pins(sys::AbstractSystem)
    [renamespace(sys.name,x) for x in pins(sys)]
end

function namespace_equations(sys::AbstractSystem)
    eqs = equations(sys)
    isempty(eqs) && return Equation[]
    iv = independent_variable(sys)
    map(eq->namespace_equation(eq,sys.name,iv), eqs)
end

function namespace_equation(eq::Equation,name,iv)
    _lhs = namespace_expr(eq.lhs,name,iv)
    _rhs = namespace_expr(eq.rhs,name,iv)
    _lhs ~ _rhs
end

function namespace_expr(O::Sym,name,iv)
    isequal(O, iv) ? O : rename(O,renamespace(name,nameof(O)))
end

_symparam(s::Symbolic{T}) where {T} = T
function namespace_expr(O,name,iv) where {T}
    if istree(O)
        renamed = map(a->namespace_expr(a,name,iv), arguments(O))
        if operation(O) isa Sym
            renamed_op = rename(operation(O),renamespace(name,operation(O).name))
            Term{_symparam(O)}(renamed_op,renamed)
        else
            similarterm(O,operation(O),renamed)
        end
    else
        O
    end
end

independent_variable(sys::AbstractSystem) = sys.iv
function states(sys::AbstractSystem)
    unique(isempty(sys.systems) ?
           sys.states :
           [sys.states;reduce(vcat,namespace_variables.(sys.systems))])
end
parameters(sys::AbstractSystem) = isempty(sys.systems) ? sys.ps : [sys.ps;reduce(vcat,namespace_parameters.(sys.systems))]
pins(sys::AbstractSystem) = isempty(sys.systems) ? sys.pins : [sys.pins;reduce(vcat,namespace_pins.(sys.systems))]
function observed(sys::AbstractSystem)
    iv = independent_variable(sys)
    [sys.observed;
     reduce(vcat,
            (map(o->namespace_equation(o, s.name, iv), observed(s)) for s in sys.systems),
            init=Equation[])]
end

function states(sys::AbstractSystem,name::Symbol)
    x = sys.states[findfirst(x->x.name==name,sys.states)]
    s = rename(x,renamespace(sys.name,x.name))
    iv = independent_variable(sys)
    iv === nothing ? iv : s(iv)
end

function parameters(sys::AbstractSystem,name::Symbol)
    x = sys.ps[findfirst(x->x.name==name,sys.ps)]
    rename(x,renamespace(sys.name,x.name))
end

function pins(sys::AbstractSystem,name::Symbol)
    x = sys.pins[findfirst(x->x.name==name,sys.ps)]
    s = rename(x,renamespace(sys.name,x.name))
    iv = independent_variable(sys)
    iv === nothing ? iv : s(iv)
end

lhss(xs) = map(x->x.lhs, xs)
rhss(xs) = map(x->x.rhs, xs)

flatten(sys::AbstractSystem) = sys

function equations(sys::ModelingToolkit.AbstractSystem)
    if isempty(sys.systems)
        return sys.eqs
    else
        eqs = Equation[sys.eqs;
               reduce(vcat,
                      namespace_equations.(sys.systems);
                      init=Equation[])]
        return eqs
    end
end

pins(sys::AbstractSystem,args...) = states(sys, args...)
function states(sys::AbstractSystem,args...)
    name = last(args)
    extra_names = reduce(Symbol,[Symbol(:₊,x.name) for x in args[1:end-1]])
    newname = renamespace(extra_names,name)
    s = rename(x,renamespace(sys.name,newname))
    iv = independent_variable(sys)
    iv === nothing ? iv : s(iv)
end

function parameters(sys::AbstractSystem,args...)
    name = last(args)
    extra_names = reduce(Symbol,[Symbol(:₊,x.name) for x in args[1:end-1]])
    newname = renamespace(extra_names,name)
    rename(x,renamespace(sys.name,newname))
end

function islinear(sys::AbstractSystem)
    rhs = [eq.rhs for eq ∈ equations(sys)]

    all(islinear(r, states(sys)) for r in rhs)
end

struct AbstractSysToExpr
    sys::AbstractSystem
    states::Vector
end
AbstractSysToExpr(sys) = AbstractSysToExpr(sys,states(sys))
function (f::AbstractSysToExpr)(O)
    !istree(O) && return toexpr(O)
    any(isequal(O), f.states) && return operation(O).name  # variables
    if isa(operation(O), Sym)
        return build_expr(:call, Any[operation(O).name; f.(arguments(O))])
    end
    return build_expr(:call, Any[operation(O); f.(arguments(O))])
end

get_default_p(sys) = sys.default_p
get_default_u0(sys) = sys.default_u0
