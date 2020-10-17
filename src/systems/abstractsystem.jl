"""
```julia
calculate_tgrad(sys::AbstractSystem)
```

Calculate the time gradient of a system.

Returns a vector of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_tgrad end

"""
```julia
calculate_gradient(sys::AbstractSystem)
```

Calculate the gradient of a scalar system.

Returns a vector of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_gradient end

"""
```julia
calculate_jacobian(sys::AbstractSystem)
```

Calculate the jacobian matrix of a system.

Returns a matrix of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_jacobian end

"""
```julia
calculate_factorized_W(sys::AbstractSystem)
```

Calculate the factorized W-matrix of a system.

Returns a matrix of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_factorized_W end

"""
```julia
calculate_hessian(sys::AbstractSystem)
```

Calculate the hessian matrix of a scalar system.

Returns a matrix of [`Expression`](@ref) instances. The result from the first
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

getname(x::Sym) = nameof(x)
getname(t::Term) = t.op isa Sym ? getname(t.op) : error("Cannot get name of $t")

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

renamespace(namespace,name) = Symbol(namespace,:₊,name)

function renamespace(namespace, x::Sym)
    Sym{symtype(x)}(renamespace(namespace,x.name))
end

function renamespace(namespace, x::Term)
    renamespace(namespace, x.op)(x.args...)
end

function namespace_variables(sys::AbstractSystem)
    [renamespace(sys.name,x) for x in states(sys)]
end

function namespace_parameters(sys::AbstractSystem)
    [makeparam(renamespace(sys.name,x)) for x in parameters(sys)]
end

function namespace_pins(sys::AbstractSystem)
    [renamespace(sys.name,x) for x in pins(sys)]
end

namespace_equations(sys::AbstractSystem) = namespace_equation.(equations(sys),sys.name,sys.iv.name)

function namespace_equation(eq::Equation,name,ivname)
    _lhs = namespace_expr(eq.lhs,name,ivname)
    _rhs = namespace_expr(eq.rhs,name,ivname)
    _lhs ~ _rhs
end

function namespace_expr(O::Sym,name,ivname)
    O.name == ivname ? O : rename(O,renamespace(name,O.name))
end

function namespace_expr(O::Term,name,ivname)
    if O.op isa Sym
        Term(rename(O.op,renamespace(name,O.op.name)),namespace_expr.(O.args,name,ivname))
    else
        Term(O.op,namespace_expr.(O.args,name,ivname))
    end
end
namespace_expr(O,name,ivname) = O

independent_variable(sys::AbstractSystem) = sys.iv
function states(sys::AbstractSystem)
    unique(isempty(sys.systems) ?
           setdiff(sys.states, value.(sys.pins)) :
           [sys.states;reduce(vcat,namespace_variables.(sys.systems))])
end
parameters(sys::AbstractSystem) = isempty(sys.systems) ? sys.ps : [sys.ps;reduce(vcat,namespace_parameters.(sys.systems))]
pins(sys::AbstractSystem) = isempty(sys.systems) ? sys.pins : [sys.pins;reduce(vcat,namespace_pins.(sys.systems))]
function observed(sys::AbstractSystem)
    [sys.observed;
     reduce(vcat,
            (namespace_equation.(observed(s), s.name, s.iv.name) for s in sys.systems),
            init=Equation[])]
end

function states(sys::AbstractSystem,name::Symbol)
    x = sys.states[findfirst(x->x.name==name,sys.states)]
    rename(x,renamespace(sys.name,x.name))(sys.iv)
end

function parameters(sys::AbstractSystem,name::Symbol)
    x = sys.ps[findfirst(x->x.name==name,sys.ps)]
    rename(x,renamespace(sys.name,x.name))()
end

function pins(sys::AbstractSystem,name::Symbol)
    x = sys.pins[findfirst(x->x.name==name,sys.ps)]
    rename(x,renamespace(sys.name,x.name))(sys.iv())
end

lhss(xs) = map(x->x.lhs, xs)
rhss(xs) = map(x->x.rhs, xs)

function equations(sys::ModelingToolkit.AbstractSystem)
    if isempty(sys.systems)
        return sys.eqs
    else
        eqs = [sys.eqs;
               reduce(vcat,
                      namespace_equations.(sys.systems);
                      init=Equation[])]

        return eqs
    end
end

function states(sys::AbstractSystem,args...)
    name = last(args)
    extra_names = reduce(Symbol,[Symbol(:₊,x.name) for x in args[1:end-1]])
    newname = renamespace(extra_names,name)
    rename(x,renamespace(sys.name,newname))(sys.iv)
end

function parameters(sys::AbstractSystem,args...)
    name = last(args)
    extra_names = reduce(Symbol,[Symbol(:₊,x.name) for x in args[1:end-1]])
    newname = renamespace(extra_names,name)
    rename(x,renamespace(sys.name,newname))()
end

function islinear(sys::AbstractSystem)
    rhs = [eq.rhs for eq ∈ equations(sys)]

    all(islinear(r, states(sys)) for r in rhs)
end

function pins(sys::AbstractSystem,args...)
    name = last(args)
    extra_names = reduce(Symbol,[Symbol(:₊,x.name) for x in args[1:end-1]])
    newname = renamespace(extra_names,name)
    rename(x,renamespace(sys.name,newname))(sys.iv())
end

struct AbstractSysToExpr
    sys::AbstractSystem
    states::Vector
end
AbstractSysToExpr(sys) = AbstractSysToExpr(sys,states(sys))
function (f::AbstractSysToExpr)(O::Term)
    any(isequal(O), f.states) && return O.op.name  # variables
    if isa(O.op, Sym)
        return build_expr(:call, Any[O.op.name; f.(O.args)])
    end
    return build_expr(:call, Any[O.op; f.(O.args)])
end
(f::AbstractSysToExpr)(x) = toexpr(x)

