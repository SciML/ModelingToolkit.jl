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

Base.nameof(sys::AbstractSystem) = getfield(sys, :name)

function getname(t)
    if istree(t)
        operation(t) isa Sym ? getname(operation(t)) : error("Cannot get name of $t")
    else
        nameof(t)
    end
end

independent_variable(sys::AbstractSystem) = getfield(sys, :iv)

function structure(sys::AbstractSystem)
    s = get_structure(sys)
    s isa SystemStructure || throw(ArgumentError("SystemStructure is not yet initialized, please run `sys = initialize_system_structure(sys)` or `sys = alias_elimination(sys)`."))
    return s
end

for prop in [:eqs, :iv, :states, :ps, :default_p, :default_u0, :observed, :tgrad, :jac, :Wfact, :Wfact_t, :systems, :structure]
    fname = Symbol(:get_, prop)
    @eval begin
        $fname(sys::AbstractSystem) = getfield(sys, $(QuoteNode(prop)))
    end
end

function Base.getproperty(sys::AbstractSystem, name::Symbol)
    sysname = nameof(sys)
    systems = get_systems(sys)
    if isdefined(sys, name)
        Base.depwarn("`sys.name` like `sys.$name` is deprecated. Use getters like `get_$name` instead.", "sys.$name")
        return getfield(sys, name)
    elseif !isempty(systems)
        i = findfirst(x->nameof(x)==name,systems)
        if i !== nothing
            return rename(systems[i],renamespace(sysname,name))
        end
    end

    sts = get_states(sys)
    i = findfirst(x->getname(x) == name, sts)

    if i !== nothing
        return rename(sts[i],renamespace(sysname,name))
    end

    if isdefined(sys, :ps)
        ps = get_ps(sys)
        i = findfirst(x->getname(x) == name,ps)
        if i !== nothing
            return rename(ps[i],renamespace(sysname,name))
        end
    end

    if isdefined(sys, :observed)
        obs = get_observed(sys)
        i = findfirst(x->getname(x.lhs)==name,obs)
        if i !== nothing
            return rename(obs[i].lhs,renamespace(sysname,name))
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

namespace_variables(sys::AbstractSystem) = states(sys, states(sys))
namespace_parameters(sys::AbstractSystem) = parameters(sys, parameters(sys))

function namespace_default_u0(sys)
    d_u0 = default_u0(sys)
    Dict(states(sys, k) => d_u0[k] for k in keys(d_u0))
end

function namespace_default_p(sys)
    d_p = default_p(sys)
    Dict(parameters(sys, k) => d_p[k] for k in keys(d_p))
end

function namespace_equations(sys::AbstractSystem)
    eqs = equations(sys)
    isempty(eqs) && return Equation[]
    iv = independent_variable(sys)
    map(eq->namespace_equation(eq,nameof(sys),iv), eqs)
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
            renamed_op = rename(operation(O),renamespace(name,nameof(operation(O))))
            Term{_symparam(O)}(renamed_op,renamed)
        else
            similarterm(O,operation(O),renamed)
        end
    else
        O
    end
end

function states(sys::AbstractSystem)
    sts = get_states(sys)
    systems = get_systems(sys)
    unique(isempty(systems) ?
           sts :
           [sts;reduce(vcat,namespace_variables.(systems))])
end
function parameters(sys::AbstractSystem)
    ps = get_ps(sys)
    systems = get_systems(sys)
    isempty(systems) ? ps : [ps;reduce(vcat,namespace_parameters.(systems))]
end
function observed(sys::AbstractSystem)
    iv = independent_variable(sys)
    obs = get_observed(sys)
    systems = get_systems(sys)
    [obs;
     reduce(vcat,
            (map(o->namespace_equation(o, nameof(s), iv), observed(s)) for s in systems),
            init=Equation[])]
end

function default_u0(sys::AbstractSystem)
    systems = get_systems(sys)
    d_u0 = get_default_u0(sys)
    isempty(systems) ? d_u0 : mapreduce(namespace_default_u0, merge, systems; init=d_u0)
end

function default_p(sys::AbstractSystem)
    systems = get_systems(sys)
    d_p = get_default_p(sys)
    isempty(systems) ? d_p : mapreduce(namespace_default_p, merge, systems; init=d_p)
end

states(sys::AbstractSystem, v) = renamespace(nameof(sys), v)
parameters(sys::AbstractSystem, v) = toparam(states(sys, v))
for f in [:states, :parameters]
    @eval $f(sys::AbstractSystem, vs::AbstractArray) = map(v->$f(sys, v), vs)
end

lhss(xs) = map(x->x.lhs, xs)
rhss(xs) = map(x->x.rhs, xs)

flatten(sys::AbstractSystem) = sys

function equations(sys::ModelingToolkit.AbstractSystem)
    eqs = get_eqs(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return eqs
    else
        eqs = Equation[eqs;
               reduce(vcat,
                      namespace_equations.(get_systems(sys));
                      init=Equation[])]
        return eqs
    end
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
    any(isequal(O), f.states) && return nameof(operation(O))  # variables
    if isa(operation(O), Sym)
        return build_expr(:call, Any[nameof(operation(O)); f.(arguments(O))])
    end
    return build_expr(:call, Any[operation(O); f.(arguments(O))])
end
