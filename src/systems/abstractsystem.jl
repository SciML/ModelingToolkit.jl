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

independent_variable(sys::AbstractSystem) = isdefined(sys, :iv) ? getfield(sys, :iv) : nothing

function structure(sys::AbstractSystem)
    s = get_structure(sys)
    s isa SystemStructure || throw(ArgumentError("SystemStructure is not yet initialized, please run `sys = initialize_system_structure(sys)` or `sys = alias_elimination(sys)`."))
    return s
end
for prop in [
             :eqs
             :noiseeqs
             :iv
             :states
             :ps
             :defaults
             :observed
             :tgrad
             :jac
             :Wfact
             :Wfact_t
             :systems
             :structure
             :op
             :equality_constraints
             :inequality_constraints
             :controls
             :loss
             :reduced_states
             :bcs
             :domain
             :depvars
             :indvars
            ]
    fname1 = Symbol(:get_, prop)
    fname2 = Symbol(:has_, prop)
    @eval begin
        $fname1(sys::AbstractSystem) = getfield(sys, $(QuoteNode(prop)))
        $fname2(sys::AbstractSystem) = isdefined(sys, $(QuoteNode(prop)))
    end
end

Setfield.get(obj::AbstractSystem, l::Setfield.PropertyLens{field}) where {field} = getfield(obj, field)
@generated function ConstructionBase.setproperties(obj::AbstractSystem, patch::NamedTuple)
    if issubset(fieldnames(patch), fieldnames(obj))
        args = map(fieldnames(obj)) do fn
            if fn in fieldnames(patch)
                :(patch.$fn)
            else
                :(getfield(obj, $(Meta.quot(fn))))
            end
        end
        return Expr(:block,
            Expr(:meta, :inline),
            Expr(:call,:(constructorof($obj)), args...)
        )
    else
        error("This should never happen. Trying to set $(typeof(obj)) with $patch.")
    end
end

rename(x::AbstractSystem, name) = @set x.name = name

function Base.propertynames(sys::AbstractSystem; private=false)
    if private
        return fieldnames(typeof(sys))
    else
        names = Symbol[]
        for s in get_systems(sys)
            push!(names, getname(s))
        end
        has_states(sys) && for s in get_states(sys)
            push!(names, getname(s))
        end
        has_ps(sys) && for s in get_ps(sys)
            push!(names, getname(s))
        end
        has_observed(sys) && for s in get_observed(sys)
            push!(names, getname(s.lhs))
        end
        return names
    end
end

function Base.getproperty(sys::AbstractSystem, name::Symbol; namespace=true)
    sysname = nameof(sys)
    systems = get_systems(sys)
    if isdefined(sys, name)
        Base.depwarn("`sys.name` like `sys.$name` is deprecated. Use getters like `get_$name` instead.", "sys.$name")
        return getfield(sys, name)
    elseif !isempty(systems)
        i = findfirst(x->nameof(x)==name,systems)
        if i !== nothing
            return namespace ? rename(systems[i],renamespace(sysname,name)) : systems[i]
        end
    end

    sts = get_states(sys)
    i = findfirst(x->getname(x) == name, sts)

    if i !== nothing
        return namespace ? rename(sts[i],renamespace(sysname,name)) : sts[i]
    end

    if has_ps(sys)
        ps = get_ps(sys)
        i = findfirst(x->getname(x) == name,ps)
        if i !== nothing
            return namespace ? rename(ps[i],renamespace(sysname,name)) : ps[i]
        end
    end

    if has_observed(sys)
        obs = get_observed(sys)
        i = findfirst(x->getname(x.lhs)==name,obs)
        if i !== nothing
            return namespace ? rename(obs[i].lhs,renamespace(sysname,name)) : obs[i]
        end
    end

    throw(error("Variable $name does not exist"))
end

function Base.setproperty!(sys::AbstractSystem, prop::Symbol, val)
    # We use this weird syntax because `parameters` and `states` calls are
    # potentially expensive.
    if (
        params = parameters(sys);
        idx = findfirst(s->getname(s) == prop, params);
        idx !== nothing;
       )
        get_defaults(sys)[params[idx]] = value(val)
    elseif (
            sts = states(sys);
            idx = findfirst(s->getname(s) == prop, sts);
            idx !== nothing;
           )
        get_defaults(sys)[sts[idx]] = value(val)
    else
        setfield!(sys, prop, val)
    end
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

function namespace_defaults(sys)
    defs = defaults(sys)
    Dict((isparameter(k) ? parameters(sys, k) : states(sys, k)) => namespace_expr(defs[k], nameof(sys), independent_variable(sys)) for k in keys(defs))
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

Base.@deprecate default_u0(x) defaults(x) false
Base.@deprecate default_p(x) defaults(x) false
function defaults(sys::AbstractSystem)
    systems = get_systems(sys)
    defs = get_defaults(sys)
    isempty(systems) ? defs : mapreduce(namespace_defaults, merge, systems; init=defs)
end

states(sys::AbstractSystem, v) = renamespace(nameof(sys), v)
parameters(sys::AbstractSystem, v) = toparam(states(sys, v))
for f in [:states, :parameters]
    @eval $f(sys::AbstractSystem, vs::AbstractArray) = map(v->$f(sys, v), vs)
end

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

function Base.show(io::IO, ::MIME"text/plain", sys::AbstractSystem)
    eqs = equations(sys)
    if eqs isa AbstractArray
        Base.printstyled(io, "Model $(nameof(sys)) with $(length(eqs)) equations\n"; bold=true)
    else
        Base.printstyled(io, "Model $(nameof(sys))\n"; bold=true)
    end
    # The reduced equations are usually very long. It's not that useful to print
    # them.
    #Base.print_matrix(io, eqs)
    #println(io)

    rows = first(displaysize(io)) ÷ 5
    limit = get(io, :limit, false)

    vars = states(sys); nvars = length(vars)
    Base.printstyled(io, "States ($nvars):"; bold=true)
    nrows = min(nvars, limit ? rows : nvars)
    limited = nrows < length(vars)
    defs = has_defaults(sys) ? defaults(sys) : nothing
    for i in 1:nrows
        s = vars[i]
        print(io, "\n  ", s)

        if defs !== nothing
            val = get(defs, s, nothing)
            if val !== nothing
                print(io, " [defaults to $val]")
            end
        end
    end
    limited && print(io, "\n⋮")
    println(io)

    vars = parameters(sys); nvars = length(vars)
    Base.printstyled(io, "Parameters ($nvars):"; bold=true)
    nrows = min(nvars, limit ? rows : nvars)
    limited = nrows < length(vars)
    for i in 1:nrows
        s = vars[i]
        print(io, "\n  ", s)

        if defs !== nothing
            val = get(defs, s, nothing)
            if val !== nothing
                print(io, " [defaults to $val]")
            end
        end
    end
    limited && print(io, "\n⋮")

    if has_structure(sys)
        s = get_structure(sys)
        if s !== nothing
            Base.printstyled(io, "\nIncidence matrix:"; color=:magenta)
            show(io, incidence_matrix(s.graph, Num(Sym{Real}(:×))))
        end
    end
    return nothing
end

function _named(expr)
    if !(expr isa Expr && expr.head === :(=) && expr.args[2].head === :call)
        throw(ArgumentError("expression should be of the form `sys = foo(a, b)`"))
    end
    name, call = expr.args

    has_kw = false
    if length(call.args) >= 2 && call.args[2] isa Expr
        # canonicalize to use `:parameters`
        if call.args[2].head === :kw
            call.args[2] = Expr(:parameters, Expr(:kw, call.args[2].args...))
            has_kw = true
        elseif call.args[2].head === :parameters
            has_kw = true
        end
    end

    if !has_kw
        param = Expr(:parameters)
        if length(call.args) == 1
            push!(call.args, param)
        else
            insert!(call.args, 2, param)
        end
    end

    kws = call.args[2].args

    if !any(kw->(kw isa Symbol ? kw : kw.args[1]) == :name, kws) # don't overwrite `name` kwarg
        pushfirst!(kws, Expr(:kw, :name, Meta.quot(name)))
    end
    :($name = $call)
end

"""
$(SIGNATURES)

Rewrite `@named y = foo(x)` to `y = foo(x; name=:y)`.
"""
macro named(expr)
    esc(_named(expr))
end

function _nonamespace(expr)
    if Meta.isexpr(expr, :.)
        return :($getproperty($(map(_nonamespace, expr.args)...); namespace=false))
    elseif expr isa Expr && !isempty(expr.args)
        return Expr(expr.head, map(_nonamespace, expr.args)...)
    else
        expr
    end
end

"""
$(SIGNATURES)

Rewrite `@nonamespace a.b.c` to
`getproperty(getproperty(a, :b; namespace = false), :c; namespace = false)`.
"""
macro nonamespace(expr)
    esc(_nonamespace(expr))
end

"""
$(SIGNATURES)

Structurally simplify algebraic equations in a system and compute the
topological sort of the observed equations.
"""
function structural_simplify(sys::AbstractSystem)
    sys = tearing(alias_elimination(sys))
    fullstates = [get_reduced_states(sys); states(sys)]
    @set! sys.observed = topsort_equations(observed(sys), fullstates)
    return sys
end

@latexrecipe function f(sys::AbstractSystem)
    return latexify(equations(sys))
end

Base.show(io::IO, ::MIME"text/latex", x::AbstractSystem) = print(io, latexify(x))

struct InvalidSystemException <: Exception
    msg::String
end
Base.showerror(io::IO, e::InvalidSystemException) = print(io, "InvalidSystemException: ", e.msg)
