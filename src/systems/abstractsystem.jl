"""
```julia
calculate_tgrad(sys::AbstractTimeDependentSystem)
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
calculate_control_jacobian(sys::AbstractSystem)
```

Calculate the jacobian matrix of a system with respect to the system's controls.

Returns a matrix of [`Num`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_control_jacobian end

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
generate_tgrad(sys::AbstractTimeDependentSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)
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

mutable struct Substitutions
    subs::Vector{Equation}
    deps::Vector{Vector{Int}}
    subed_eqs::Union{Nothing, Vector{Equation}}
end
Substitutions(subs, deps) = Substitutions(subs, deps, nothing)

Base.nameof(sys::AbstractSystem) = getfield(sys, :name)

#Deprecated
function independent_variable(sys::AbstractSystem)
    Base.depwarn("`independent_variable` is deprecated. Use `get_iv` or `independent_variables` instead.",
                 :independent_variable)
    isdefined(sys, :iv) ? getfield(sys, :iv) : nothing
end

#Treat the result as a vector of symbols always
function independent_variables(sys::AbstractSystem)
    systype = typeof(sys)
    @warn "Please declare ($systype) as a subtype of `AbstractTimeDependentSystem`, `AbstractTimeIndependentSystem` or `AbstractMultivariateSystem`."
    if isdefined(sys, :iv)
        return [getfield(sys, :iv)]
    elseif isdefined(sys, :ivs)
        return getfield(sys, :ivs)
    else
        return []
    end
end

independent_variables(sys::AbstractTimeDependentSystem) = [getfield(sys, :iv)]
independent_variables(sys::AbstractTimeIndependentSystem) = []
independent_variables(sys::AbstractMultivariateSystem) = getfield(sys, :ivs)

const NULL_AFFECT = Equation[]
struct SymbolicContinuousCallback
    eqs::Vector{Equation}
    affect::Vector{Equation}
    function SymbolicContinuousCallback(eqs::Vector{Equation}, affect = NULL_AFFECT)
        new(eqs, affect)
    end # Default affect to nothing
end

function Base.:(==)(e1::SymbolicContinuousCallback, e2::SymbolicContinuousCallback)
    isequal(e1.eqs, e2.eqs) && isequal(e1.affect, e2.affect)
end
Base.isempty(cb::SymbolicContinuousCallback) = isempty(cb.eqs)
function Base.hash(cb::SymbolicContinuousCallback, s::UInt)
    s = foldr(hash, cb.eqs, init = s)
    foldr(hash, cb.affect, init = s)
end

to_equation_vector(eq::Equation) = [eq]
to_equation_vector(eqs::Vector{Equation}) = eqs
function to_equation_vector(eqs::Vector{Any})
    isempty(eqs) || error("This should never happen")
    Equation[]
end

function SymbolicContinuousCallback(args...)
    SymbolicContinuousCallback(to_equation_vector.(args)...)
end # wrap eq in vector
SymbolicContinuousCallback(p::Pair) = SymbolicContinuousCallback(p[1], p[2])
SymbolicContinuousCallback(cb::SymbolicContinuousCallback) = cb # passthrough

SymbolicContinuousCallbacks(cb::SymbolicContinuousCallback) = [cb]
SymbolicContinuousCallbacks(cbs::Vector{<:SymbolicContinuousCallback}) = cbs
SymbolicContinuousCallbacks(cbs::Vector) = SymbolicContinuousCallback.(cbs)
function SymbolicContinuousCallbacks(ve::Vector{Equation})
    SymbolicContinuousCallbacks(SymbolicContinuousCallback(ve))
end
function SymbolicContinuousCallbacks(others)
    SymbolicContinuousCallbacks(SymbolicContinuousCallback(others))
end
SymbolicContinuousCallbacks(::Nothing) = SymbolicContinuousCallbacks(Equation[])

equations(cb::SymbolicContinuousCallback) = cb.eqs
function equations(cbs::Vector{<:SymbolicContinuousCallback})
    reduce(vcat, [equations(cb) for cb in cbs])
end
affect_equations(cb::SymbolicContinuousCallback) = cb.affect
function affect_equations(cbs::Vector{SymbolicContinuousCallback})
    reduce(vcat, [affect_equations(cb) for cb in cbs])
end
namespace_equation(cb::SymbolicContinuousCallback, s)::SymbolicContinuousCallback = SymbolicContinuousCallback(namespace_equation.(equations(cb),
                                                                                                                                   (s,)),
                                                                                                               namespace_equation.(affect_equations(cb),
                                                                                                                                   (s,)))

for prop in [:eqs
             :noiseeqs
             :iv
             :states
             :ps
             :var_to_name
             :ctrls
             :defaults
             :observed
             :tgrad
             :jac
             :ctrl_jac
             :Wfact
             :Wfact_t
             :systems
             :structure
             :op
             :equality_constraints
             :inequality_constraints
             :controls
             :loss
             :bcs
             :domain
             :ivs
             :dvs
             :connector_type
             :connections
             :preface
             :torn_matching
             :tearing_state
             :substitutions]
    fname1 = Symbol(:get_, prop)
    fname2 = Symbol(:has_, prop)
    @eval begin
        $fname1(sys::AbstractSystem) = getfield(sys, $(QuoteNode(prop)))
        $fname2(sys::AbstractSystem) = isdefined(sys, $(QuoteNode(prop)))
    end
end

const EMPTY_TGRAD = Vector{Num}(undef, 0)
const EMPTY_JAC = Matrix{Num}(undef, 0, 0)
function invalidate_cache!(sys::AbstractSystem)
    if has_tgrad(sys)
        get_tgrad(sys)[] = EMPTY_TGRAD
    end
    if has_jac(sys)
        get_jac(sys)[] = EMPTY_JAC
    end
    if has_ctrl_jac(sys)
        get_ctrl_jac(sys)[] = EMPTY_JAC
    end
    if has_Wfact(sys)
        get_Wfact(sys)[] = EMPTY_JAC
    end
    if has_Wfact_t(sys)
        get_Wfact_t(sys)[] = EMPTY_JAC
    end
    return sys
end

function Setfield.get(obj::AbstractSystem, ::Setfield.PropertyLens{field}) where {field}
    getfield(obj, field)
end
@generated function ConstructionBase.setproperties(obj::AbstractSystem, patch::NamedTuple)
    if issubset(fieldnames(patch), fieldnames(obj))
        args = map(fieldnames(obj)) do fn
            if fn in fieldnames(patch)
                :(patch.$fn)
            else
                :(getfield(obj, $(Meta.quot(fn))))
            end
        end
        kwarg = :($(Expr(:kw, :checks, false))) # Inputs should already be checked
        return Expr(:block,
                    Expr(:meta, :inline),
                    Expr(:call, :(constructorof($obj)), args..., kwarg))
    else
        error("This should never happen. Trying to set $(typeof(obj)) with $patch.")
    end
end

rename(x::AbstractSystem, name) = @set x.name = name

function Base.propertynames(sys::AbstractSystem; private = false)
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

function Base.getproperty(sys::AbstractSystem, name::Symbol; namespace = true)
    wrap(getvar(sys, name; namespace = namespace))
end
function getvar(sys::AbstractSystem, name::Symbol; namespace = false)
    systems = get_systems(sys)
    if isdefined(sys, name)
        Base.depwarn("`sys.name` like `sys.$name` is deprecated. Use getters like `get_$name` instead.",
                     "sys.$name")
        return getfield(sys, name)
    elseif !isempty(systems)
        i = findfirst(x -> nameof(x) == name, systems)
        if i !== nothing
            return namespace ? renamespace(sys, systems[i]) : systems[i]
        end
    end

    if has_var_to_name(sys)
        avs = get_var_to_name(sys)
        v = get(avs, name, nothing)
        v === nothing || return namespace ? renamespace(sys, v) : v
    else
        sts = get_states(sys)
        i = findfirst(x -> getname(x) == name, sts)
        if i !== nothing
            return namespace ? renamespace(sys, sts[i]) : sts[i]
        end

        if has_ps(sys)
            ps = get_ps(sys)
            i = findfirst(x -> getname(x) == name, ps)
            if i !== nothing
                return namespace ? renamespace(sys, ps[i]) : ps[i]
            end
        end
    end

    sts = get_states(sys)
    i = findfirst(x -> getname(x) == name, sts)

    if has_observed(sys)
        obs = get_observed(sys)
        i = findfirst(x -> getname(x.lhs) == name, obs)
        if i !== nothing
            return namespace ? renamespace(sys, obs[i]) : obs[i]
        end
    end

    throw(ArgumentError("System $(nameof(sys)): variable $name does not exist"))
end

function Base.setproperty!(sys::AbstractSystem, prop::Symbol, val)
    # We use this weird syntax because `parameters` and `states` calls are
    # potentially expensive.
    if (params = parameters(sys);
        idx = findfirst(s -> getname(s) == prop, params);
        idx !== nothing)
        get_defaults(sys)[params[idx]] = value(val)
    elseif (sts = states(sys);
            idx = findfirst(s -> getname(s) == prop, sts);
            idx !== nothing)
        get_defaults(sys)[sts[idx]] = value(val)
    else
        setfield!(sys, prop, val)
    end
end

abstract type SymScope end

struct LocalScope <: SymScope end
LocalScope(sym::Union{Num, Symbolic}) = setmetadata(sym, SymScope, LocalScope())

struct ParentScope <: SymScope
    parent::SymScope
end
function ParentScope(sym::Union{Num, Symbolic})
    setmetadata(sym, SymScope, ParentScope(getmetadata(value(sym), SymScope, LocalScope())))
end

struct GlobalScope <: SymScope end
GlobalScope(sym::Union{Num, Symbolic}) = setmetadata(sym, SymScope, GlobalScope())

renamespace(sys, eq::Equation) = namespace_equation(eq, sys)

renamespace(names::AbstractVector, x) = foldr(renamespace, names, init = x)
function renamespace(sys, x)
    sys === nothing && return x
    x = unwrap(x)
    if x isa Symbolic
        if isdifferential(x)
            return similarterm(x, operation(x), Any[renamespace(sys, only(arguments(x)))])
        end
        let scope = getmetadata(x, SymScope, LocalScope())
            if scope isa LocalScope
                rename(x, renamespace(getname(sys), getname(x)))
            elseif scope isa ParentScope
                setmetadata(x, SymScope, scope.parent)
            else # GlobalScope
                x
            end
        end
    elseif x isa AbstractSystem
        rename(x, renamespace(sys, nameof(x)))
    else
        Symbol(getname(sys), :₊, x)
    end
end

namespace_variables(sys::AbstractSystem) = states(sys, states(sys))
namespace_parameters(sys::AbstractSystem) = parameters(sys, parameters(sys))
namespace_controls(sys::AbstractSystem) = controls(sys, controls(sys))

function namespace_defaults(sys)
    defs = defaults(sys)
    Dict((isparameter(k) ? parameters(sys, k) : states(sys, k)) => namespace_expr(defs[k],
                                                                                  sys)
         for k in keys(defs))
end

function namespace_equations(sys::AbstractSystem)
    eqs = equations(sys)
    isempty(eqs) && return Equation[]
    map(eq -> namespace_equation(eq, sys), eqs)
end

function namespace_equation(eq::Equation, sys, n = nameof(sys))
    _lhs = namespace_expr(eq.lhs, sys, n)
    _rhs = namespace_expr(eq.rhs, sys, n)
    _lhs ~ _rhs
end

function namespace_assignment(eq::Assignment, sys)
    _lhs = namespace_expr(eq.lhs, sys)
    _rhs = namespace_expr(eq.rhs, sys)
    Assignment(_lhs, _rhs)
end

function namespace_expr(O, sys, n = nameof(sys)) where {T}
    ivs = independent_variables(sys)
    O = unwrap(O)
    if any(isequal(O), ivs)
        return O
    elseif isvariable(O)
        renamespace(n, O)
    elseif istree(O)
        renamed = map(a -> namespace_expr(a, sys, n), arguments(O))
        if symtype(operation(O)) <: FnType
            renamespace(n, O)
        else
            similarterm(O, operation(O), renamed)
        end
    elseif O isa Array
        map(o -> namespace_expr(o, sys, n), O)
    else
        O
    end
end

function states(sys::AbstractSystem)
    sts = get_states(sys)
    systems = get_systems(sys)
    unique(isempty(systems) ?
           sts :
           [sts; reduce(vcat, namespace_variables.(systems))])
end

function parameters(sys::AbstractSystem)
    ps = get_ps(sys)
    systems = get_systems(sys)
    unique(isempty(systems) ? ps : [ps; reduce(vcat, namespace_parameters.(systems))])
end

function controls(sys::AbstractSystem)
    ctrls = get_ctrls(sys)
    systems = get_systems(sys)
    isempty(systems) ? ctrls : [ctrls; reduce(vcat, namespace_controls.(systems))]
end

function observed(sys::AbstractSystem)
    obs = get_observed(sys)
    systems = get_systems(sys)
    [obs;
     reduce(vcat,
            (map(o -> namespace_equation(o, s), observed(s)) for s in systems),
            init = Equation[])]
end

function continuous_events(sys::AbstractSystem)
    obs = get_continuous_events(sys)
    filter(!isempty, obs)
    systems = get_systems(sys)
    cbs = [obs;
           reduce(vcat,
                  (map(o -> namespace_equation(o, s), continuous_events(s))
                   for s in systems),
                  init = SymbolicContinuousCallback[])]
    filter(!isempty, cbs)
end

Base.@deprecate default_u0(x) defaults(x) false
Base.@deprecate default_p(x) defaults(x) false
function defaults(sys::AbstractSystem)
    systems = get_systems(sys)
    defs = get_defaults(sys)
    # `mapfoldr` is really important!!! We should prefer the base model for
    # defaults, because people write:
    #
    # `compose(ODESystem(...; defaults=defs), ...)`
    #
    # Thus, right associativity is required and crucial for correctness.
    isempty(systems) ? defs : mapfoldr(namespace_defaults, merge, systems; init = defs)
end

states(sys::AbstractSystem, v) = renamespace(sys, v)
parameters(sys::AbstractSystem, v) = toparam(states(sys, v))
for f in [:states, :parameters]
    @eval $f(sys::AbstractSystem, vs::AbstractArray) = map(v -> $f(sys, v), vs)
end

flatten(sys::AbstractSystem, args...) = sys

function equations(sys::ModelingToolkit.AbstractSystem)
    eqs = get_eqs(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return eqs
    else
        eqs = Equation[eqs;
                       reduce(vcat,
                              namespace_equations.(get_systems(sys));
                              init = Equation[])]
        return eqs
    end
end

function preface(sys::ModelingToolkit.AbstractSystem)
    has_preface(sys) || return nothing
    pre = get_preface(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return pre
    else
        pres = pre === nothing ? [] : pre
        for sys in systems
            pre = get_preface(sys)
            pre === nothing && continue
            for eq in pre
                push!(pres, namespace_assignment(eq, sys))
            end
        end
        return isempty(pres) ? nothing : pres
    end
end

function islinear(sys::AbstractSystem)
    rhs = [eq.rhs for eq in equations(sys)]

    all(islinear(r, states(sys)) for r in rhs)
end

function isaffine(sys::AbstractSystem)
    rhs = [eq.rhs for eq in equations(sys)]

    all(isaffine(r, states(sys)) for r in rhs)
end

struct AbstractSysToExpr
    sys::AbstractSystem
    states::Vector
end
AbstractSysToExpr(sys) = AbstractSysToExpr(sys, states(sys))
function (f::AbstractSysToExpr)(O)
    !istree(O) && return toexpr(O)
    any(isequal(O), f.states) && return nameof(operation(O))  # variables
    if isa(operation(O), Sym)
        return build_expr(:call, Any[nameof(operation(O)); f.(arguments(O))])
    end
    return build_expr(:call, Any[operation(O); f.(arguments(O))])
end

###
### System utils
###
function push_vars!(stmt, name, typ, vars)
    isempty(vars) && return
    vars_expr = Expr(:macrocall, typ, nothing)
    for s in vars
        if istree(s)
            f = nameof(operation(s))
            args = arguments(s)
            ex = :($f($(args...)))
        else
            ex = nameof(s)
        end
        push!(vars_expr.args, ex)
    end
    push!(stmt, :($name = $collect($vars_expr)))
    return
end

function round_trip_expr(t, var2name)
    name = get(var2name, t, nothing)
    name !== nothing && return name
    t isa Sym && return nameof(t)
    istree(t) || return t
    f = round_trip_expr(operation(t), var2name)
    args = map(Base.Fix2(round_trip_expr, var2name), arguments(t))
    return :($f($(args...)))
end

function round_trip_eq(eq::Equation, var2name)
    if eq.lhs isa Connection
        syss = get_systems(eq.rhs)
        call = Expr(:call, connect)
        for sys in syss
            strs = split(string(nameof(sys)), "₊")
            s = Symbol(strs[1])
            for st in strs[2:end]
                s = Expr(:., s, Meta.quot(Symbol(st)))
            end
            push!(call.args, s)
        end
        call
    else
        Expr(:call, (~), round_trip_expr(eq.lhs, var2name),
             round_trip_expr(eq.rhs, var2name))
    end
end

function push_eqs!(stmt, eqs, var2name)
    eqs_name = gensym(:eqs)
    eqs_expr = Expr(:vcat)
    eqs_blk = Expr(:(=), eqs_name, eqs_expr)
    for eq in eqs
        push!(eqs_expr.args, round_trip_eq(eq, var2name))
    end

    push!(stmt, eqs_blk)
    return eqs_name
end

function push_defaults!(stmt, defs, var2name)
    defs_name = gensym(:defs)
    defs_expr = Expr(:call, Dict)
    defs_blk = Expr(:(=), defs_name, defs_expr)
    for d in defs
        n = round_trip_expr(d.first, var2name)
        v = round_trip_expr(d.second, var2name)
        push!(defs_expr.args, :($(=>)($n, $v)))
    end

    push!(stmt, defs_blk)
    return defs_name
end

###
### System I/O
###
function toexpr(sys::AbstractSystem)
    sys = flatten(sys)
    expr = Expr(:block)
    stmt = expr.args

    name = Meta.quot(nameof(sys))
    ivs = independent_variables(sys)
    ivname = gensym(:iv)
    for iv in ivs
        ivname = gensym(:iv)
        push!(stmt, :($ivname = (@variables $(getname(iv)))[1]))
    end

    stsname = gensym(:sts)
    sts = states(sys)
    push_vars!(stmt, stsname, Symbol("@variables"), sts)
    psname = gensym(:ps)
    ps = parameters(sys)
    push_vars!(stmt, psname, Symbol("@parameters"), ps)

    var2name = Dict{Any, Symbol}()
    for v in Iterators.flatten((sts, ps))
        var2name[v] = getname(v)
    end

    eqs_name = push_eqs!(stmt, equations(sys), var2name)
    defs_name = push_defaults!(stmt, defaults(sys), var2name)

    if sys isa ODESystem
        iv = get_iv(sys)
        ivname = gensym(:iv)
        push!(stmt, :($ivname = (@variables $(getname(iv)))[1]))
        push!(stmt,
              :($ODESystem($eqs_name, $ivname, $stsname, $psname; defaults = $defs_name,
                           name = $name, checks = false)))
    elseif sys isa NonlinearSystem
        push!(stmt,
              :($NonlinearSystem($eqs_name, $stsname, $psname; defaults = $defs_name,
                                 name = $name, checks = false)))
    end

    striplines(expr) # keeping the line numbers is never helpful
end

Base.write(io::IO, sys::AbstractSystem) = write(io, readable_code(toexpr(sys)))

function get_or_construct_tearing_state(sys)
    if has_tearing_state(sys)
        state = get_tearing_state(sys)
        if state === nothing
            state = TearingState(sys)
        end
    else
        state = nothing
    end
    state
end

# TODO: what about inputs?
function n_extra_equations(sys::AbstractSystem)
    isconnector(sys) && return length(get_states(sys))
    sys, csets = generate_connection_set(sys)
    ceqs, instream_csets = generate_connection_equations_and_stream_connections(csets)
    n_outer_stream_variables = 0
    for cset in instream_csets
        n_outer_stream_variables += count(x -> x.isouter, cset.set)
    end

    #n_toplevel_unused_flows = 0
    #toplevel_flows = Set()
    #for cset in csets
    #    e1 = first(cset.set)
    #    e1.sys.namespace === nothing || continue
    #    for e in cset.set
    #        get_connection_type(e.v) === Flow || continue
    #        push!(toplevel_flows, e.v)
    #    end
    #end
    #for m in get_systems(sys)
    #    isconnector(m) || continue
    #    n_toplevel_unused_flows += count(x->get_connection_type(x) === Flow && !(x in toplevel_flows), get_states(m))
    #end

    nextras = n_outer_stream_variables + length(ceqs)
end

function Base.show(io::IO, mime::MIME"text/plain", sys::AbstractSystem)
    eqs = equations(sys)
    vars = states(sys)
    nvars = length(vars)
    if eqs isa AbstractArray && eltype(eqs) <: Equation
        neqs = count(eq -> !(eq.lhs isa Connection), eqs)
        Base.printstyled(io, "Model $(nameof(sys)) with $neqs "; bold = true)
        nextras = n_extra_equations(sys)
        if nextras > 0
            Base.printstyled(io, "("; bold = true)
            Base.printstyled(io, neqs + nextras; bold = true, color = :magenta)
            Base.printstyled(io, ") "; bold = true)
        end
        Base.printstyled(io, "equations\n"; bold = true)
    else
        Base.printstyled(io, "Model $(nameof(sys))\n"; bold = true)
    end
    # The reduced equations are usually very long. It's not that useful to print
    # them.
    #Base.print_matrix(io, eqs)
    #println(io)

    rows = first(displaysize(io)) ÷ 5
    limit = get(io, :limit, false)

    Base.printstyled(io, "States ($nvars):"; bold = true)
    nrows = min(nvars, limit ? rows : nvars)
    limited = nrows < length(vars)
    defs = has_defaults(sys) ? defaults(sys) : nothing
    for i in 1:nrows
        s = vars[i]
        print(io, "\n  ", s)

        if defs !== nothing
            val = get(defs, s, nothing)
            if val !== nothing
                print(io, " [defaults to ")
                show(IOContext(io, :compact => true, :limit => true,
                               :displaysize => (1, displaysize(io)[2])), val)
                print(io, "]")
            end
        end
    end
    limited && print(io, "\n⋮")
    println(io)

    vars = parameters(sys)
    nvars = length(vars)
    Base.printstyled(io, "Parameters ($nvars):"; bold = true)
    nrows = min(nvars, limit ? rows : nvars)
    limited = nrows < length(vars)
    for i in 1:nrows
        s = vars[i]
        print(io, "\n  ", s)

        if defs !== nothing
            val = get(defs, s, nothing)
            if val !== nothing
                print(io, " [defaults to ")
                show(IOContext(io, :compact => true, :limit => true,
                               :displaysize => (1, displaysize(io)[2])), val)
                print(io, "]")
            end
        end
    end
    limited && print(io, "\n⋮")

    if has_torn_matching(sys) && has_tearing_state(sys)
        # If the system can take a torn matching, then we can initialize a tearing
        # state on it. Do so and get show the structure.
        state = get_tearing_state(sys)
        if state !== nothing
            Base.printstyled(io, "\nIncidence matrix:"; color = :magenta)
            show(io, mime, incidence_matrix(state.structure.graph, Num(Sym{Real}(:×))))
        end
    end
    return nothing
end

function split_assign(expr)
    if !(expr isa Expr && expr.head === :(=) && expr.args[2].head === :call)
        throw(ArgumentError("expression should be of the form `sys = foo(a, b)`"))
    end
    name, call = expr.args
end

function _named(name, call, runtime = false)
    has_kw = false
    call isa Expr || throw(Meta.ParseError("The rhs must be an Expr. Got $call."))
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

    if !any(kw -> (kw isa Symbol ? kw : kw.args[1]) == :name, kws) # don't overwrite `name` kwarg
        pushfirst!(kws, Expr(:kw, :name, runtime ? name : Meta.quot(name)))
    end
    call
end

function _named_idxs(name::Symbol, idxs, call)
    if call.head !== :->
        throw(ArgumentError("Not an anonymous function"))
    end
    if !isa(call.args[1], Symbol)
        throw(ArgumentError("not a single-argument anonymous function"))
    end
    sym, ex = call.args
    ex = Base.Cartesian.poplinenum(ex)
    ex = _named(:(Symbol($(Meta.quot(name)), :_, $sym)), ex, true)
    ex = Base.Cartesian.poplinenum(ex)
    :($name = $map($sym -> $ex, $idxs))
end

function check_name(name)
    name isa Symbol ||
        throw(Meta.ParseError("The lhs must be a symbol (a) or a ref (a[1:10]). Got $name."))
end

"""
    @named y = foo(x)
    @named y[1:10] = foo(x)
    @named y 1:10 i -> foo(x*i)

Rewrite `@named y = foo(x)` to `y = foo(x; name=:y)`.

Rewrite `@named y[1:10] = foo(x)` to `y = map(i′->foo(x; name=Symbol(:y_, i′)), 1:10)`.

Rewrite `@named y 1:10 i -> foo(x*i)` to `y = map(i->foo(x*i; name=Symbol(:y_, i)), 1:10)`.

Examples:
```julia
julia> using ModelingToolkit

julia> foo(i; name) = i, name
foo (generic function with 1 method)

julia> x = 41
41

julia> @named y = foo(x)
(41, :y)

julia> @named y[1:3] = foo(x)
3-element Vector{Tuple{Int64, Symbol}}:
 (41, :y_1)
 (41, :y_2)
 (41, :y_3)

julia> @named y 1:3 i -> foo(x*i)
3-element Vector{Tuple{Int64, Symbol}}:
 (41, :y_1)
 (82, :y_2)
 (123, :y_3)
```
"""
macro named(expr)
    name, call = split_assign(expr)
    if Meta.isexpr(name, :ref)
        name, idxs = name.args
        check_name(name)
        esc(_named_idxs(name, idxs, :($(gensym()) -> $call)))
    else
        check_name(name)
        esc(:($name = $(_named(name, call))))
    end
end

macro named(name::Symbol, idxs, call)
    esc(_named_idxs(name, idxs, call))
end

function _config(expr, namespace)
    cn = Base.Fix2(_config, namespace)
    if Meta.isexpr(expr, :.)
        return :($getproperty($(map(cn, expr.args)...); namespace = $namespace))
    elseif Meta.isexpr(expr, :function)
        def = splitdef(expr)
        def[:args] = map(cn, def[:args])
        def[:body] = cn(def[:body])
        combinedef(def)
    elseif expr isa Expr && !isempty(expr.args)
        return Expr(expr.head, map(cn, expr.args)...)
    elseif Meta.isexpr(expr, :(=))
        return Expr(:(=), map(cn, expr.args)...)
    else
        expr
    end
end

"""
$(SIGNATURES)

Rewrite `@nonamespace a.b.c` to
`getvar(getvar(a, :b; namespace = false), :c; namespace = false)`.

This is the default behavior of `getvar`. This should be used when inheriting states from a model.
"""
macro nonamespace(expr)
    esc(_config(expr, false))
end

"""
$(SIGNATURES)

Rewrite `@namespace a.b.c` to
`getvar(getvar(a, :b; namespace = true), :c; namespace = true)`.
"""
macro namespace(expr)
    esc(_config(expr, true))
end

"""
$(SIGNATURES)

Structurally simplify algebraic equations in a system and compute the
topological sort of the observed equations. When `simplify=true`, the `simplify`
function will be applied during the tearing process. It also takes kwargs
`allow_symbolic=false` and `allow_parameter=true` which limits the coefficient
types during tearing.
"""
function structural_simplify(sys::AbstractSystem; simplify = false, kwargs...)
    sys = expand_connections(sys)
    sys = alias_elimination(sys)
    state = TearingState(sys)
    state = inputs_to_parameters!(state)
    sys = state.sys
    check_consistency(state)
    if sys isa ODESystem
        sys = dae_order_lowering(dummy_derivative(sys, state))
    end
    state = TearingState(sys)
    find_solvables!(state; kwargs...)
    sys = tearing_reassemble(state, tearing(state), simplify = simplify)
    fullstates = [map(eq -> eq.lhs, observed(sys)); states(sys)]
    @set! sys.observed = topsort_equations(observed(sys), fullstates)
    invalidate_cache!(sys)
    return sys
end

@latexrecipe function f(sys::AbstractSystem)
    return latexify(equations(sys))
end

function Base.show(io::IO, ::MIME"text/latex", x::AbstractSystem)
    print(io, "\$\$ " * latexify(x) * " \$\$")
end

struct InvalidSystemException <: Exception
    msg::String
end
function Base.showerror(io::IO, e::InvalidSystemException)
    print(io, "InvalidSystemException: ", e.msg)
end

struct ExtraVariablesSystemException <: Exception
    msg::String
end
function Base.showerror(io::IO, e::ExtraVariablesSystemException)
    print(io, "ExtraVariablesSystemException: ", e.msg)
end

struct ExtraEquationsSystemException <: Exception
    msg::String
end
function Base.showerror(io::IO, e::ExtraEquationsSystemException)
    print(io, "ExtraEquationsSystemException: ", e.msg)
end

function AbstractTrees.children(sys::ModelingToolkit.AbstractSystem)
    ModelingToolkit.get_systems(sys)
end
function AbstractTrees.printnode(io::IO, sys::ModelingToolkit.AbstractSystem)
    print(io, nameof(sys))
end
function Base.IteratorEltype(::Type{<:TreeIterator{ModelingToolkit.AbstractSystem}})
    Base.HasEltype()
end
function Base.eltype(::Type{<:TreeIterator{ModelingToolkit.AbstractSystem}})
    ModelingToolkit.AbstractSystem
end

function check_eqs_u0(eqs, dvs, u0; check_length = true, kwargs...)
    if u0 !== nothing
        if check_length
            if !(length(eqs) == length(dvs) == length(u0))
                throw(ArgumentError("Equations ($(length(eqs))), states ($(length(dvs))), and initial conditions ($(length(u0))) are of different lengths. To allow a different number of equations than states use kwarg check_length=false."))
            end
        elseif length(dvs) != length(u0)
            throw(ArgumentError("States ($(length(dvs))) and initial conditions ($(length(u0))) are of different lengths."))
        end
    elseif check_length && (length(eqs) != length(dvs))
        throw(ArgumentError("Equations ($(length(eqs))) and states ($(length(dvs))) are of different lengths. To allow these to differ use kwarg check_length=false."))
    end
    return nothing
end

###
### Inheritance & composition
###
function Base.hash(sys::AbstractSystem, s::UInt)
    s = hash(nameof(sys), s)
    s = foldr(hash, get_systems(sys), init = s)
    s = foldr(hash, get_states(sys), init = s)
    s = foldr(hash, get_ps(sys), init = s)
    if sys isa OptimizationSystem
        s = hash(get_op(sys), s)
    else
        s = foldr(hash, get_eqs(sys), init = s)
    end
    s = foldr(hash, get_observed(sys), init = s)
    s = foldr(hash, get_continuous_events(sys), init = s)
    s = hash(independent_variables(sys), s)
    return s
end

"""
$(TYPEDSIGNATURES)

extend the `basesys` with `sys`, the resulting system would inherit `sys`'s name
by default.
"""
function extend(sys::AbstractSystem, basesys::AbstractSystem; name::Symbol = nameof(sys))
    T = SciMLBase.parameterless_type(basesys)
    ivs = independent_variables(basesys)
    if !(sys isa T)
        if length(ivs) == 0
            sys = convert_system(T, sys)
        elseif length(ivs) == 1
            sys = convert_system(T, sys, ivs[1])
        else
            throw("Extending multivariate systems is not supported")
        end
    end

    eqs = union(get_eqs(basesys), get_eqs(sys))
    sts = union(get_states(basesys), get_states(sys))
    ps = union(get_ps(basesys), get_ps(sys))
    obs = union(get_observed(basesys), get_observed(sys))
    evs = union(get_continuous_events(basesys), get_continuous_events(sys))
    defs = merge(get_defaults(basesys), get_defaults(sys)) # prefer `sys`
    syss = union(get_systems(basesys), get_systems(sys))

    if length(ivs) == 0
        T(eqs, sts, ps, observed = obs, defaults = defs, name = name, systems = syss,
          continuous_events = evs)
    elseif length(ivs) == 1
        T(eqs, ivs[1], sts, ps, observed = obs, defaults = defs, name = name,
          systems = syss, continuous_events = evs)
    end
end

function Base.:(&)(sys::AbstractSystem, basesys::AbstractSystem; name::Symbol = nameof(sys))
    extend(sys, basesys; name = name)
end

"""
$(SIGNATURES)

compose multiple systems together. The resulting system would inherit the first
system's name.
"""
function compose(sys::AbstractSystem, systems::AbstractArray; name = nameof(sys))
    nsys = length(systems)
    nsys == 0 && return sys
    @set! sys.name = name
    @set! sys.systems = [get_systems(sys); systems]
    return sys
end
function compose(syss...; name = nameof(first(syss)))
    compose(first(syss), collect(syss[2:end]); name = name)
end
Base.:(∘)(sys1::AbstractSystem, sys2::AbstractSystem) = compose(sys1, sys2)

function UnPack.unpack(sys::ModelingToolkit.AbstractSystem, ::Val{p}) where {p}
    getproperty(sys, p; namespace = false)
end
