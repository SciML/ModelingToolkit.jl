const SYSTEM_COUNT = Threads.Atomic{UInt}(0)

get_component_type(x::AbstractSystem) = get_gui_metadata(x).type
struct GUIMetadata
    type::GlobalRef
    layout::Any
end

GUIMetadata(type) = GUIMetadata(type, nothing)

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

Calculate the Jacobian matrix of a system.

Returns a matrix of [`Num`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_jacobian end

"""
```julia
calculate_control_jacobian(sys::AbstractSystem)
```

Calculate the Jacobian matrix of a system with respect to the system's controls.

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
generate_tgrad(sys::AbstractTimeDependentSystem, dvs = states(sys), ps = parameters(sys),
               expression = Val{true}; kwargs...)
```

Generates a function for the time gradient of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_tgrad end

"""
```julia
generate_gradient(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys),
                  expression = Val{true}; kwargs...)
```

Generates a function for the gradient of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_gradient end

"""
```julia
generate_jacobian(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys),
                  expression = Val{true}; sparse = false, kwargs...)
```

Generates a function for the Jacobian matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_jacobian end

"""
```julia
generate_factorized_W(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys),
                      expression = Val{true}; sparse = false, kwargs...)
```

Generates a function for the factorized W matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_factorized_W end

"""
```julia
generate_hessian(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys),
                 expression = Val{true}; sparse = false, kwargs...)
```

Generates a function for the hessian matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_hessian end

"""
```julia
generate_function(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys),
                  expression = Val{true}; kwargs...)
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
function SymbolicIndexingInterface.independent_variables(sys::AbstractSystem)
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

function SymbolicIndexingInterface.independent_variables(sys::AbstractTimeDependentSystem)
    [getfield(sys, :iv)]
end
SymbolicIndexingInterface.independent_variables(sys::AbstractTimeIndependentSystem) = []
function SymbolicIndexingInterface.independent_variables(sys::AbstractMultivariateSystem)
    getfield(sys, :ivs)
end

iscomplete(sys::AbstractSystem) = isdefined(sys, :complete) && getfield(sys, :complete)

"""
$(TYPEDSIGNATURES)

Mark a system as completed. If a system is complete, the system will no longer
namespace its subsystems or variables, i.e. `isequal(complete(sys).v.i, v.i)`.
"""
function complete(sys::AbstractSystem)
    isdefined(sys, :complete) ? (@set! sys.complete = true) : sys
end

for prop in [:eqs
    :tag
    :noiseeqs
    :iv
    :states
    :ps
    :tspan
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
    :constraints
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
    :substitutions
    :metadata
    :gui_metadata
    :discrete_subsystems
    :unknown_states
    :split_idxs
    :parent]
    fname1 = Symbol(:get_, prop)
    fname2 = Symbol(:has_, prop)
    @eval begin
        $fname1(sys::AbstractSystem) = getfield(sys, $(QuoteNode(prop)))
        $fname2(sys::AbstractSystem) = isdefined(sys, $(QuoteNode(prop)))
    end
end

has_equations(::AbstractSystem) = true

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

rename(x, name) = @set x.name = name

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

function Base.getproperty(sys::AbstractSystem, name::Symbol; namespace = !iscomplete(sys))
    if has_parent(sys) && (parent = get_parent(sys); parent !== nothing)
        sys = parent
    end
    wrap(getvar(sys, name; namespace = namespace))
end
function getvar(sys::AbstractSystem, name::Symbol; namespace = !iscomplete(sys))
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
    if i !== nothing
        return namespace ? renamespace(sys, sts[i]) : sts[i]
    end

    if has_observed(sys)
        obs = get_observed(sys)
        i = findfirst(x -> getname(x.lhs) == name, obs)
        if i !== nothing
            return namespace ? renamespace(sys, obs[i].lhs) : obs[i].lhs
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

apply_to_variables(f::F, ex) where {F} = _apply_to_variables(f, ex)
apply_to_variables(f::F, ex::Num) where {F} = wrap(_apply_to_variables(f, unwrap(ex)))
function _apply_to_variables(f::F, ex) where {F}
    if isvariable(ex)
        return f(ex)
    end
    istree(ex) || return ex
    similarterm(ex, _apply_to_variables(f, operation(ex)),
        map(Base.Fix1(_apply_to_variables, f), arguments(ex)),
        metadata = metadata(ex))
end

abstract type SymScope end

struct LocalScope <: SymScope end
function LocalScope(sym::Union{Num, Symbolic})
    apply_to_variables(sym) do sym
        setmetadata(sym, SymScope, LocalScope())
    end
end

struct ParentScope <: SymScope
    parent::SymScope
end
function ParentScope(sym::Union{Num, Symbolic})
    apply_to_variables(sym) do sym
        setmetadata(sym, SymScope,
            ParentScope(getmetadata(value(sym), SymScope, LocalScope())))
    end
end

struct DelayParentScope <: SymScope
    parent::SymScope
    N::Int
end
function DelayParentScope(sym::Union{Num, Symbolic}, N)
    apply_to_variables(sym) do sym
        setmetadata(sym, SymScope,
            DelayParentScope(getmetadata(value(sym), SymScope, LocalScope()), N))
    end
end
DelayParentScope(sym::Union{Num, Symbolic}) = DelayParentScope(sym, 1)

struct GlobalScope <: SymScope end
function GlobalScope(sym::Union{Num, Symbolic})
    apply_to_variables(sym) do sym
        setmetadata(sym, SymScope, GlobalScope())
    end
end

renamespace(sys, eq::Equation) = namespace_equation(eq, sys)

renamespace(names::AbstractVector, x) = foldr(renamespace, names, init = x)
function renamespace(sys, x)
    sys === nothing && return x
    x = unwrap(x)
    if x isa Symbolic
        T = typeof(x)
        if istree(x) && operation(x) isa Operator
            return similarterm(x, operation(x),
                Any[renamespace(sys, only(arguments(x)))])::T
        end
        let scope = getmetadata(x, SymScope, LocalScope())
            if scope isa LocalScope
                rename(x, renamespace(getname(sys), getname(x)))::T
            elseif scope isa ParentScope
                setmetadata(x, SymScope, scope.parent)::T
            elseif scope isa DelayParentScope
                if scope.N > 0
                    x = setmetadata(x, SymScope,
                        DelayParentScope(scope.parent, scope.N - 1))
                    rename(x, renamespace(getname(sys), getname(x)))::T
                else
                    #rename(x, renamespace(getname(sys), getname(x)))
                    setmetadata(x, SymScope, scope.parent)::T
                end
            else # GlobalScope
                x::T
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
    Dict((isparameter(k) ? parameters(sys, k) : states(sys, k)) => namespace_expr(v, sys)
         for (k, v) in pairs(defs))
end

function namespace_equations(sys::AbstractSystem, ivs = independent_variables(sys))
    eqs = equations(sys)
    isempty(eqs) && return Equation[]
    map(eq -> namespace_equation(eq, sys; ivs), eqs)
end

function namespace_equation(eq::Equation,
        sys,
        n = nameof(sys);
        ivs = independent_variables(sys))
    _lhs = namespace_expr(eq.lhs, sys, n; ivs)
    _rhs = namespace_expr(eq.rhs, sys, n; ivs)
    _lhs ~ _rhs
end

function namespace_assignment(eq::Assignment, sys)
    _lhs = namespace_expr(eq.lhs, sys)
    _rhs = namespace_expr(eq.rhs, sys)
    Assignment(_lhs, _rhs)
end

function namespace_expr(O, sys, n = nameof(sys); ivs = independent_variables(sys))
    O = unwrap(O)
    if any(isequal(O), ivs)
        return O
    elseif istree(O)
        T = typeof(O)
        renamed = let sys = sys, n = n, T = T
            map(a -> namespace_expr(a, sys, n; ivs)::Any, arguments(O))
        end
        if isvariable(O)
            # Use renamespace so the scope is correct, and make sure to use the
            # metadata from the rescoped variable
            rescoped = renamespace(n, O)
            similarterm(O, operation(rescoped), renamed,
                metadata = metadata(rescoped))::T
        else
            similarterm(O, operation(O), renamed, metadata = metadata(O))::T
        end
    elseif isvariable(O)
        renamespace(n, O)
    elseif O isa Array
        let sys = sys, n = n
            map(o -> namespace_expr(o, sys, n; ivs), O)
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
           [sts; reduce(vcat, namespace_variables.(systems))])
end

function SymbolicIndexingInterface.parameters(sys::AbstractSystem)
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

states(sys::Union{AbstractSystem, Nothing}, v) = renamespace(sys, v)
parameters(sys::Union{AbstractSystem, Nothing}, v) = toparam(states(sys, v))
for f in [:states, :parameters]
    @eval function $f(sys::AbstractSystem, vs::AbstractArray)
        map(v -> $f(sys, v), vs)
    end
end

flatten(sys::AbstractSystem, args...) = sys

function equations(sys::AbstractSystem)
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

function preface(sys::AbstractSystem)
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

function time_varying_as_func(x, sys::AbstractTimeDependentSystem)
    # if something is not x(t) (the current state)
    # but is `x(t-1)` or something like that, pass in `x` as a callable function rather
    # than pass in a value in place of x(t).
    #
    # This is done by just making `x` the argument of the function.
    if istree(x) &&
       issym(operation(x)) &&
       !(length(arguments(x)) == 1 && isequal(arguments(x)[1], get_iv(sys)))
        return operation(x)
    end
    return x
end

SymbolicIndexingInterface.is_indep_sym(sys::AbstractSystem, sym) = isequal(sym, get_iv(sys))

"""
$(SIGNATURES)

Return a list of actual states needed to be solved by solvers.
"""
function unknown_states(sys::AbstractSystem)
    sts = states(sys)
    if has_unknown_states(sys)
        sts = something(get_unknown_states(sys), sts)
    end
    return sts
end

function SymbolicIndexingInterface.state_sym_to_index(sys::AbstractSystem, sym)
    findfirst(isequal(sym), unknown_states(sys))
end
function SymbolicIndexingInterface.is_state_sym(sys::AbstractSystem, sym)
    !isnothing(SymbolicIndexingInterface.state_sym_to_index(sys, sym))
end

function SymbolicIndexingInterface.param_sym_to_index(sys::AbstractSystem, sym)
    findfirst(isequal(sym), SymbolicIndexingInterface.parameters(sys))
end
function SymbolicIndexingInterface.is_param_sym(sys::AbstractSystem, sym)
    !isnothing(SymbolicIndexingInterface.param_sym_to_index(sys, sym))
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
    issym(t) && return nameof(t)
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
    obs = observed(sys)
    obsvars = [o.lhs for o in obs]
    obsvarsname = gensym(:obs)
    push_vars!(stmt, obsvarsname, Symbol("@variables"), obsvars)

    var2name = Dict{Any, Symbol}()
    for v in Iterators.flatten((sts, ps, obsvars))
        var2name[v] = getname(v)
    end

    eqs_name = push_eqs!(stmt, full_equations(sys), var2name)
    defs_name = push_defaults!(stmt, defaults(sys), var2name)
    obs_name = push_eqs!(stmt, obs, var2name)

    if sys isa ODESystem
        iv = get_iv(sys)
        ivname = gensym(:iv)
        push!(stmt, :($ivname = (@variables $(getname(iv)))[1]))
        push!(stmt,
            :($ODESystem($eqs_name, $ivname, $stsname, $psname; defaults = $defs_name,
                observed = $obs_name,
                name = $name, checks = false)))
    elseif sys isa NonlinearSystem
        push!(stmt,
            :($NonlinearSystem($eqs_name, $stsname, $psname; defaults = $defs_name,
                observed = $obs_name,
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
    sys, (csets, _) = generate_connection_set(sys)
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
            description = getdescription(s)
            if description !== nothing && description != ""
                print(io, ": ", description)
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
            description = getdescription(s)
            if description !== nothing && description != ""
                print(io, ": ", description)
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

varname_fix!(s) = return

function varname_fix!(expr::Expr)
    for arg in expr.args
        MLStyle.@match arg begin
            ::Symbol => continue
            Expr(:kw, a...) || Expr(:kw, a) => varname_sanitization!(arg)
            Expr(:parameters, a...) => begin
                for _arg in arg.args
                    varname_sanitization!(_arg)
                end
            end
            _ => @debug "skipping variable sanitization of $arg"
        end
    end
end

varname_sanitization!(a) = return

function varname_sanitization!(expr::Expr)
    var_splits = split(string(expr.args[1]), ".")
    if length(var_splits) > 1
        expr.args[1] = Symbol(join(var_splits, "__"))
    end
end

function _named(name, call, runtime = false)
    has_kw = false
    call isa Expr || throw(Meta.ParseError("The rhs must be an Expr. Got $call."))
    if length(call.args) >= 2 && call.args[2] isa Expr
        # canonicalize to use `:parameters`
        if call.args[2].head === :kw
            call = Expr(call.head, call.args[1], Expr(:parameters, call.args[2:end]...))
            has_kw = true
        elseif call.args[2].head === :parameters
            has_kw = true
        end
    end

    varname_fix!(call)

    if !has_kw
        param = Expr(:parameters)
        if length(call.args) == 1
            push!(call.args, param)
        else
            insert!(call.args, 2, param)
        end
    end

    is_sys_construction = gensym("###__is_system_construction###")
    kws = call.args[2].args
    for (i, kw) in enumerate(kws)
        if Meta.isexpr(kw, (:(=), :kw))
            kw.args[2] = :($is_sys_construction ? $(kw.args[2]) :
                           $default_to_parentscope($(kw.args[2])))
        elseif kw isa Symbol
            rhs = :($is_sys_construction ? $(kw) : $default_to_parentscope($(kw)))
            kws[i] = Expr(:kw, kw, rhs)
        end
    end

    if !any(kw -> (kw isa Symbol ? kw : kw.args[1]) == :name, kws) # don't overwrite `name` kwarg
        pushfirst!(kws, Expr(:kw, :name, runtime ? name : Meta.quot(name)))
    end
    op = call.args[1]
    quote
        $is_sys_construction = ($op isa $DataType) && ($op <: $AbstractSystem)
        $call
    end
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

function single_named_expr(expr)
    name, call = split_assign(expr)
    if Meta.isexpr(name, :ref)
        name, idxs = name.args
        check_name(name)
        var = gensym(name)
        ex = quote
            $var = $(_named(name, call))
            $name = map(i -> $rename($var, Symbol($(Meta.quot(name)), :_, i)), $idxs)
        end
        ex
    else
        check_name(name)
        :($name = $(_named(name, call)))
    end
end

function named_expr(expr)
    if Meta.isexpr(expr, :block)
        newexpr = Expr(:block)
        names = Expr(:vcat)
        for ex in expr.args
            ex isa LineNumberNode && continue
            push!(newexpr.args, single_named_expr(ex))
            push!(names.args, ex.args[1])
        end
        push!(newexpr.args, names)
        newexpr
    else
        single_named_expr(expr)
    end
end

function check_name(name)
    name isa Symbol ||
        throw(Meta.ParseError("The lhs must be a symbol (a) or a ref (a[1:10]). Got $name."))
end

"""
    @named y = foo(x)
    @named y[1:10] = foo(x)
    @named begin
        y[1:10] = foo(x)
        z = foo(x)
    end # returns `[y; z]`
    @named y 1:10 i -> foo(x*i)  # This is not recommended

Pass the LHS name to the model. When it's calling anything that's not an
AbstractSystem, it wraps all keyword arguments in `default_to_parentscope` so
that namespacing works intuitively when passing a symbolic default into a
component.

Examples:

```julia-repl
julia> using ModelingToolkit

julia> foo(i; name) = (; i, name)
foo (generic function with 1 method)

julia> x = 41
41

julia> @named y = foo(x)
(i = 41, name = :y)

julia> @named y[1:3] = foo(x)
3-element Vector{NamedTuple{(:i, :name), Tuple{Int64, Symbol}}}:
 (i = 41, name = :y_1)
 (i = 41, name = :y_2)
 (i = 41, name = :y_3)
```
"""
macro named(expr)
    esc(named_expr(expr))
end

macro named(name::Symbol, idxs, call)
    esc(_named_idxs(name, idxs, call))
end

function default_to_parentscope(v)
    uv = unwrap(v)
    uv isa Symbolic || return v
    apply_to_variables(v) do sym
        if !hasmetadata(uv, SymScope)
            setmetadata(sym, SymScope,
                ParentScope(getmetadata(value(sym), SymScope, LocalScope())))
        else
            sym
        end
    end
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

function component_post_processing(expr, isconnector)
    @assert expr isa Expr && (expr.head == :function || (expr.head == :(=) &&
              expr.args[1] isa Expr &&
              expr.args[1].head == :call))

    sig = expr.args[1]
    body = expr.args[2]

    fname = sig.args[1]
    args = sig.args[2:end]

    quote
        function $fname($(args...))
            # we need to create a closure to escape explicit return in `body`.
            res = (() -> $body)()
            if $isdefined(res, :gui_metadata) && $getfield(res, :gui_metadata) === nothing
                name = $(Meta.quot(fname))
                if $isconnector
                    $Setfield.@set!(res.connector_type=$connector_type(res))
                end
                $Setfield.@set!(res.gui_metadata=$GUIMetadata($GlobalRef(@__MODULE__, name)))
            else
                res
            end
        end
    end
end

macro component(expr)
    esc(component_post_processing(expr, false))
end

macro mtkbuild(expr)
    named_expr = ModelingToolkit.named_expr(expr)
    name = named_expr.args[1]
    esc(quote
        $named_expr
        $name = $structural_simplify($name)
    end)
end

"""
$(SIGNATURES)

Replace functions with singularities with a function that errors with symbolic
information. E.g.

```julia-repl
julia> sys = debug_system(sys);

julia> prob = ODEProblem(sys, [], (0, 1.0));

julia> du = zero(prob.u0);

julia> prob.f(du, prob.u0, prob.p, 0.0)
ERROR: DomainError with (-1.0,):
log errors with input(s): -cos(Q(t)) => -1.0
Stacktrace:
  [1] (::ModelingToolkit.LoggedFun{typeof(log)})(args::Float64)
  ...
```
"""
function debug_system(sys::AbstractSystem)
    if has_systems(sys) && !isempty(get_systems(sys))
        error("debug_system only works on systems with no sub-systems!")
    end
    if has_eqs(sys)
        @set! sys.eqs = debug_sub.(equations(sys))
    end
    if has_observed(sys)
        @set! sys.observed = debug_sub.(observed(sys))
    end
    return sys
end

function eliminate_constants(sys::AbstractSystem)
    if has_eqs(sys)
        eqs = get_eqs(sys)
        eq_cs = collect_constants(eqs)
        if !isempty(eq_cs)
            new_eqs = eliminate_constants(eqs, eq_cs)
            @set! sys.eqs = new_eqs
        end
    end
    return sys
end

function io_preprocessing(sys::AbstractSystem, inputs,
        outputs; simplify = false, kwargs...)
    sys, input_idxs = structural_simplify(sys, (inputs, outputs); simplify, kwargs...)

    eqs = equations(sys)
    alg_start_idx = findfirst(!isdiffeq, eqs)
    if alg_start_idx === nothing
        alg_start_idx = length(eqs) + 1
    end
    diff_idxs = 1:(alg_start_idx - 1)
    alge_idxs = alg_start_idx:length(eqs)

    sys, diff_idxs, alge_idxs, input_idxs
end

"""
    lin_fun, simplified_sys = linearization_function(sys::AbstractSystem, inputs, outputs; simplify = false, initialize = true, kwargs...)

Return a function that linearizes the system `sys`. The function [`linearize`](@ref) provides a higher-level and easier to use interface.

`lin_fun` is a function `(variables, p, t) -> (; f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u)`, i.e., it returns a NamedTuple with the Jacobians of `f,g,h` for the nonlinear `sys` (technically for `simplified_sys`) on the form

```math
\\begin{aligned}
ẋ &= f(x, z, u) \\\\
0 &= g(x, z, u) \\\\
y &= h(x, z, u)
\\end{aligned}
```

where `x` are differential state variables, `z` algebraic variables, `u` inputs and `y` outputs. To obtain a linear statespace representation, see [`linearize`](@ref). The input argument `variables` is a vector defining the operating point, corresponding to `states(simplified_sys)` and `p` is a vector corresponding to the parameters of `simplified_sys`. Note: all variables in `inputs` have been converted to parameters in `simplified_sys`.

The `simplified_sys` has undergone [`structural_simplify`](@ref) and had any occurring input or output variables replaced with the variables provided in arguments `inputs` and `outputs`. The states of this system also indicate the order of the states that holds for the linearized matrices.

# Arguments:

  - `sys`: An [`ODESystem`](@ref). This function will automatically apply simplification passes on `sys` and return the resulting `simplified_sys`.
  - `inputs`: A vector of variables that indicate the inputs of the linearized input-output model.
  - `outputs`: A vector of variables that indicate the outputs of the linearized input-output model.
  - `simplify`: Apply simplification in tearing.
  - `initialize`: If true, a check is performed to ensure that the operating point is consistent (satisfies algebraic equations). If the op is not consistent, initialization is performed.
  - `kwargs`: Are passed on to `find_solvables!`

See also [`linearize`](@ref) which provides a higher-level interface.
"""
function linearization_function(sys::AbstractSystem, inputs,
        outputs; simplify = false,
        initialize = true,
        op = Dict(),
        p = DiffEqBase.NullParameters(),
        zero_dummy_der = false,
        kwargs...)
    inputs isa AbstractVector || (inputs = [inputs])
    outputs isa AbstractVector || (outputs = [outputs])
    ssys, diff_idxs, alge_idxs, input_idxs = io_preprocessing(sys, inputs, outputs;
        simplify,
        kwargs...)
    if zero_dummy_der
        dummyder = setdiff(states(ssys), states(sys))
        defs = Dict(x => 0.0 for x in dummyder)
        @set! ssys.defaults = merge(defs, defaults(ssys))
        op = merge(defs, op)
    end
    sys = ssys
    x0 = merge(defaults(sys), Dict(missing_variable_defaults(sys)), op)
    u0, p, _ = get_u0_p(sys, x0, p; use_union = false, tofloat = true)
    p, split_idxs = split_parameters_by_type(p)
    ps = parameters(sys)
    if p isa Tuple
        ps = Base.Fix1(getindex, ps).(split_idxs)
        ps = (ps...,) #if p is Tuple, ps should be Tuple
    end

    lin_fun = let diff_idxs = diff_idxs,
        alge_idxs = alge_idxs,
        input_idxs = input_idxs,
        sts = states(sys),
        fun = ODEFunction{true, SciMLBase.FullSpecialize}(sys, states(sys), ps; p = p),
        h = build_explicit_observed_function(sys, outputs),
        chunk = ForwardDiff.Chunk(input_idxs)

        function (u, p, t)
            if u !== nothing # Handle systems without states
                length(sts) == length(u) ||
                    error("Number of state variables ($(length(sts))) does not match the number of input states ($(length(u)))")
                if initialize && !isempty(alge_idxs) # This is expensive and can be omitted if the user knows that the system is already initialized
                    residual = fun(u, p, t)
                    if norm(residual[alge_idxs]) > √(eps(eltype(residual)))
                        prob = ODEProblem(fun, u, (t, t + 1), p)
                        integ = init(prob, OrdinaryDiffEq.Rodas5P())
                        u = integ.u
                    end
                end
                uf = SciMLBase.UJacobianWrapper(fun, t, p)
                fg_xz = ForwardDiff.jacobian(uf, u)
                h_xz = ForwardDiff.jacobian(let p = p, t = t
                        xz -> h(xz, p, t)
                    end, u)
                pf = SciMLBase.ParamJacobianWrapper(fun, t, u)
                fg_u = jacobian_wrt_vars(pf, p, input_idxs, chunk)
            else
                length(sts) == 0 ||
                    error("Number of state variables (0) does not match the number of input states ($(length(u)))")
                fg_xz = zeros(0, 0)
                h_xz = fg_u = zeros(0, length(inputs))
            end
            hp = let u = u, t = t
                p -> h(u, p, t)
            end
            h_u = jacobian_wrt_vars(hp, p, input_idxs, chunk)
            (f_x = fg_xz[diff_idxs, diff_idxs],
                f_z = fg_xz[diff_idxs, alge_idxs],
                g_x = fg_xz[alge_idxs, diff_idxs],
                g_z = fg_xz[alge_idxs, alge_idxs],
                f_u = fg_u[diff_idxs, :],
                g_u = fg_u[alge_idxs, :],
                h_x = h_xz[:, diff_idxs],
                h_z = h_xz[:, alge_idxs],
                h_u = h_u)
        end
    end
    return lin_fun, sys
end

"""
    (; A, B, C, D), simplified_sys = linearize_symbolic(sys::AbstractSystem, inputs, outputs; simplify = false, allow_input_derivatives = false, kwargs...)

Similar to [`linearize`](@ref), but returns symbolic matrices `A,B,C,D` rather than numeric. While `linearize` uses ForwardDiff to perform the linearization, this function uses `Symbolics.jacobian`.

See [`linearize`](@ref) for a description of the arguments.

# Extended help
The named tuple returned as the first argument additionally contains the jacobians `f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u` of
```math
\\begin{aligned}
ẋ &= f(x, z, u) \\\\
0 &= g(x, z, u) \\\\
y &= h(x, z, u)
\\end{aligned}
```
where `x` are differential state variables, `z` algebraic variables, `u` inputs and `y` outputs.
"""
function linearize_symbolic(sys::AbstractSystem, inputs,
        outputs; simplify = false, allow_input_derivatives = false,
        kwargs...)
    sys, diff_idxs, alge_idxs, input_idxs = io_preprocessing(sys, inputs, outputs; simplify,
        kwargs...)
    sts = states(sys)
    t = get_iv(sys)
    p = parameters(sys)

    fun = generate_function(sys, sts, p; expression = Val{false})[1]
    dx = fun(sts, p, t)

    h = build_explicit_observed_function(sys, outputs)
    y = h(sts, p, t)

    fg_xz = Symbolics.jacobian(dx, sts)
    fg_u = Symbolics.jacobian(dx, inputs)
    h_xz = Symbolics.jacobian(y, sts)
    h_u = Symbolics.jacobian(y, inputs)
    f_x = fg_xz[diff_idxs, diff_idxs]
    f_z = fg_xz[diff_idxs, alge_idxs]
    g_x = fg_xz[alge_idxs, diff_idxs]
    g_z = fg_xz[alge_idxs, alge_idxs]
    f_u = fg_u[diff_idxs, :]
    g_u = fg_u[alge_idxs, :]
    h_x = h_xz[:, diff_idxs]
    h_z = h_xz[:, alge_idxs]

    nx, nu = size(f_u)
    nz = size(f_z, 2)
    ny = size(h_x, 1)

    D = h_u

    if isempty(g_z) # ODE
        A = f_x
        B = f_u
        C = h_x
    else
        gz = lu(g_z; check = false)
        issuccess(gz) ||
            error("g_z not invertible, this indicates that the DAE is of index > 1.")
        gzgx = -(gz \ g_x)
        A = [f_x f_z
            gzgx*f_x gzgx*f_z]
        B = [f_u
            gzgx * f_u] # The cited paper has zeros in the bottom block, see derivation in https://github.com/SciML/ModelingToolkit.jl/pull/1691 for the correct formula

        C = [h_x h_z]
        Bs = -(gz \ g_u) # This equation differ from the cited paper, the paper is likely wrong since their equaiton leads to a dimension mismatch.
        if !iszero(Bs)
            if !allow_input_derivatives
                der_inds = findall(vec(any(!iszero, Bs, dims = 1)))
                @show typeof(der_inds)
                error("Input derivatives appeared in expressions (-g_z\\g_u != 0), the following inputs appeared differentiated: $(ModelingToolkit.inputs(sys)[der_inds]). Call `linearize_symbolic` with keyword argument `allow_input_derivatives = true` to allow this and have the returned `B` matrix be of double width ($(2nu)), where the last $nu inputs are the derivatives of the first $nu inputs.")
            end
            B = [B [zeros(nx, nu); Bs]]
            D = [D zeros(ny, nu)]
        end
    end

    (; A, B, C, D, f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u), sys
end

function markio!(state, orig_inputs, inputs, outputs; check = true)
    fullvars = get_fullvars(state)
    inputset = Dict{Any, Bool}(i => false for i in inputs)
    outputset = Dict{Any, Bool}(o => false for o in outputs)
    for (i, v) in enumerate(fullvars)
        if v in keys(inputset)
            if v in keys(outputset)
                v = setio(v, true, true)
                outputset[v] = true
            else
                v = setio(v, true, false)
            end
            inputset[v] = true
            fullvars[i] = v
        elseif v in keys(outputset)
            v = setio(v, false, true)
            outputset[v] = true
            fullvars[i] = v
        else
            if isinput(v)
                push!(orig_inputs, v)
            end
            v = setio(v, false, false)
            fullvars[i] = v
        end
    end
    if check
        ikeys = keys(filter(!last, inputset))
        if !isempty(ikeys)
            error("Some specified inputs were not found in system. The following variables were not found ",
                ikeys)
        end
    end
    check && (all(values(outputset)) ||
     error("Some specified outputs were not found in system. The following Dict indicates the found variables ",
        outputset))
    state, orig_inputs
end

"""
    (; A, B, C, D), simplified_sys = linearize(sys, inputs, outputs;    t=0.0, op = Dict(), allow_input_derivatives = false, zero_dummy_der=false, kwargs...)
    (; A, B, C, D)                 = linearize(simplified_sys, lin_fun; t=0.0, op = Dict(), allow_input_derivatives = false, zero_dummy_der=false)

Linearize `sys` between `inputs` and `outputs`, both vectors of variables. Return a NamedTuple with the matrices of a linear statespace representation
on the form

```math
\\begin{aligned}
ẋ &= Ax + Bu\\\\
y &= Cx + Du
\\end{aligned}
```

The first signature automatically calls [`linearization_function`](@ref) internally,
while the second signature expects the outputs of [`linearization_function`](@ref) as input.

`op` denotes the operating point around which to linearize. If none is provided,
the default values of `sys` are used.

If `allow_input_derivatives = false`, an error will be thrown if input derivatives (``u̇``) appear as inputs in the linearized equations. If input derivatives are allowed, the returned `B` matrix will be of double width, corresponding to the input `[u; u̇]`.

`zero_dummy_der` can be set to automatically set the operating point to zero for all dummy derivatives.

See also [`linearization_function`](@ref) which provides a lower-level interface, [`linearize_symbolic`](@ref) and [`ModelingToolkit.reorder_states`](@ref).

See extended help for an example.

The implementation and notation follows that of
["Linear Analysis Approach for Modelica Models", Allain et al. 2009](https://ep.liu.se/ecp/043/075/ecp09430097.pdf)

# Extended help

This example builds the following feedback interconnection and linearizes it from the input of `F` to the output of `P`.

```

  r ┌─────┐       ┌─────┐     ┌─────┐
───►│     ├──────►│     │  u  │     │
    │  F  │       │  C  ├────►│  P  │ y
    └─────┘     ┌►│     │     │     ├─┬─►
                │ └─────┘     └─────┘ │
                │                     │
                └─────────────────────┘
```

```julia
using ModelingToolkit
@variables t
function plant(; name)
    @variables x(t) = 1
    @variables u(t)=0 y(t)=0
    D = Differential(t)
    eqs = [D(x) ~ -x + u
           y ~ x]
    ODESystem(eqs, t; name = name)
end

function ref_filt(; name)
    @variables x(t)=0 y(t)=0
    @variables u(t)=0 [input = true]
    D = Differential(t)
    eqs = [D(x) ~ -2 * x + u
           y ~ x]
    ODESystem(eqs, t, name = name)
end

function controller(kp; name)
    @variables y(t)=0 r(t)=0 u(t)=0
    @parameters kp = kp
    eqs = [
        u ~ kp * (r - y),
    ]
    ODESystem(eqs, t; name = name)
end

@named f = ref_filt()
@named c = controller(1)
@named p = plant()

connections = [f.y ~ c.r # filtered reference to controller reference
               c.u ~ p.u # controller output to plant input
               p.y ~ c.y]

@named cl = ODESystem(connections, t, systems = [f, c, p])

lsys0, ssys = linearize(cl, [f.u], [p.x])
desired_order = [f.x, p.x]
lsys = ModelingToolkit.reorder_states(lsys0, states(ssys), desired_order)

@assert lsys.A == [-2 0; 1 -2]
@assert lsys.B == [1; 0;;]
@assert lsys.C == [0 1]
@assert lsys.D[] == 0

## Symbolic linearization
lsys_sym, _ = ModelingToolkit.linearize_symbolic(cl, [f.u], [p.x])

@assert substitute(lsys_sym.A, ModelingToolkit.defaults(cl)) == lsys.A
```
"""
function linearize(sys, lin_fun; t = 0.0, op = Dict(), allow_input_derivatives = false,
        p = DiffEqBase.NullParameters())
    x0 = merge(defaults(sys), op)
    u0, p2, _ = get_u0_p(sys, x0, p; use_union = false, tofloat = true)

    linres = lin_fun(u0, p2, t)
    f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u = linres

    nx, nu = size(f_u)
    nz = size(f_z, 2)
    ny = size(h_x, 1)

    D = h_u

    if isempty(g_z)
        A = f_x
        B = f_u
        C = h_x
        @assert iszero(g_x)
        @assert iszero(g_z)
        @assert iszero(g_u)
    else
        gz = lu(g_z; check = false)
        issuccess(gz) ||
            error("g_z not invertible, this indicates that the DAE is of index > 1.")
        gzgx = -(gz \ g_x)
        A = [f_x f_z
            gzgx*f_x gzgx*f_z]
        B = [f_u
            gzgx * f_u] # The cited paper has zeros in the bottom block, see derivation in https://github.com/SciML/ModelingToolkit.jl/pull/1691 for the correct formula

        C = [h_x h_z]
        Bs = -(gz \ g_u) # This equation differ from the cited paper, the paper is likely wrong since their equaiton leads to a dimension mismatch.
        if !iszero(Bs)
            if !allow_input_derivatives
                der_inds = findall(vec(any(!=(0), Bs, dims = 1)))
                error("Input derivatives appeared in expressions (-g_z\\g_u != 0), the following inputs appeared differentiated: $(inputs(sys)[der_inds]). Call `linearize` with keyword argument `allow_input_derivatives = true` to allow this and have the returned `B` matrix be of double width ($(2nu)), where the last $nu inputs are the derivatives of the first $nu inputs.")
            end
            B = [B [zeros(nx, nu); Bs]]
            D = [D zeros(ny, nu)]
        end
    end

    (; A, B, C, D)
end

function linearize(sys, inputs, outputs; op = Dict(), t = 0.0,
        allow_input_derivatives = false,
        zero_dummy_der = false,
        kwargs...)
    lin_fun, ssys = linearization_function(sys,
        inputs,
        outputs;
        zero_dummy_der,
        op,
        kwargs...)
    linearize(ssys, lin_fun; op, t, allow_input_derivatives), ssys
end

"""
    (; Ã, B̃, C̃, D̃) = similarity_transform(sys, T; unitary=false)

Perform a similarity transform `T : Tx̃ = x` on linear system represented by matrices in NamedTuple `sys` such that

```
Ã = T⁻¹AT
B̃ = T⁻¹ B
C̃ = CT
D̃ = D
```

If `unitary=true`, `T` is assumed unitary and the matrix adjoint is used instead of the inverse.
"""
function similarity_transform(sys::NamedTuple, T; unitary = false)
    if unitary
        A = T'sys.A * T
        B = T'sys.B
    else
        Tf = lu(T)
        A = Tf \ sys.A * T
        B = Tf \ sys.B
    end
    C = sys.C * T
    D = sys.D
    (; A, B, C, D)
end

"""
    reorder_states(sys::NamedTuple, old, new)

Permute the state representation of `sys` obtained from [`linearize`](@ref) so that the state order is changed from `old` to `new`
Example:

```
lsys, ssys = linearize(pid, [reference.u, measurement.u], [ctr_output.u])
desired_order = [int.x, der.x] # States that are present in states(ssys)
lsys = ModelingToolkit.reorder_states(lsys, states(ssys), desired_order)
```

See also [`ModelingToolkit.similarity_transform`](@ref)
"""
function reorder_states(sys::NamedTuple, old, new)
    nx = length(old)
    length(new) == nx || error("old and new must have the same length")
    perm = [findfirst(isequal(n), old) for n in new]
    issorted(perm) && return sys # shortcut return, no reordering
    P = zeros(Int, nx, nx)
    for i in 1:nx # Build permutation matrix
        P[i, perm[i]] = 1
    end
    similarity_transform(sys, P; unitary = true)
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
    s = foldr(hash, get_discrete_events(sys), init = s)
    s = hash(independent_variables(sys), s)
    return s
end

"""
$(TYPEDSIGNATURES)

extend the `basesys` with `sys`, the resulting system would inherit `sys`'s name
by default.
"""
function extend(sys::AbstractSystem, basesys::AbstractSystem; name::Symbol = nameof(sys),
        gui_metadata = get_gui_metadata(sys))
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
    cevs = union(get_continuous_events(basesys), get_continuous_events(sys))
    devs = union(get_discrete_events(basesys), get_discrete_events(sys))
    defs = merge(get_defaults(basesys), get_defaults(sys)) # prefer `sys`
    syss = union(get_systems(basesys), get_systems(sys))

    if length(ivs) == 0
        T(eqs, sts, ps, observed = obs, defaults = defs, name = name, systems = syss,
            continuous_events = cevs, discrete_events = devs, gui_metadata = gui_metadata)
    elseif length(ivs) == 1
        T(eqs, ivs[1], sts, ps, observed = obs, defaults = defs, name = name,
            systems = syss, continuous_events = cevs, discrete_events = devs,
            gui_metadata = gui_metadata)
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

"""
    missing_variable_defaults(sys::AbstractSystem, default = 0.0)

returns a `Vector{Pair}` of variables set to `default` which are missing from `get_defaults(sys)`.  The `default` argument can be a single value or vector to set the missing defaults respectively.
"""
function missing_variable_defaults(sys::AbstractSystem, default = 0.0)
    varmap = get_defaults(sys)
    varmap = Dict(Symbolics.diff2term(value(k)) => value(varmap[k]) for k in keys(varmap))
    missingvars = setdiff(states(sys), keys(varmap))
    ds = Pair[]

    n = length(missingvars)

    if default isa Vector
        @assert length(default)==n "`default` size ($(length(default))) should match the number of missing variables: $n"
    end

    for (i, missingvar) in enumerate(missingvars)
        if default isa Vector
            push!(ds, missingvar => default[i])
        else
            push!(ds, missingvar => default)
        end
    end

    return ds
end

keytype(::Type{<:Pair{T, V}}) where {T, V} = T
function Symbolics.substitute(sys::AbstractSystem, rules::Union{Vector{<:Pair}, Dict})
    if keytype(eltype(rules)) <: Symbol
        dict = todict(rules)
        systems = get_systems(sys)
        # post-walk to avoid infinite recursion
        @set! sys.systems = map(Base.Fix2(substitute, dict), systems)
        something(get(rules, nameof(sys), nothing), sys)
    elseif sys isa ODESystem
        rules = todict(map(r -> Symbolics.unwrap(r[1]) => Symbolics.unwrap(r[2]),
            collect(rules)))
        eqs = fast_substitute(equations(sys), rules)
        ODESystem(eqs, get_iv(sys); name = nameof(sys))
    else
        error("substituting symbols is not supported for $(typeof(sys))")
    end
end
