const AnyDict = Dict{Any, Any}

"""
    $(TYPEDSIGNATURES)

If called without arguments, return `Dict{Any, Any}`. Otherwise, interpret the input
as a symbolic map and turn it into a `Dict{Any, Any}`. Handles `SciMLBase.NullParameters`,
`missing` and `nothing`.
"""
anydict() = AnyDict()
anydict(::SciMLBase.NullParameters) = AnyDict()
anydict(::Nothing) = AnyDict()
anydict(::Missing) = AnyDict()
anydict(x::AnyDict) = x
anydict(x) = AnyDict(x)

"""
    $(TYPEDSIGNATURES)

Check if `x` is a symbolic with known size. Assumes `Symbolics.shape(unwrap(x))`
is a valid operation.
"""
is_sized_array_symbolic(x) = Symbolics.shape(unwrap(x)) != Symbolics.Unknown()

"""
    $(TYPEDSIGNATURES)

Check if the system is in split form (has an `IndexCache`).
"""
is_split(sys::AbstractSystem) = has_index_cache(sys) && get_index_cache(sys) !== nothing

"""
    $(TYPEDSIGNATURES)

Given a variable-value mapping, add mappings for the `toterm` of each of the keys.
"""
function add_toterms!(varmap::AbstractDict; toterm = default_toterm)
    for k in collect(keys(varmap))
        varmap[toterm(k)] = varmap[k]
    end
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Out-of-place version of [`add_toterms!`](@ref).
"""
function add_toterms(varmap::AbstractDict; toterm = default_toterm)
    cp = copy(varmap)
    add_toterms!(cp; toterm)
    return cp
end

"""
    $(TYPEDSIGNATURES)

Turn any `Symbol` keys in `varmap` to the appropriate symbolic variables in `sys`. Any
symbols that cannot be converted are ignored.
"""
function symbols_to_symbolics!(sys::AbstractSystem, varmap::AbstractDict)
    if is_split(sys)
        ic = get_index_cache(sys)
        for k in collect(keys(varmap))
            k isa Symbol || continue
            newk = get(ic.symbol_to_variable, k, nothing)
            newk === nothing && continue
            varmap[newk] = varmap[k]
            delete!(varmap, k)
        end
    else
        syms = all_symbols(sys)
        for k in collect(keys(varmap))
            k isa Symbol || continue
            idx = findfirst(syms) do sym
                hasname(sym) || return false
                name = getname(sym)
                return name == k
            end
            idx === nothing && continue
            newk = syms[idx]
            if iscall(newk) && operation(newk) === getindex
                newk = arguments(newk)[1]
            end
            varmap[newk] = varmap[k]
            delete!(varmap, k)
        end
    end
end

"""
    $(TYPEDSIGNATURES)

Utility function to get the value `val` corresponding to key `var` in `varmap`, and
return `getindex(val, idx)` if it exists or `nothing` otherwise.
"""
function get_and_getindex(varmap, var, idx)
    val = get(varmap, var, nothing)
    val === nothing && return nothing
    return val[idx]
end

"""
    $(TYPEDSIGNATURES)

Ensure `varmap` contains entries for all variables in `vars` by using values from
`fallbacks` if they don't already exist in `varmap`. Return the set of all variables in
`vars` not present in `varmap` or `fallbacks`. If an array variable in `vars` does not
exist in `varmap` or `fallbacks`, each of its scalarized elements will be searched for.
In case none of the scalarized elements exist, the array variable will be reported as
missing. In case some of the scalarized elements exist, the missing elements will be
reported as missing. If `fallbacks` contains both the scalarized and non-scalarized forms,
the latter will take priority.

Variables as they are specified in `vars` will take priority over their `toterm` forms.
"""
function add_fallbacks!(
        varmap::AnyDict, vars::Vector, fallbacks::Dict; toterm = default_toterm)
    missingvars = Set()
    arrvars = Set()
    for var in vars
        haskey(varmap, var) && continue
        ttvar = toterm(var)
        haskey(varmap, ttvar) && continue

        # array symbolics with a defined size may be present in the scalarized form
        if Symbolics.isarraysymbolic(var) && is_sized_array_symbolic(var)
            val = map(eachindex(var)) do idx
                # @something is lazy and saves from writing a massive if-elseif-else
                @something(get(varmap, var[idx], nothing),
                    get(varmap, ttvar[idx], nothing), get_and_getindex(fallbacks, var, idx),
                    get_and_getindex(fallbacks, ttvar, idx), get(
                        fallbacks, var[idx], nothing),
                    get(fallbacks, ttvar[idx], nothing), Some(nothing))
            end
            # only push the missing entries
            mask = map(x -> x === nothing, val)
            if all(mask)
                push!(missingvars, var)
            elseif any(mask)
                for i in eachindex(var)
                    if mask[i]
                        push!(missingvars, var)
                    else
                        varmap[var[i]] = val[i]
                    end
                end
            else
                varmap[var] = val
            end
        else
            if iscall(var) && operation(var) == getindex
                args = arguments(var)
                arrvar = args[1]
                ttarrvar = toterm(arrvar)
                idxs = args[2:end]
                val = @something get(varmap, arrvar, nothing) get(varmap, ttarrvar, nothing) get(
                    fallbacks, arrvar, nothing) get(fallbacks, ttarrvar, nothing) Some(nothing)
                if val !== nothing
                    val = val[idxs...]
                    is_sized_array_symbolic(arrvar) && push!(arrvars, arrvar)
                end
            else
                val = nothing
            end
            val = @something val get(fallbacks, var, nothing) get(fallbacks, ttvar, nothing) Some(nothing)
            if val === nothing
                push!(missingvars, var)
            else
                varmap[var] = val
            end
        end
    end

    for arrvar in arrvars
        varmap[arrvar] = collect(arrvar)
    end

    return missingvars
end

"""
    $(TYPEDSIGNATURES)

Return the list of variables in `varlist` not present in `varmap`. Uses the same criteria
for missing array variables and `toterm` forms as [`add_fallbacks!`](@ref).
"""
function missingvars(
        varmap::AbstractDict, varlist::Vector; toterm = default_toterm)
    missingvars = Set()
    for var in varlist
        haskey(varmap, var) && continue
        ttsym = toterm(var)
        haskey(varmap, ttsym) && continue

        if Symbolics.isarraysymbolic(var) && is_sized_array_symbolic(var)
            mask = map(eachindex(var)) do idx
                !haskey(varmap, var[idx]) && !haskey(varmap, ttsym[idx])
            end
            if all(mask)
                push!(missingvars, var)
            else
                for i in eachindex(var)
                    mask[i] && push!(missingvars, var[i])
                end
            end
        else
            push!(missingvars, var)
        end
    end
    return missingvars
end

"""
    $(TYPEDSIGNATURES)

Attempt to interpret `vals` as a symbolic map of variables in `varlist` to values. Return
the result as a `Dict{Any, Any}`. In case `vals` is already an iterable of pairs, convert
it to a `Dict{Any, Any}` and return. If `vals` is an array (whose `eltype` is not `Pair`)
with the same length as `varlist`, assume the `i`th element of `varlist` is mapped to the
`i`th element of `vals`. Automatically `unwrap`s all keys and values in the mapping. Also
handles `SciMLBase.NullParameters` and `nothing`, both of which are interpreted as empty
maps.
"""
function to_varmap(vals, varlist::Vector)
    if vals isa AbstractArray && !(eltype(vals) <: Pair) && !isempty(vals)
        check_eqs_u0(varlist, varlist, vals)
        vals = vec(varlist) .=> vec(vals)
    end
    return recursive_unwrap(anydict(vals))
end

"""
    $(TYPEDSIGNATURES)

Recursively call `Symbolics.unwrap` on `x`. Useful when `x` is an array of (potentially)
symbolic values, all of which need to be unwrapped. Specializes when `x isa AbstractDict`
to unwrap keys and values, returning an `AnyDict`.
"""
function recursive_unwrap(x::AbstractArray)
    symbolic_type(x) == ArraySymbolic() ? unwrap(x) : recursive_unwrap.(x)
end

recursive_unwrap(x) = unwrap(x)

function recursive_unwrap(x::AbstractDict)
    return anydict(unwrap(k) => recursive_unwrap(v) for (k, v) in x)
end

"""
    $(TYPEDSIGNATURES)

Add equations `eqs` to `varmap`. Assumes each element in `eqs` maps a single symbolic
variable to an expression representing its value. In case `varmap` already contains an
entry for `eq.lhs`, insert the reverse mapping if `eq.rhs` is not a number.
"""
function add_observed_equations!(varmap::AbstractDict, eqs)
    for eq in eqs
        if var_in_varlist(eq.lhs, keys(varmap), nothing)
            eq.rhs isa Number && continue
            var_in_varlist(eq.rhs, keys(varmap), nothing) && continue
            !iscall(eq.rhs) || issym(operation(eq.rhs)) || continue
            varmap[eq.rhs] = eq.lhs
        else
            varmap[eq.lhs] = eq.rhs
        end
    end
end

"""
    $(TYPEDSIGNATURES)

Add all equations in `observed(sys)` to `varmap` using [`add_observed_equations!`](@ref).
"""
function add_observed!(sys::AbstractSystem, varmap::AbstractDict)
    add_observed_equations!(varmap, observed(sys))
end

"""
    $(TYPEDSIGNATURES)

Add all equations in `parameter_dependencies(sys)` to `varmap` using
[`add_observed_equations!`](@ref).
"""
function add_parameter_dependencies!(sys::AbstractSystem, varmap::AbstractDict)
    has_parameter_dependencies(sys) || return nothing
    add_observed_equations!(varmap, parameter_dependencies(sys))
end

struct UnexpectedSymbolicValueInVarmap <: Exception
    sym::Any
    val::Any
end

function Base.showerror(io::IO, err::UnexpectedSymbolicValueInVarmap)
    println(io,
        """
        Found symbolic value $(err.val) for variable $(err.sym). You may be missing an \
        initial condition or have cyclic initial conditions. If this is intended, pass \
        `symbolic_u0 = true`. In case the initial conditions are not cyclic but \
        require more substitutions to resolve, increase `substitution_limit`. To report \
        cycles in initial conditions of unknowns/parameters, pass \
        `warn_cyclic_dependency = true`. If the cycles are still not reported, you \
        may need to pass a larger value for `circular_dependency_max_cycle_length` \
        or `circular_dependency_max_cycles`.
        """)
end

struct MissingGuessError <: Exception
    syms::Vector{Any}
    vals::Vector{Any}
end

function Base.showerror(io::IO, err::MissingGuessError)
    println(io,
        """
        Cyclic guesses detected in the system. Symbolic values were found for the following variables/parameters in the map: \
        """)
    for (sym, val) in zip(err.syms, err.vals)
        println(io, "$sym  => $val")
    end
    println(io,
        """
        In order to resolve this, please provide additional numeric guesses so that the chain can be resolved to assign numeric values to each variable.            """)
end

"""
    $(TYPEDSIGNATURES)

Return an array of values where the `i`th element corresponds to the value of `vars[i]`
in `varmap`. Does not perform symbolic substitution in the values of `varmap`.

Keyword arguments:
- `tofloat`: Convert values to floating point numbers using `float`.
- `container_type`: The type of container to use for the values.
- `toterm`: The `toterm` method to use for converting symbolics.
- `promotetoconcrete`: whether the promote to a concrete buffer (respecting
  `tofloat`). Defaults to `container_type <: AbstractArray`.
- `check`: Error if any variables in `vars` do not have a mapping in `varmap`. Uses
  [`missingvars`](@ref) to perform the check.
- `allow_symbolic` allows the returned array to contain symbolic values. If this is `true`,
  `promotetoconcrete` is set to `false`.
- `is_initializeprob, guesses`: Used to determine whether the system is missing guesses.
"""
function better_varmap_to_vars(varmap::AbstractDict, vars::Vector;
        tofloat = true, container_type = Array, floatT = Nothing,
        toterm = default_toterm, promotetoconcrete = nothing, check = true,
        allow_symbolic = false, is_initializeprob = false)
    isempty(vars) && return nothing

    if check
        missing_vars = missingvars(varmap, vars; toterm)
        isempty(missing_vars) || throw(MissingVariablesError(missing_vars))
    end
    vals = map(x -> varmap[x], vars)
    if !allow_symbolic
        missingsyms = Any[]
        missingvals = Any[]
        for (sym, val) in zip(vars, vals)
            symbolic_type(val) == NotSymbolic() && continue
            push!(missingsyms, sym)
            push!(missingvals, val)
        end

        if !isempty(missingsyms)
            is_initializeprob ? throw(MissingGuessError(missingsyms, missingvals)) :
            throw(UnexpectedSymbolicValueInVarmap(missingsyms[1], missingvals[1]))
        end
        if tofloat && !(floatT == Nothing)
            vals = floatT.(vals)
        end
    end

    if container_type <: Union{AbstractDict, Tuple, Nothing, SciMLBase.NullParameters}
        container_type = Array
    end

    promotetoconcrete === nothing && (promotetoconcrete = container_type <: AbstractArray)
    if promotetoconcrete && !allow_symbolic
        vals = promote_to_concrete(vals; tofloat = tofloat, use_union = false)
    end

    if container_type <: Tuple
        return (vals...,)
    else
        return SymbolicUtils.Code.create_array(container_type, eltype(vals), Val{1}(),
            Val(length(vals)), vals...)
    end
end

"""
    $(TYPEDSIGNATURES)

Check if any of the substitution rules in `varmap` lead to cycles involving
variables in `vars`. Return a vector of vectors containing all the variables
in each cycle.

Keyword arguments:
- `max_cycle_length`: The maximum length (number of variables) of detected cycles.
- `max_cycles`: The maximum number of cycles to report.
"""
function check_substitution_cycles(
        varmap::AbstractDict, vars; max_cycle_length = length(varmap), max_cycles = 10)
    # ordered set so that `vars` are the first `k` in the list
    allvars = OrderedSet{Any}(vars)
    union!(allvars, keys(varmap))
    allvars = collect(allvars)
    var_to_idx = Dict(allvars .=> eachindex(allvars))
    graph = SimpleDiGraph(length(allvars))

    buffer = Set()
    for (k, v) in varmap
        kidx = var_to_idx[k]
        if symbolic_type(v) != NotSymbolic()
            vars!(buffer, v)
            for var in buffer
                haskey(var_to_idx, var) || continue
                add_edge!(graph, kidx, var_to_idx[var])
            end
        elseif v isa AbstractArray
            for val in v
                vars!(buffer, val)
            end
            for var in buffer
                haskey(var_to_idx, var) || continue
                add_edge!(graph, kidx, var_to_idx[var])
            end
        end
        empty!(buffer)
    end

    # detect at most 100 cycles involving at most `length(varmap)` vertices
    cycles = Graphs.simplecycles_limited_length(graph, max_cycle_length, max_cycles)
    # only count those which contain variables in `vars`
    filter!(Base.Fix1(any, <=(length(vars))), cycles)

    map(cycles) do cycle
        map(Base.Fix1(getindex, allvars), cycle)
    end
end

"""
    $(TYPEDSIGNATURES)

Performs symbolic substitution on the values in `varmap` for the keys in `vars`, using
`varmap` itself as the set of substitution rules. If an entry in `vars` is not a key
in `varmap`, it is ignored.
"""
function evaluate_varmap!(varmap::AbstractDict, vars; limit = 100)
    for k in vars
        haskey(varmap, k) || continue
        varmap[k] = fixpoint_sub(varmap[k], varmap; maxiters = limit)
    end
end

"""
    $(TYPEDSIGNATURES)

Remove keys in `varmap` whose values are `nothing`.

If `missing_values` is not `nothing`, it is assumed to be a collection and all removed
keys will be added to it.
"""
function filter_missing_values!(varmap::AbstractDict; missing_values = nothing)
    filter!(varmap) do kvp
        if kvp[2] !== nothing
            return true
        end
        if missing_values !== nothing
            push!(missing_values, kvp[1])
        end
        return false
    end
end

"""
    $(TYPEDSIGNATURES)

For each `k => v` in `varmap` where `k` is an array (or array symbolic) add
`k[i] => v[i]` for all `i  in eachindex(k)`. Return the modified `varmap`.
"""
function scalarize_varmap!(varmap::AbstractDict)
    for k in collect(keys(varmap))
        symbolic_type(k) == ArraySymbolic() || continue
        for i in eachindex(k)
            varmap[k[i]] = varmap[k][i]
        end
    end
    return varmap
end

struct GetUpdatedMTKParameters{G, S}
    # `getu` functor which gets parameters that are unknowns during initialization
    getpunknowns::G
    # `setu` functor which returns a modified MTKParameters using those parameters
    setpunknowns::S
end

function (f::GetUpdatedMTKParameters)(prob, initializesol)
    p = parameter_values(prob)
    p === nothing && return nothing
    mtkp = copy(p)
    f.setpunknowns(mtkp, f.getpunknowns(initializesol))
    mtkp
end

struct UpdateInitializeprob{G, S}
    # `getu` functor which gets all values from prob
    getvals::G
    # `setu` functor which updates initializeprob with values
    setvals::S
end

function (f::UpdateInitializeprob)(initializeprob, prob)
    f.setvals(initializeprob, f.getvals(prob))
end

function get_temporary_value(p, floatT = Float64)
    stype = symtype(unwrap(p))
    return if stype == Real
        zero(floatT)
    elseif stype <: AbstractArray{Real}
        zeros(floatT, size(p))
    elseif stype <: Real
        zero(stype)
    elseif stype <: AbstractArray
        zeros(eltype(stype), size(p))
    else
        error("Nonnumeric parameter $p with symtype $stype cannot be solved for during initialization")
    end
end

"""
    $(TYPEDEF)

A simple utility meant to be used as the `constructor` passed to `process_SciMLProblem` in
case constructing a SciMLFunction is not required. The arguments passed to it are available
in the `args` field, and the keyword arguments in the `kwargs` field.
"""
struct EmptySciMLFunction{A, K}
    args::A
    kwargs::K
end

function EmptySciMLFunction(args...; kwargs...)
    return EmptySciMLFunction{typeof(args), typeof(kwargs)}(args, kwargs)
end

"""
    $(TYPEDSIGNATURES)

Construct the operating point of the system from the user-provided `u0map` and `pmap`, system
defaults `defs`, constant equations `cmap` (from `get_cmap(sys)`), unknowns `dvs` and
parameters `ps`. Return the operating point as a dictionary, the list of unknowns for which
no values can be determined, and the list of parameters for which no values can be determined.

Also updates `u0map` and `pmap` in-place to contain all the initial conditions in `op`, split
by unknowns and parameters respectively.
"""
function build_operating_point!(sys::AbstractSystem,
        u0map::AbstractDict, pmap::AbstractDict, defs::AbstractDict, cmap, dvs, ps)
    op = add_toterms(u0map)
    missing_unknowns = add_fallbacks!(op, dvs, defs)
    for (k, v) in defs
        haskey(op, k) && continue
        op[k] = v
    end
    filter_missing_values!(op; missing_values = missing_unknowns)

    merge!(op, pmap)
    missing_pars = add_fallbacks!(op, ps, defs)
    filter_missing_values!(op; missing_values = missing_pars)
    for eq in cmap
        op[eq.lhs] = eq.rhs
    end

    filter!(kvp -> kvp[2] === nothing, u0map)
    filter!(kvp -> kvp[2] === nothing, pmap)
    neithermap = anydict()

    for (k, v) in op
        k = unwrap(k)
        if is_parameter(sys, k)
            pmap[k] = v
        elseif has_parameter_dependency_with_lhs(sys, k) && is_variable_floatingpoint(k) &&
               v !== nothing && !isequal(v, Initial(k))
            op[Initial(k)] = v
            pmap[Initial(k)] = v
            op[k] = Initial(k)
            pmap[k] = Initial(k)
        elseif is_variable(sys, k) || has_observed_with_lhs(sys, k) ||
               iscall(k) &&
               operation(k) isa Differential && is_variable(sys, arguments(k)[1])
            if symbolic_type(v) == NotSymbolic() && !is_array_of_symbolics(v) &&
               v !== nothing
                op[Initial(k)] = v
                pmap[Initial(k)] = v
                op[k] = Initial(k)
                v = Initial(k)
            end
            u0map[k] = v
        else
            neithermap[k] = v
        end
    end

    for k in keys(u0map)
        v = fixpoint_sub(u0map[k], neithermap; operator = Symbolics.Operator)
        isequal(k, v) && continue
        u0map[k] = v
    end
    for k in keys(pmap)
        v = fixpoint_sub(pmap[k], neithermap; operator = Symbolics.Operator)
        isequal(k, v) && continue
        pmap[k] = v
    end

    return op, missing_unknowns, missing_pars
end

"""
    $(TYPEDEF)

A callable struct used to reconstruct the `u0` and `p` of the initialization problem
with promoted types.

# Fields

$(TYPEDFIELDS)
"""
struct ReconstructInitializeprob{G}
    """
    A function which when given the original problem and initialization problem, returns
    the parameter object of the initialization problem with values copied from the
    original.
    """
    getter::G
end

"""
    $(TYPEDSIGNATURES)

Given an index provider `indp` and a vector of symbols `syms` return a type-stable getter
function by splitting `syms` into contiguous buffers where the getter of each buffer
is type-stable and constructing a function that calls and concatenates the results.
"""
function concrete_getu(indp, syms::AbstractVector)
    # a list of contiguous buffer
    split_syms = [Any[syms[1]]]
    # the type of the getter of the last buffer
    current = typeof(getu(indp, syms[1]))
    for sym in syms[2:end]
        getter = getu(indp, sym)
        if typeof(getter) != current
            # if types don't match, build a new buffer
            push!(split_syms, [])
            current = typeof(getter)
        end
        push!(split_syms[end], sym)
    end
    split_syms = Tuple(split_syms)
    # the getter is now type-stable, and we can vcat it to get the full buffer
    return Base.Fix1(reduce, vcat) ∘ getu(indp, split_syms)
end

"""
    $(TYPEDSIGNATURES)

Construct a `ReconstructInitializeprob` which reconstructs the `u0` and `p` of `dstsys`
with values from `srcsys`.
"""
function ReconstructInitializeprob(
        srcsys::AbstractSystem, dstsys::AbstractSystem)
    @assert is_initializesystem(dstsys)
    if is_split(dstsys)
        # if we call `getu` on this (and it were able to handle empty tuples) we get the
        # fields of `MTKParameters` except caches.
        syms = reorder_parameters(dstsys, parameters(dstsys); flatten = false)
        # `dstsys` is an initialization system, do basically everything is a tunable
        # and tunables are a mix of different types in `srcsys`. No initials. Constants
        # are going to be constants in `srcsys`, as are `nonnumeric`.

        # `syms[1]` is always the tunables because `srcsys` will have initials.
        tunable_syms = syms[1]
        tunable_getter = concrete_getu(srcsys, tunable_syms)
        rest_getters = map(Base.tail(Base.tail(syms))) do buf
            if buf == ()
                return Returns(())
            else
                return getu(srcsys, buf)
            end
        end
        getters = (tunable_getter, Returns(SizedVector{0, Float64}()), rest_getters...)
        getter = let getters = getters
            function _getter(valp, initprob)
                oldcache = parameter_values(initprob).caches
                MTKParameters(getters[1](valp), getters[2](valp), getters[3](valp),
                    getters[4](valp), getters[5](valp), oldcache isa Tuple{} ? () :
                                                        copy.(oldcache))
            end
        end
    else
        syms = parameters(dstsys)
        getter = let inner = concrete_getu(srcsys, syms)
            function _getter2(valp, initprob)
                inner(valp)
            end
        end
    end
    return ReconstructInitializeprob(getter)
end

"""
    $(TYPEDSIGNATURES)

Copy values from `srcvalp` to `dstvalp`. Returns the new `u0` and `p`.
"""
function (rip::ReconstructInitializeprob)(srcvalp, dstvalp)
    # copy parameters
    newp = rip.getter(srcvalp, dstvalp)
    # no `u0`, so no type-promotion
    if state_values(dstvalp) === nothing
        return nothing, newp
    end
    # the `eltype` of the `u0` of the source
    srcu0 = state_values(srcvalp)
    T = srcu0 === nothing ? Union{} : eltype(srcu0)
    # promote with the tunable eltype
    if parameter_values(dstvalp) isa MTKParameters
        if !isempty(newp.tunable)
            T = promote_type(eltype(newp.tunable), T)
        end
    elseif !isempty(newp)
        T = promote_type(eltype(newp), T)
    end
    # and the eltype of the destination u0
    if T == eltype(state_values(dstvalp))
        u0 = state_values(dstvalp)
    elseif T != Union{}
        u0 = T.(state_values(dstvalp))
    end
    # apply the promotion to tunables portion
    buf, repack, alias = SciMLStructures.canonicalize(SciMLStructures.Tunable(), newp)
    if eltype(buf) != T
        # only do a copy if the eltype doesn't match
        newbuf = similar(buf, T)
        copyto!(newbuf, buf)
        newp = repack(newbuf)
    end
    # and initials portion
    buf, repack, alias = SciMLStructures.canonicalize(SciMLStructures.Initials(), newp)
    if eltype(buf) != T
        newbuf = similar(buf, T)
        copyto!(newbuf, buf)
        newp = repack(newbuf)
    end
    return u0, newp
end

"""
    $(TYPEDEF)

Metadata attached to `OverrideInitData` used in `remake` hooks for handling initialization
properly.

# Fields

$(TYPEDFIELDS)
"""
struct InitializationMetadata{R <: ReconstructInitializeprob, GUU, SIU}
    """
    The `u0map` used to construct the initialization.
    """
    u0map::Dict{Any, Any}
    """
    The `pmap` used to construct the initialization.
    """
    pmap::Dict{Any, Any}
    """
    The `guesses` used to construct the initialization.
    """
    guesses::Dict{Any, Any}
    """
    The `initialization_eqs` in addition to those of the system that were used to construct
    the initialization.
    """
    additional_initialization_eqs::Vector{Equation}
    """
    Whether to use `SCCNonlinearProblem` if possible.
    """
    use_scc::Bool
    """
    `ReconstructInitializeprob` for this initialization problem.
    """
    oop_reconstruct_u0_p::R
    """
    A function which takes `(prob, initializeprob)` and return the `u0` to use for the problem.
    """
    get_updated_u0::GUU
    """
    A function which takes the `u0` of the problem and sets
    `Initial.(unknowns(sys))`.
    """
    set_initial_unknowns!::SIU
end

"""
    $(TYPEDEF)

A callable struct to use as the `get_updated_u0` field of `InitializationMetadata`.
Returns the value to use for the `u0` of the problem. 

# Fields

$(TYPEDFIELDS)
"""
struct GetUpdatedU0{GG, GIU}
    """
    Mask with length `length(unknowns(sys))` denoting indices of variables which should
    take the guess value from `initializeprob`.
    """
    guessvars::BitVector
    """
    Function which returns the values of variables in `initializeprob` for which
    `guessvars` is `true`, in the order they occur in `unknowns(sys)`.
    """
    get_guessvars::GG
    """
    Function which returns `Initial.(unknowns(sys))` as a `Vector`.
    """
    get_initial_unknowns::GIU
end

function GetUpdatedU0(sys::AbstractSystem, initsys::AbstractSystem, op::AbstractDict)
    dvs = unknowns(sys)
    eqs = equations(sys)
    guessvars = trues(length(dvs))
    for (i, var) in enumerate(dvs)
        guessvars[i] = !isequal(get(op, var, nothing), Initial(var))
    end
    get_guessvars = getu(initsys, dvs[guessvars])
    get_initial_unknowns = getu(sys, Initial.(dvs))
    return GetUpdatedU0(guessvars, get_guessvars, get_initial_unknowns)
end

function (guu::GetUpdatedU0)(prob, initprob)
    buffer = guu.get_initial_unknowns(prob)
    algebuf = view(buffer, guu.guessvars)
    copyto!(algebuf, guu.get_guessvars(initprob))
    return buffer
end

"""
    $(TYPEDSIGNATURES)

Build and return the initialization problem and associated data as a `NamedTuple` to be passed
to the `SciMLFunction` constructor. Requires the system `sys`, operating point `op`,
user-provided `u0map` and `pmap`, initial time `t`, system defaults `defs`, user-provided
`guesses`, and list of unknowns which don't have a value in `op`. The keyword `implicit_dae`
denotes whether the `SciMLProblem` being constructed is in implicit DAE form (`DAEProblem`).
All other keyword arguments are forwarded to `InitializationProblem`.
"""
function maybe_build_initialization_problem(
        sys::AbstractSystem, op::AbstractDict, u0map, pmap, t, defs,
        guesses, missing_unknowns; implicit_dae = false, u0_constructor = identity,
        floatT = Float64, initialization_eqs = [], use_scc = true, kwargs...)
    guesses = merge(ModelingToolkit.guesses(sys), todict(guesses))

    if t === nothing && is_time_dependent(sys)
        t = zero(floatT)
    end

    initializeprob = ModelingToolkit.InitializationProblem{true, SciMLBase.FullSpecialize}(
        sys, t, u0map, pmap; guesses, initialization_eqs, use_scc, kwargs...)
    if state_values(initializeprob) !== nothing
        initializeprob = remake(initializeprob; u0 = floatT.(state_values(initializeprob)))
    end
    initp = parameter_values(initializeprob)
    if is_split(sys)
        buffer, repack, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), initp)
        initp = repack(floatT.(buffer))
        buffer, repack, _ = SciMLStructures.canonicalize(SciMLStructures.Initials(), initp)
        initp = repack(floatT.(buffer))
    elseif initp isa AbstractArray
        initp′ = similar(initp, floatT)
        copyto!(initp′, initp)
        initp = initp′
    end
    initializeprob = remake(initializeprob; p = initp)

    get_initial_unknowns = if is_time_dependent(sys)
        GetUpdatedU0(sys, initializeprob.f.sys, op)
    else
        nothing
    end
    meta = InitializationMetadata(
        u0map, pmap, guesses, Vector{Equation}(initialization_eqs),
        use_scc, ReconstructInitializeprob(sys, initializeprob.f.sys),
        get_initial_unknowns, setp(sys, Initial.(unknowns(sys))))

    if is_time_dependent(sys)
        all_init_syms = Set(all_symbols(initializeprob))
        solved_unknowns = filter(var -> var in all_init_syms, unknowns(sys))
        initializeprobmap = u0_constructor ∘ getu(initializeprob, solved_unknowns)
    else
        initializeprobmap = nothing
    end

    punknowns = [p
                 for p in all_variable_symbols(initializeprob)
                 if is_parameter(sys, p)]
    if isempty(punknowns)
        initializeprobpmap = nothing
    else
        getpunknowns = getu(initializeprob, punknowns)
        setpunknowns = setp(sys, punknowns)
        initializeprobpmap = GetUpdatedMTKParameters(getpunknowns, setpunknowns)
    end

    reqd_syms = parameter_symbols(initializeprob)
    # we still want the `initialization_data` because it helps with `remake`
    if initializeprobmap === nothing && initializeprobpmap === nothing
        update_initializeprob! = nothing
    else
        update_initializeprob! = UpdateInitializeprob(
            getu(sys, reqd_syms), setu(initializeprob, reqd_syms))
    end

    for p in punknowns
        is_parameter_solvable(p, pmap, defs, guesses) || continue
        get(op, p, missing) === missing || continue
        p = unwrap(p)
        stype = symtype(p)
        op[p] = get_temporary_value(p, floatT)
        if iscall(p) && operation(p) === getindex
            arrp = arguments(p)[1]
            op[arrp] = collect(arrp)
        end
    end

    if is_time_dependent(sys)
        for v in missing_unknowns
            op[v] = get_temporary_value(v, floatT)
        end
        empty!(missing_unknowns)
    end

    return (;
        initialization_data = SciMLBase.OverrideInitData(
            initializeprob, update_initializeprob!, initializeprobmap,
            initializeprobpmap; metadata = meta))
end

"""
    $(TYPEDSIGNATURES)

Calculate the floating point type to use from the given `varmap` by looking at variables
with a constant value.
"""
function float_type_from_varmap(varmap, floatT = Bool)
    for (k, v) in varmap
        symbolic_type(v) == NotSymbolic() || continue
        is_array_of_symbolics(v) && continue

        if v isa AbstractArray
            floatT = promote_type(floatT, eltype(v))
        elseif v isa Real
            floatT = promote_type(floatT, typeof(v))
        end
    end
    return float(floatT)
end

"""
    $(TYPEDSIGNATURES)

Return the SciMLFunction created via calling `constructor`, the initial conditions `u0`
and parameter object `p` given the system `sys`, and user-provided initial values `u0map`
and `pmap`. `u0map` and `pmap` are converted into variable maps via [`to_varmap`](@ref).

The order of unknowns is determined by `unknowns(sys)`. If the system is split
[`is_split`](@ref) create an [`MTKParameters`](@ref) object. Otherwise, a parameter vector.
Initial values provided in terms of other variables will be symbolically evaluated using
[`evaluate_varmap!`](@ref). The type of `u0map` and `pmap` will be used to determine the
type of the containers (if parameters are not in an `MTKParameters` object). `Dict`s will be
turned into `Array`s.

If `sys isa ODESystem`, this will also build the initialization problem and related objects
and pass them to the SciMLFunction as keyword arguments.

Keyword arguments:
- `build_initializeprob`: If `false`, avoids building the initialization problem for an
  `ODESystem`.
- `t`: The initial time of the `ODEProblem`. If this is not provided, the initialization
  problem cannot be built.
- `implicit_dae`: Also build a mapping of derivatives of states to values for implicit DAEs,
  using `du0map`. Changes the return value of this function to `(f, du0, u0, p)` instead of
  `(f, u0, p)`.
- `guesses`: The guesses for variables in the system, used as initial values for the
  initialization problem.
- `warn_initialize_determined`: Warn if the initialization system is under/over-determined.
- `initialization_eqs`: Extra equations to use in the initialization problem.
- `eval_expression`: Whether to compile any functions via `eval` or `RuntimeGeneratedFunctions`.
- `eval_module`: If `eval_expression == true`, the module to `eval` into. Otherwise, the module
  in which to generate the `RuntimeGeneratedFunction`.
- `fully_determined`: Override whether the initialization system is fully determined.
- `check_initialization_units`: Enable or disable unit checks when constructing the
  initialization problem.
- `tofloat`, `is_initializeprob`: Passed to [`better_varmap_to_vars`](@ref) for building `u0` (and possibly `p`).
- `u0_constructor`: A function to apply to the `u0` value returned from `better_varmap_to_vars`
  to construct the final `u0` value.
- `du0map`: A map of derivatives to values. See `implicit_dae`.
- `check_length`: Whether to check the number of equations along with number of unknowns and
  length of `u0` vector for consistency. If `false`, do not check with equations. This is
  forwarded to `check_eqs_u0`
- `symbolic_u0` allows the returned `u0` to be an array of symbolics.
- `warn_cyclic_dependency`: Whether to emit a warning listing out cycles in initial
  conditions provided for unknowns and parameters.
- `circular_dependency_max_cycle_length`: Maximum length of cycle to check for.
  Only applicable if `warn_cyclic_dependency == true`.
- `circular_dependency_max_cycles`: Maximum number of cycles to check for.
  Only applicable if `warn_cyclic_dependency == true`.
- `substitution_limit`: The number times to substitute initial conditions into each
  other to attempt to arrive at a numeric value.
- `use_scc`: Whether to use `SCCNonlinearProblem` for initialization if the system is fully
  determined.
- `force_initialization_time_independent`: Whether to force the initialization to not use
  the independent variable of `sys`.
- `algebraic_only`: Whether to build the initialization problem using only algebraic equations.
- `allow_incomplete`: Whether to allow incomplete initialization problems.

All other keyword arguments are passed as-is to `constructor`.
"""
function process_SciMLProblem(
        constructor, sys::AbstractSystem, u0map, pmap; build_initializeprob = true,
        implicit_dae = false, t = nothing, guesses = AnyDict(),
        warn_initialize_determined = true, initialization_eqs = [],
        eval_expression = false, eval_module = @__MODULE__, fully_determined = nothing,
        check_initialization_units = false, tofloat = true,
        u0_constructor = identity, du0map = nothing, check_length = true,
        symbolic_u0 = false, warn_cyclic_dependency = false,
        circular_dependency_max_cycle_length = length(all_symbols(sys)),
        circular_dependency_max_cycles = 10,
        substitution_limit = 100, use_scc = true,
        force_initialization_time_independent = false, algebraic_only = false,
        allow_incomplete = false, is_initializeprob = false, kwargs...)
    dvs = unknowns(sys)
    ps = parameters(sys; initial_parameters = true)
    iv = has_iv(sys) ? get_iv(sys) : nothing
    eqs = equations(sys)

    check_array_equations_unknowns(eqs, dvs)

    u0Type = typeof(u0map)
    pType = typeof(pmap)

    u0map = to_varmap(u0map, dvs)
    symbols_to_symbolics!(sys, u0map)
    pmap = to_varmap(pmap, parameters(sys))
    symbols_to_symbolics!(sys, pmap)

    check_inputmap_keys(sys, u0map, pmap)

    defs = add_toterms(recursive_unwrap(defaults(sys)))
    cmap, cs = get_cmap(sys)
    kwargs = NamedTuple(kwargs)

    if eltype(eqs) <: Equation
        obs, eqs = unhack_observed(observed(sys), eqs)
    else
        obs, _ = unhack_observed(observed(sys), Equation[x for x in eqs if x isa Equation])
    end

    op, missing_unknowns, missing_pars = build_operating_point!(sys,
        u0map, pmap, defs, cmap, dvs, ps)

    floatT = Bool
    if u0Type <: AbstractArray && eltype(u0Type) <: Real
        floatT = float(eltype(u0Type))
    else
        floatT = float_type_from_varmap(op, floatT)
    end

    if !is_time_dependent(sys) || is_initializesystem(sys)
        add_observed_equations!(u0map, obs)
    end
    if u0_constructor === identity && u0Type <: StaticArray
        u0_constructor = vals -> SymbolicUtils.Code.create_array(
            u0Type, floatT, Val(1), Val(length(vals)), vals...)
    end
    if build_initializeprob
        kws = maybe_build_initialization_problem(
            sys, op, u0map, pmap, t, defs, guesses, missing_unknowns;
            implicit_dae, warn_initialize_determined, initialization_eqs,
            eval_expression, eval_module, fully_determined,
            warn_cyclic_dependency, check_units = check_initialization_units,
            circular_dependency_max_cycle_length, circular_dependency_max_cycles, use_scc,
            force_time_independent = force_initialization_time_independent, algebraic_only, allow_incomplete,
            u0_constructor, floatT)

        kwargs = merge(kwargs, kws)
    end

    if t !== nothing && !(constructor <: Union{DDEFunction, SDDEFunction})
        op[iv] = t
    end

    add_observed_equations!(op, obs)
    add_parameter_dependencies!(sys, op)

    if warn_cyclic_dependency
        cycles = check_substitution_cycles(
            op, dvs; max_cycle_length = circular_dependency_max_cycle_length,
            max_cycles = circular_dependency_max_cycles)
        if !isempty(cycles)
            buffer = IOBuffer()
            for cycle in cycles
                println(buffer, cycle)
            end
            msg = String(take!(buffer))
            @warn "Cycles in unknowns:\n$msg"
        end
    end
    evaluate_varmap!(op, dvs; limit = substitution_limit)

    u0 = better_varmap_to_vars(
        op, dvs; tofloat, floatT,
        container_type = u0Type, allow_symbolic = symbolic_u0, is_initializeprob)

    if u0 !== nothing
        u0 = u0_constructor(u0)
    end

    check_eqs_u0(eqs, dvs, u0; check_length, kwargs...)

    if warn_cyclic_dependency
        cycles = check_substitution_cycles(
            op, ps; max_cycle_length = circular_dependency_max_cycle_length,
            max_cycles = circular_dependency_max_cycles)
        if !isempty(cycles)
            buffer = IOBuffer()
            for cycle in cycles
                println(buffer, cycle)
            end
            msg = String(take!(buffer))
            @warn "Cycles in parameters:\n$msg"
        end
    end
    evaluate_varmap!(op, ps; limit = substitution_limit)
    if is_split(sys)
        p = MTKParameters(sys, op; floatT = floatT)
    else
        p = better_varmap_to_vars(op, ps; tofloat, container_type = pType)
    end

    if implicit_dae && du0map !== nothing
        ddvs = map(Differential(iv), dvs)
        du0map = to_varmap(du0map, ddvs)
        merge!(op, du0map)
        du0 = varmap_to_vars(op, ddvs; toterm = identity,
            tofloat)
        kwargs = merge(kwargs, (; ddvs))
    else
        du0 = nothing
    end

    if build_initializeprob
        t0 = t
        if is_time_dependent(sys) && t0 === nothing
            t0 = zero(floatT)
        end
        initialization_data = SciMLBase.remake_initialization_data(
            kwargs.initialization_data, kwargs, u0, t0, p, u0, p)
        kwargs = merge(kwargs,)
    end

    f = constructor(sys, dvs, ps, u0; p = p,
        eval_expression = eval_expression,
        eval_module = eval_module,
        kwargs...)
    implicit_dae ? (f, du0, u0, p) : (f, u0, p)
end

# Check that the keys of a u0map or pmap are valid
# (i.e. are symbolic keys, and are defined for the system.)
function check_inputmap_keys(sys, u0map, pmap)
    badvarkeys = Any[]
    for k in keys(u0map)
        if symbolic_type(k) === NotSymbolic()
            push!(badvarkeys, k)
        end
    end

    badparamkeys = Any[]
    for k in keys(pmap)
        if symbolic_type(k) === NotSymbolic()
            push!(badparamkeys, k)
        end
    end
    (isempty(badvarkeys) && isempty(badparamkeys)) ||
        throw(InvalidKeyError(collect(badvarkeys), collect(badparamkeys)))
end

const BAD_KEY_MESSAGE = """
                        Undefined keys found in the parameter or initial condition maps. Check if symbolic variable names have been reassigned. 
                        The following keys are invalid:
                        """

struct InvalidKeyError <: Exception
    vars::Any
    params::Any
end

function Base.showerror(io::IO, e::InvalidKeyError)
    println(io, BAD_KEY_MESSAGE)
    println(io, "u0map: $(join(e.vars, ", "))")
    println(io, "pmap: $(join(e.params, ", "))")
end

function SciMLBase.detect_cycles(sys::AbstractSystem, varmap::Dict{Any, Any}, vars)
    varmap = AnyDict(unwrap(k) => unwrap(v) for (k, v) in varmap)
    vars = map(unwrap, vars)
    cycles = check_substitution_cycles(varmap, vars)
    return !isempty(cycles)
end

##############
# Legacy functions for backward compatibility
##############

"""
    u0, p, defs = get_u0_p(sys, u0map, parammap; use_union=true, tofloat=true)

Take dictionaries with initial conditions and parameters and convert them to numeric arrays `u0` and `p`. Also return the merged dictionary `defs` containing the entire operating point.
"""
function get_u0_p(sys,
        u0map,
        parammap = nothing;
        t0 = nothing,
        tofloat = true,
        use_union = true,
        symbolic_u0 = false)
    dvs = unknowns(sys)
    ps = parameters(sys; initial_parameters = true)

    defs = defaults(sys)
    if t0 !== nothing
        defs[get_iv(sys)] = t0
    end
    if parammap !== nothing
        defs = mergedefaults(defs, parammap, ps)
    end
    if u0map isa Vector && eltype(u0map) <: Pair
        u0map = Dict(u0map)
    end
    if u0map isa Dict
        allobs = Set(observables(sys))
        if any(in(allobs), keys(u0map))
            u0s_in_obs = filter(in(allobs), keys(u0map))
            @warn "Observed variables cannot be assigned initial values. Initial values for $u0s_in_obs will be ignored."
        end
    end
    obs = filter!(x -> !(x[1] isa Number), map(x -> x.rhs => x.lhs, observed(sys)))
    observedmap = isempty(obs) ? Dict() : todict(obs)
    defs = mergedefaults(defs, observedmap, u0map, dvs)
    for (k, v) in defs
        if Symbolics.isarraysymbolic(k)
            ks = scalarize(k)
            length(ks) == length(v) || error("$k has default value $v with unmatched size")
            for (kk, vv) in zip(ks, v)
                if !haskey(defs, kk)
                    defs[kk] = vv
                end
            end
        end
    end

    if symbolic_u0
        u0 = varmap_to_vars(u0map, dvs; defaults = defs, tofloat = false, use_union = false)
    else
        u0 = varmap_to_vars(u0map, dvs; defaults = defs, tofloat, use_union)
    end
    p = varmap_to_vars(parammap, ps; defaults = defs, tofloat, use_union)
    p = p === nothing ? SciMLBase.NullParameters() : p
    t0 !== nothing && delete!(defs, get_iv(sys))
    u0, p, defs
end

function get_u0(
        sys, u0map, parammap = nothing; symbolic_u0 = false,
        toterm = default_toterm, t0 = nothing, use_union = true)
    dvs = unknowns(sys)
    ps = parameters(sys)
    defs = defaults(sys)
    if t0 !== nothing
        defs[get_iv(sys)] = t0
    end
    if parammap !== nothing
        defs = mergedefaults(defs, parammap, ps)
    end

    # Convert observed equations "lhs ~ rhs" into defaults.
    # Use the order "lhs => rhs" by default, but flip it to "rhs => lhs"
    # if "lhs" is known by other means (parameter, another default, ...)
    # TODO: Is there a better way to determine which equations to flip?
    obs = map(x -> x.lhs => x.rhs, observed(sys))
    obs = map(x -> x[1] in keys(defs) ? reverse(x) : x, obs)
    obs = filter!(x -> !(x[1] isa Number), obs) # exclude e.g. "0 => x^2 + y^2 - 25"
    obsmap = isempty(obs) ? Dict() : todict(obs)

    defs = mergedefaults(defs, obsmap, u0map, dvs)
    if symbolic_u0
        u0 = varmap_to_vars(
            u0map, dvs; defaults = defs, tofloat = false, use_union = false, toterm)
    else
        u0 = varmap_to_vars(u0map, dvs; defaults = defs, tofloat = true, use_union, toterm)
    end
    t0 !== nothing && delete!(defs, get_iv(sys))
    return u0, defs
end
