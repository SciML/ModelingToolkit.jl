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

Given a variable-value mapping, add mappings for the `toterm` of each of the keys. `replace` controls whether
the old value should be removed.
"""
function add_toterms!(varmap::AbstractDict; toterm = default_toterm, replace = false)
    for k in collect(keys(varmap))
        ttk = toterm(unwrap(k))
        haskey(varmap, ttk) && continue
        varmap[ttk] = varmap[k]
        !isequal(k, ttk) && replace && delete!(varmap, k)
    end
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Out-of-place version of [`add_toterms!`](@ref).
"""
function add_toterms(varmap::AbstractDict; kwargs...)
    cp = copy(varmap)
    add_toterms!(cp; kwargs...)
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

const MISSING_VARIABLES_MESSAGE = """
                                Initial condition underdefined. Some are missing from the variable map.
                                Please provide a default (`u0`), initialization equation, or guess
                                for the following variables:
                                """

struct MissingVariablesError <: Exception
    vars::Any
end

function Base.showerror(io::IO, e::MissingVariablesError)
    println(io, MISSING_VARIABLES_MESSAGE)
    println(io, join(e.vars, ", "))
end

"""
    $(TYPEDSIGNATURES)

Return an array of values where the `i`th element corresponds to the value of `vars[i]`
in `varmap`. Will mutate `varmap` by symbolically substituting it into itself.

Keyword arguments:
- `container_type`: The type of the returned container.
- `allow_symbolic`: Whether the returned container of values can have symbolic expressions.
- `buffer_eltype`: The `eltype` of the returned container if `!allow_symbolic`. If
  `Nothing`, automatically promotes the values in the container to a common `eltype`.
- `tofloat`: Whether to promote values to floating point numbers if
  `buffer_eltype == Nothing`.
- `use_union`: Whether to allow using a `Union` as the `eltype` if
  `buffer_eltype == Nothing`.
- `toterm`: The `toterm` function for canonicalizing keys of `varmap`. A value of `nothing`
  disables this process.
- `check`: Whether to check if all of `vars` are keys of `varmap`.
- `is_initializeprob`: Whether an initialization problem is being constructed. Used for
  better error messages.
- `substitution_limit`: The maximum number of times to recursively substitute `varmap` into
  itself to get a numeric value for each variable in `vars`.
"""
function varmap_to_vars(varmap::AbstractDict, vars::Vector;
        tofloat = true, use_union = false, container_type = Array, buffer_eltype = Nothing,
        toterm = default_toterm, check = true, allow_symbolic = false,
        is_initializeprob = false, substitution_limit = 100)
    isempty(vars) && return nothing

    varmap = recursive_unwrap(varmap)
    if toterm !== nothing
        add_toterms!(varmap; toterm)
    end
    if check && !allow_symbolic
        missing_vars = missingvars(varmap, vars; toterm)
        if !isempty(missing_vars)
            if is_initializeprob
                throw(MissingGuessError(collect(missing_vars), collect(missing_vars)))
            else
                throw(MissingVariablesError(missing_vars))
            end
        end
    end
    evaluate_varmap!(varmap, vars; limit = substitution_limit)
    vals = map(x -> get(varmap, x, x), vars)
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
        if buffer_eltype == Nothing
            vals = promote_to_concrete(vals; tofloat, use_union)
        else
            vals = Vector{buffer_eltype}(vals)
        end
    end

    if container_type <: Union{AbstractDict, Nothing, SciMLBase.NullParameters}
        container_type = Array
    end

    if isempty(vals)
        return nothing
    elseif container_type <: Tuple
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
        v = get(varmap, k, nothing)
        v === nothing && continue
        symbolic_type(v) == NotSymbolic() && !is_array_of_symbolics(v) && continue
        haskey(varmap, k) || continue
        varmap[k] = fixpoint_sub(v, varmap; maxiters = limit)
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

"""
    $(TYPEDSIGNATURES)

For each array variable in `vars`, scalarize the corresponding entry in `varmap`.
If a scalarized entry already exists, it is not overridden.
"""
function scalarize_vars_in_varmap!(varmap::AbstractDict, vars)
    for var in vars
        symbolic_type(var) == ArraySymbolic() || continue
        is_sized_array_symbolic(var) || continue
        haskey(varmap, var) || continue
        for i in eachindex(var)
            haskey(varmap, var[i]) && continue
            varmap[var[i]] = varmap[var][i]
        end
    end
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
struct EmptySciMLFunction{iip, A, K} <: SciMLBase.AbstractSciMLFunction{iip}
    args::A
    kwargs::K
end

function EmptySciMLFunction{iip}(args...; kwargs...) where {iip}
    return EmptySciMLFunction{iip, typeof(args), typeof(kwargs)}(args, kwargs)
end

"""
    $(TYPEDSIGNATURES)

Construct the operating point of the system from the user-provided `u0map` and `pmap`, system
defaults `defs`, unknowns `dvs` and parameters `ps`. Return the operating point as a dictionary,
the list of unknowns for which no values can be determined, and the list of parameters for which
no values can be determined.

Also updates `u0map` and `pmap` in-place to contain all the initial conditions in `op`, split
by unknowns and parameters respectively.
"""
function build_operating_point!(sys::AbstractSystem,
        op::AbstractDict, u0map::AbstractDict, pmap::AbstractDict, defs::AbstractDict, dvs, ps)
    add_toterms!(op)
    missing_unknowns = add_fallbacks!(op, dvs, defs)
    for (k, v) in defs
        haskey(op, k) && continue
        op[k] = v
    end
    filter_missing_values!(op; missing_values = missing_unknowns)

    merge!(op, pmap)
    missing_pars = add_fallbacks!(op, ps, defs)
    filter_missing_values!(op; missing_values = missing_pars)

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

    if !isempty(neithermap)
        for (k, v) in u0map
            symbolic_type(v) == NotSymbolic() && !is_array_of_symbolics(v) && continue
            v = fixpoint_sub(v, neithermap; operator = Symbolics.Operator)
            isequal(k, v) && continue
            u0map[k] = v
        end
        for (k, v) in pmap
            symbolic_type(v) == NotSymbolic() && !is_array_of_symbolics(v) && continue
            v = fixpoint_sub(v, neithermap; operator = Symbolics.Operator)
            isequal(k, v) && continue
            pmap[k] = v
        end
    end

    return missing_unknowns, missing_pars
end

"""
    $(TYPEDEF)

A callable struct used to reconstruct the `u0` and `p` of the initialization problem
with promoted types.

# Fields

$(TYPEDFIELDS)
"""
struct ReconstructInitializeprob{GP, GU}
    """
    A function which when given the original problem and initialization problem, returns
    the parameter object of the initialization problem with values copied from the
    original.
    """
    pgetter::GP
    """
    Given the original problem, return the `u0` of the initialization problem.
    """
    ugetter::GU
end

"""
    $(TYPEDEF)

A wrapper over an observed function which allows calling it on a problem-like object.
`TD` determines whether the getter function is `(u, p, t)` (if `true`) or `(u, p)` (if
`false`).
"""
struct ObservedWrapper{TD, F}
    f::F
end

ObservedWrapper{TD}(f::F) where {TD, F} = ObservedWrapper{TD, F}(f)

function (ow::ObservedWrapper{true})(prob)
    ow.f(state_values(prob), parameter_values(prob), current_time(prob))
end

function (ow::ObservedWrapper{false})(prob)
    ow.f(state_values(prob), parameter_values(prob))
end

"""
    $(TYPEDSIGNATURES)

Given an index provider `indp` and a vector of symbols `syms` return a type-stable getter
function.

Note that the getter ONLY works for problem-like objects, since it generates an observed
function. It does NOT work for solutions.
"""
Base.@nospecializeinfer function concrete_getu(indp, syms; eval_expression, eval_module)
    @nospecialize
    obsfn = build_explicit_observed_function(
        indp, syms; wrap_delays = false, eval_expression, eval_module)
    return ObservedWrapper{is_time_dependent(indp)}(obsfn)
end

"""
    $(TYPEDEF)

A callable struct which applies `p_constructor` to possibly nested arrays. It also
ensures that views (including nested ones) are concretized. This is implemented manually
of using `narrow_buffer_type` to preserve type-stability.
"""
struct PConstructorApplicator{F}
    p_constructor::F
end

function (pca::PConstructorApplicator)(x::AbstractArray)
    pca.p_constructor(x)
end

function (pca::PConstructorApplicator)(x::AbstractArray{Bool})
    pca.p_constructor(BitArray(x))
end

function (pca::PConstructorApplicator{typeof(identity)})(x::SubArray)
    collect(x)
end

function (pca::PConstructorApplicator{typeof(identity)})(x::SubArray{Bool})
    BitArray(x)
end

function (pca::PConstructorApplicator{typeof(identity)})(x::SubArray{<:AbstractArray})
    collect(pca.(x))
end

function (pca::PConstructorApplicator)(x::AbstractArray{<:AbstractArray})
    pca.p_constructor(pca.(x))
end

"""
    $(TYPEDSIGNATURES)

Given a source system `srcsys` and destination system `dstsys`, return a function that
takes a value provider of `srcsys` and a value provider of `dstsys` and returns the
`MTKParameters` object of the latter with values from the former.

# Keyword Arguments
- `initials`: Whether to include the `Initial` parameters of `dstsys` among the values
  to be transferred.
- `unwrap_initials`: Whether initials in `dstsys` corresponding to unknowns in `srcsys` are
  unwrapped.
- `p_constructor`: The `p_constructor` argument to `process_SciMLProblem`.
"""
function get_mtkparameters_reconstructor(srcsys::AbstractSystem, dstsys::AbstractSystem;
        initials = false, unwrap_initials = false, p_constructor = identity,
        eval_expression = false, eval_module = @__MODULE__)
    _p_constructor = p_constructor
    p_constructor = PConstructorApplicator(p_constructor)
    # if we call `getu` on this (and it were able to handle empty tuples) we get the
    # fields of `MTKParameters` except caches.
    syms = reorder_parameters(
        dstsys, parameters(dstsys; initial_parameters = initials); flatten = false)
    # `dstsys` is an initialization system, do basically everything is a tunable
    # and tunables are a mix of different types in `srcsys`. No initials. Constants
    # are going to be constants in `srcsys`, as are `nonnumeric`.

    # `syms[1]` is always the tunables because `srcsys` will have initials.
    tunable_syms = syms[1]
    tunable_getter = if isempty(tunable_syms)
        Returns(SizedVector{0, Float64}())
    else
        p_constructor ∘ concrete_getu(srcsys, tunable_syms; eval_expression, eval_module)
    end
    initials_getter = if initials && !isempty(syms[2])
        initsyms = Vector{Any}(syms[2])
        allsyms = Set(variable_symbols(srcsys))
        if unwrap_initials
            for i in eachindex(initsyms)
                sym = initsyms[i]
                innersym = if operation(sym) === getindex
                    sym, idxs... = arguments(sym)
                    only(arguments(sym))[idxs...]
                else
                    only(arguments(sym))
                end
                if innersym in allsyms
                    initsyms[i] = innersym
                end
            end
        end
        p_constructor ∘ concrete_getu(srcsys, initsyms; eval_expression, eval_module)
    else
        Returns(SizedVector{0, Float64}())
    end
    discs_getter = if isempty(syms[3])
        Returns(())
    else
        ic = get_index_cache(dstsys)
        blockarrsizes = Tuple(map(ic.discrete_buffer_sizes) do bufsizes
            p_constructor(map(x -> x.length, bufsizes))
        end)
        # discretes need to be blocked arrays
        # the `getu` returns a tuple of arrays corresponding to `p.discretes`
        # `Base.Fix1(...)` applies `p_constructor` to each of the arrays in the tuple
        # `Base.Fix2(...)` does `BlockedArray.(tuple_of_arrs, blockarrsizes)` returning a
        # tuple of `BlockedArray`s
        Base.Fix2(Broadcast.BroadcastFunction(BlockedArray), blockarrsizes) ∘
        Base.Fix1(broadcast, p_constructor) ∘
        concrete_getu(srcsys, syms[3]; eval_expression, eval_module)
    end
    const_getter = if syms[4] == ()
        Returns(())
    else
        Base.Fix1(broadcast, p_constructor) ∘ getu(srcsys, syms[4])
    end
    nonnumeric_getter = if syms[5] == ()
        Returns(())
    else
        ic = get_index_cache(dstsys)
        buftypes = Tuple(map(ic.nonnumeric_buffer_sizes) do bufsize
            Vector{bufsize.type}
        end)
        # nonnumerics retain the assigned buffer type without narrowing
        Base.Fix1(broadcast, _p_constructor) ∘
        Base.Fix1(Broadcast.BroadcastFunction(call), buftypes) ∘
        concrete_getu(srcsys, syms[5]; eval_expression, eval_module)
    end
    getters = (
        tunable_getter, initials_getter, discs_getter, const_getter, nonnumeric_getter)
    getter = let getters = getters
        function _getter(valp, initprob)
            oldcache = parameter_values(initprob).caches
            MTKParameters(getters[1](valp), getters[2](valp), getters[3](valp),
                getters[4](valp), getters[5](valp), oldcache isa Tuple{} ? () :
                                                    copy.(oldcache))
        end
    end

    return getter
end

function call(f, args...)
    f(args...)
end

"""
    $(TYPEDSIGNATURES)

Construct a `ReconstructInitializeprob` which reconstructs the `u0` and `p` of `dstsys`
with values from `srcsys`.
"""
function ReconstructInitializeprob(
        srcsys::AbstractSystem, dstsys::AbstractSystem; u0_constructor = identity, p_constructor = identity,
        eval_expression = false, eval_module = @__MODULE__)
    @assert is_initializesystem(dstsys)
    ugetter = u0_constructor ∘
              concrete_getu(srcsys, unknowns(dstsys); eval_expression, eval_module)
    if is_split(dstsys)
        pgetter = get_mtkparameters_reconstructor(
            srcsys, dstsys; p_constructor, eval_expression, eval_module)
    else
        syms = parameters(dstsys)
        pgetter = let inner = concrete_getu(srcsys, syms; eval_expression, eval_module),
            p_constructor = p_constructor

            function _getter2(valp, initprob)
                p_constructor(inner(valp))
            end
        end
    end
    return ReconstructInitializeprob(pgetter, ugetter)
end

"""
    $(TYPEDSIGNATURES)

Copy values from `srcvalp` to `dstvalp`. Returns the new `u0` and `p`.
"""
function (rip::ReconstructInitializeprob)(srcvalp, dstvalp)
    # copy parameters
    newp = rip.pgetter(srcvalp, dstvalp)
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
    u0 = rip.ugetter(srcvalp)
    # and the eltype of the destination u0
    if T != eltype(u0) && T != Union{}
        u0 = T.(u0)
    end
    # apply the promotion to tunables portion
    buf, repack, alias = SciMLStructures.canonicalize(SciMLStructures.Tunable(), newp)
    if eltype(buf) != T
        # only do a copy if the eltype doesn't match
        newbuf = similar(buf, T)
        copyto!(newbuf, buf)
        newp = repack(newbuf)
    end
    if newp isa MTKParameters
        # and initials portion
        buf, repack, alias = SciMLStructures.canonicalize(SciMLStructures.Initials(), newp)
        if eltype(buf) != T
            newbuf = similar(buf, T)
            copyto!(newbuf, buf)
            newp = repack(newbuf)
        end
    end
    return u0, newp
end

"""
    $(TYPEDSIGNATURES)

Given `sys` and its corresponding initialization system `initsys`, return the
`initializeprobpmap` function in `OverrideInitData` for the systems.
"""
function construct_initializeprobpmap(
        sys::AbstractSystem, initsys::AbstractSystem; p_constructor = identity, eval_expression, eval_module)
    @assert is_initializesystem(initsys)
    if is_split(sys)
        return let getter = get_mtkparameters_reconstructor(
                initsys, sys; initials = true, unwrap_initials = true, p_constructor,
                eval_expression, eval_module)
            function initprobpmap_split(prob, initsol)
                getter(initsol, prob)
            end
        end
    else
        return let getter = concrete_getu(
                initsys, parameters(sys; initial_parameters = true);
                eval_expression, eval_module), p_constructor = p_constructor

            function initprobpmap_nosplit(prob, initsol)
                return p_constructor(getter(initsol))
            end
        end
    end
end

function get_scimlfn(valp)
    valp isa SciMLBase.AbstractSciMLFunction && return valp
    if hasmethod(symbolic_container, Tuple{typeof(valp)}) &&
       (sc = symbolic_container(valp)) !== valp
        return get_scimlfn(sc)
    end
    throw(ArgumentError("SciMLFunction not found. This should never happen."))
end

"""
    $(TYPEDSIGNATURES)

A function to be used as `update_initializeprob!` in `OverrideInitData`. Requires
`is_update_oop = Val(true)` to be passed to `update_initializeprob!`.
"""
function update_initializeprob!(initprob, prob)
    pgetter = ChainRulesCore.@ignore_derivatives get_scimlfn(prob).initialization_data.metadata.oop_reconstruct_u0_p.pgetter
    p = pgetter(prob, initprob)
    return remake(initprob; p)
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
    The operating point used to construct the initialization.
    """
    op::Dict{Any, Any}
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
    Whether the initialization uses the independent variable.
    """
    time_dependent_init::Bool
    """
    `ReconstructInitializeprob` for this initialization problem.
    """
    oop_reconstruct_u0_p::R
    """
    A function which takes `(prob, initializeprob)` and return the `u0` to use for the problem.
    """
    get_updated_u0::GUU
    """
    A function which takes parameter object and `u0` of the problem and sets
    `Initial.(unknowns(sys))` in the former, returning the updated parameter object.
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

function GetUpdatedU0(sys::AbstractSystem, initprob::SciMLBase.AbstractNonlinearProblem, op::AbstractDict)
    dvs = unknowns(sys)
    eqs = equations(sys)
    guessvars = trues(length(dvs))
    for (i, var) in enumerate(dvs)
        guessvars[i] = !isequal(get(op, var, nothing), Initial(var))
    end
    get_guessvars = getu(initprob, dvs[guessvars])
    get_initial_unknowns = getu(sys, Initial.(dvs))
    return GetUpdatedU0(guessvars, get_guessvars, get_initial_unknowns)
end

function (guu::GetUpdatedU0)(prob, initprob)
    buffer = guu.get_initial_unknowns(prob)
    algebuf = view(buffer, guu.guessvars)
    copyto!(algebuf, guu.get_guessvars(initprob))
    return buffer
end

struct SetInitialUnknowns{S}
    setter!::S
end

function SetInitialUnknowns(sys::AbstractSystem)
    return SetInitialUnknowns(setu(sys, Initial.(unknowns(sys))))
end

function (siu::SetInitialUnknowns)(p::MTKParameters, u0)
    if ArrayInterface.ismutable(p.initials)
        siu.setter!(p, u0)
    else
        originalT = similar_type(p.initials)
        @set! p.initials = MVector{length(p.initials), eltype(p.initials)}(p.initials)
        siu.setter!(p, u0)
        @set! p.initials = originalT(p.initials)
    end
    return p
end

function (siu::SetInitialUnknowns)(p::AbstractVector, u0)
    if ArrayInterface.ismutable(p)
        siu.setter!(p, u0)
    else
        originalT = similar_type(p)
        p = MVector{length(p), eltype(p)}(p)
        siu.setter!(p, u0)
        p = originalT(p)
    end
    return p
end

safe_float(x) = x
safe_float(x::AbstractArray) = isempty(x) ? x : float(x)

"""
    $(TYPEDSIGNATURES)

Build and return the initialization problem and associated data as a `NamedTuple` to be passed
to the `SciMLFunction` constructor. Requires the system `sys`, operating point `op`, initial
time `t`, system defaults `defs`, user-provided `guesses`, and list of unknowns which don't
have a value in `op`. The keyword `implicit_dae` denotes whether the `SciMLProblem` being
constructed is in implicit DAE form (`DAEProblem`). All other keyword arguments are forwarded
to `InitializationProblem`.
"""
function maybe_build_initialization_problem(
        sys::AbstractSystem, iip, op::AbstractDict, t, defs,
        guesses, missing_unknowns; implicit_dae = false,
        time_dependent_init = is_time_dependent(sys), u0_constructor = identity,
        p_constructor = identity, floatT = Float64, initialization_eqs = [],
        use_scc = true, eval_expression = false, eval_module = @__MODULE__, kwargs...)
    guesses = merge(ModelingToolkit.guesses(sys), todict(guesses))

    if t === nothing && is_time_dependent(sys)
        t = zero(floatT)
    end

    initializeprob = ModelingToolkit.InitializationProblem{iip}(
        sys, t, op; guesses, time_dependent_init, initialization_eqs,
        use_scc, u0_constructor, p_constructor, eval_expression, eval_module, kwargs...)
    if state_values(initializeprob) !== nothing
        _u0 = state_values(initializeprob)
        if ArrayInterface.ismutable(_u0)
            _u0 = floatT.(_u0)
        else
            _u0 = similar_type(_u0, floatT)(_u0)
        end
        initializeprob = remake(initializeprob; u0 = _u0)
    end
    initp = parameter_values(initializeprob)
    if is_split(sys)
        buffer, repack, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), initp)
        initp = repack(floatT.(buffer))
        buffer, repack, _ = SciMLStructures.canonicalize(SciMLStructures.Initials(), initp)
        initp = repack(floatT.(buffer))
    elseif initp isa AbstractArray
        if ArrayInterface.ismutable(initp)
            initp′ = similar(initp, floatT)
            copyto!(initp′, initp)
            initp = initp′
        else
            initp = similar_type(initp, floatT)(initp)
        end
    end
    initializeprob = remake(initializeprob; p = initp)

    get_initial_unknowns = if time_dependent_init
        GetUpdatedU0(sys, initializeprob, op)
    else
        nothing
    end
    meta = InitializationMetadata(
        copy(op), copy(guesses), Vector{Equation}(initialization_eqs),
        use_scc, time_dependent_init,
        ReconstructInitializeprob(
            sys, initializeprob.f.sys; u0_constructor,
            p_constructor, eval_expression, eval_module),
        get_initial_unknowns, SetInitialUnknowns(sys))

    if time_dependent_init
        all_init_syms = Set(all_symbols(initializeprob))
        solved_unknowns = filter(var -> var in all_init_syms, unknowns(sys))
        initializeprobmap = u0_constructor ∘ safe_float ∘
                            getu(initializeprob, solved_unknowns)
    else
        initializeprobmap = nothing
    end

    punknowns = [p
                 for p in all_variable_symbols(initializeprob)
                 if is_parameter(sys, p)]
    if initializeprobmap === nothing && isempty(punknowns)
        initializeprobpmap = nothing
    else
        initializeprobpmap = construct_initializeprobpmap(
            sys, initializeprob.f.sys; p_constructor, eval_expression, eval_module)
    end

    # we still want the `initialization_data` because it helps with `remake`
    if initializeprobmap === nothing && initializeprobpmap === nothing
        update_initializeprob! = nothing
    else
        update_initializeprob! = ModelingToolkit.update_initializeprob!
    end

    filter!(punknowns) do p
        is_parameter_solvable(p, op, defs, guesses) && get(op, p, missing) === missing
    end
    # See comment below for why `getu` is not used here.
    _pgetter = build_explicit_observed_function(initializeprob.f.sys, punknowns)
    pvals = _pgetter(state_values(initializeprob), parameter_values(initializeprob))
    for (p, pval) in zip(punknowns, pvals)
        p = unwrap(p)
        op[p] = pval
        if iscall(p) && operation(p) === getindex
            arrp = arguments(p)[1]
            get(op, arrp, nothing) !== missing && continue
            op[arrp] = collect(arrp)
        end
    end

    if time_dependent_init
        # We can't use `getu` here because that goes to `SII.observed`, which goes to
        # `ObservedFunctionCache` which uses `eval_expression` and `eval_module`. If
        # `eval_expression == true`, this then runs into world-age issues. Building an
        # RGF here is fine since it is always discarded. We can't use `eval_module` for
        # the RGF since the user may not have run RGF's init.
        _ugetter = build_explicit_observed_function(initializeprob.f.sys, collect(missing_unknowns))
        uvals = _ugetter(state_values(initializeprob), parameter_values(initializeprob))
        for (v, val) in zip(missing_unknowns, uvals)
            op[v] = val
        end
        empty!(missing_unknowns)
    end

    return (;
        initialization_data = SciMLBase.OverrideInitData(
        initializeprob, update_initializeprob!, initializeprobmap,
        initializeprobpmap; metadata = meta, is_update_oop = Val(true)))
end

"""
    $(TYPEDSIGNATURES)

Calculate the floating point type to use from the given `varmap` by looking at variables
with a constant value.
"""
function float_type_from_varmap(varmap, floatT = Bool)
    for (k, v) in varmap
        is_variable_floatingpoint(k) || continue
        symbolic_type(v) == NotSymbolic() || continue
        is_array_of_symbolics(v) && continue

        if v isa AbstractArray
            floatT = promote_type(floatT, eltype(v))
        elseif v isa Number
            floatT = promote_type(floatT, typeof(v))
        end
    end
    return float(floatT)
end

"""
    $(TYPEDSIGNATURES)

Calculate the floating point type to use from the given `varmap` by looking at variables
with a constant value. `u0Type` takes priority if it is a real-valued array type.
"""
function calculate_float_type(varmap, u0Type::Type, floatT = Bool)
    if u0Type <: AbstractArray && eltype(u0Type) <: Real && eltype(u0Type) != Union{}
        return float(eltype(u0Type))
    else
        return float_type_from_varmap(varmap, floatT)
    end
end

"""
    $(TYPEDSIGNATURES)

Calculate the `resid_prototype` for a `NonlinearFunction` with `N` equations and the
provided `u0` and `p`.
"""
function calculate_resid_prototype(N::Int, u0, p)
    u0ElType = u0 === nothing ? Float64 : eltype(u0)
    if SciMLStructures.isscimlstructure(p)
        u0ElType = promote_type(
            eltype(SciMLStructures.canonicalize(SciMLStructures.Tunable(), p)[1]),
            u0ElType)
    end
    return zeros(u0ElType, N)
end

"""
    $(TYPEDSIGNATURES)

Given the user-provided value of `u0_constructor`, the container type of user-provided
`op`, the desired floating point type and whether a symbolic `u0` is allowed, return the
updated `u0_constructor`.
"""
function get_u0_constructor(u0_constructor, u0Type::Type, floatT::Type, symbolic_u0::Bool)
    u0_constructor === identity || return u0_constructor
    u0Type <: StaticArray || return u0_constructor
    return function (vals)
        elT = if symbolic_u0 && any(x -> symbolic_type(x) != NotSymbolic(), vals)
            nothing
        else
            floatT
        end
        SymbolicUtils.Code.create_array(u0Type, elT, Val(1), Val(length(vals)), vals...)
    end
end

"""
    $(TYPEDSIGNATURES)

Given the user-provided value of `p_constructor`, the container type of user-provided `op`,
ans the desired floating point type, return the updated `p_constructor`.
"""
function get_p_constructor(p_constructor, pType::Type, floatT::Type)
    p_constructor === identity || return p_constructor
    pType <: StaticArray || return p_constructor
    return function (vals)
        SymbolicUtils.Code.create_array(
            pType, floatT, Val(ndims(vals)), Val(size(vals)), vals...)
    end
end

abstract type ProblemConstructionHook end

"""
    $(TYPEDSIGNATURES)

Return the SciMLFunction created via calling `constructor`, the initial conditions `u0`
and parameter object `p` given the system `sys`, and user-provided initial values `u0map`
and `pmap`. `u0map` and `pmap` are converted into variable maps via [`to_varmap`](@ref).

$U0_P_DOCS

This will also build the initialization problem and related objects and pass them to the
SciMLFunction as keyword arguments.

Keyword arguments:
$PROBLEM_KWARGS
$PROBLEM_INTERNAL_KWARGS
- `t`: The initial time of the `SciMLProblem`. This does not need to be provided for time-
  independent problems. If not provided for time-dependent problems, will be assumed as
  zero.
- `implicit_dae`: Also build a mapping of derivatives of states to values for implicit DAEs.
  Changes the return value of this function to `(f, du0, u0, p)` instead of `(f, u0, p)`.
- `symbolic_u0` allows the returned `u0` to be an array of symbolics.

All other keyword arguments are passed as-is to `constructor`.
"""
function process_SciMLProblem(
        constructor, sys::AbstractSystem, op;
        build_initializeprob = supports_initialization(sys),
        implicit_dae = false, t = nothing, guesses = AnyDict(),
        warn_initialize_determined = true, initialization_eqs = [],
        eval_expression = false, eval_module = @__MODULE__, fully_determined = nothing,
        check_initialization_units = false, u0_eltype = nothing, tofloat = true,
        u0_constructor = identity, p_constructor = identity,
        check_length = true, symbolic_u0 = false, warn_cyclic_dependency = false,
        circular_dependency_max_cycle_length = length(all_symbols(sys)),
        circular_dependency_max_cycles = 10,
        substitution_limit = 100, use_scc = true, time_dependent_init = is_time_dependent(sys),
        algebraic_only = false,
        allow_incomplete = false, is_initializeprob = false, kwargs...)
    dvs = unknowns(sys)
    ps = parameters(sys; initial_parameters = true)
    iv = has_iv(sys) ? get_iv(sys) : nothing
    eqs = equations(sys)

    check_array_equations_unknowns(eqs, dvs)

    u0Type = pType = typeof(op)

    op = to_varmap(op, dvs)
    symbols_to_symbolics!(sys, op)

    check_inputmap_keys(sys, op)

    op = getmetadata(sys, ProblemConstructionHook, identity)(op)

    defs = add_toterms(recursive_unwrap(defaults(sys)); replace = is_discrete_system(sys))
    kwargs = NamedTuple(kwargs)

    if eltype(eqs) <: Equation
        obs, eqs = unhack_observed(observed(sys), eqs)
    else
        obs, _ = unhack_observed(observed(sys), Equation[x for x in eqs if x isa Equation])
    end

    u0map = anydict()
    pmap = anydict()
    missing_unknowns,
    missing_pars = build_operating_point!(sys, op,
        u0map, pmap, defs, dvs, ps)

    floatT = calculate_float_type(op, u0Type)
    u0_eltype = something(u0_eltype, floatT)

    if !is_time_dependent(sys) || is_initializesystem(sys)
        add_observed_equations!(op, obs)
    end

    u0_constructor = get_u0_constructor(u0_constructor, u0Type, u0_eltype, symbolic_u0)
    p_constructor = get_p_constructor(p_constructor, pType, floatT)

    if build_initializeprob
        kws = maybe_build_initialization_problem(
            sys, constructor <: SciMLBase.AbstractSciMLFunction{true},
            op, t, defs, guesses, missing_unknowns;
            implicit_dae, warn_initialize_determined, initialization_eqs,
            eval_expression, eval_module, fully_determined,
            warn_cyclic_dependency, check_units = check_initialization_units,
            circular_dependency_max_cycle_length, circular_dependency_max_cycles, use_scc,
            algebraic_only, allow_incomplete, u0_constructor, p_constructor, floatT,
            time_dependent_init)

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

    u0 = varmap_to_vars(
        op, dvs; buffer_eltype = u0_eltype, container_type = u0Type,
        allow_symbolic = symbolic_u0, is_initializeprob, substitution_limit)

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

    if is_split(sys)
        # `pType` is usually `Dict` when the user passes key-value pairs.
        if !(pType <: AbstractArray)
            pType = Array
        end
        p = MTKParameters(sys, op; floatT = floatT, p_constructor, fast_path = true)
    else
        p = p_constructor(varmap_to_vars(op, ps; tofloat, container_type = pType))
    end

    if implicit_dae
        ddvs = map(Differential(iv), dvs)
        du0 = varmap_to_vars(op, ddvs; toterm = default_toterm,
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
        initialization_data = @invokelatest SciMLBase.remake_initialization_data(
            sys, kwargs, u0, t0, p, u0, p)
        kwargs = merge(kwargs, (; initialization_data))
    end

    if constructor <: NonlinearFunction && length(dvs) != length(eqs)
        kwargs = merge(kwargs,
            (;
                resid_prototype = u0_constructor(calculate_resid_prototype(
                length(eqs), u0, p))))
    end

    f = constructor(sys; u0 = u0, p = p,
        eval_expression = eval_expression,
        eval_module = eval_module,
        kwargs...)
    implicit_dae ? (f, du0, u0, p) : (f, u0, p)
end

# Check that the keys of a u0map or pmap are valid
# (i.e. are symbolic keys, and are defined for the system.)
function check_inputmap_keys(sys, op)
    badvarkeys = Any[]
    for k in keys(op)
        if symbolic_type(k) === NotSymbolic()
            push!(badvarkeys, k)
        end
    end

    if !isempty(badvarkeys)
        throw(InvalidKeyError(collect(badvarkeys)))
    end
end

const BAD_KEY_MESSAGE = """
                        Undefined keys found in the parameter or initial condition maps. Check if symbolic variable names have been reassigned.
                        The following keys are invalid:
                        """

struct InvalidKeyError <: Exception
    vars::Any
end

function Base.showerror(io::IO, e::InvalidKeyError)
    println(io, BAD_KEY_MESSAGE)
    println(io, join(e.vars, ", "))
end

function SciMLBase.detect_cycles(sys::AbstractSystem, varmap::Dict{Any, Any}, vars)
    varmap = AnyDict(unwrap(k) => unwrap(v) for (k, v) in varmap)
    vars = map(unwrap, vars)
    cycles = check_substitution_cycles(varmap, vars)
    return !isempty(cycles)
end

function process_kwargs(sys::System; expression = Val{false}, callback = nothing,
        eval_expression = false, eval_module = @__MODULE__, kwargs...)
    kwargs = filter_kwargs(kwargs)
    kwargs1 = (;)

    if is_time_dependent(sys)
        if expression == Val{false}
            cbs = process_events(sys; callback, eval_expression, eval_module, kwargs...)
            if cbs !== nothing
                kwargs1 = merge(kwargs1, (callback = cbs,))
            end
        end

        tstops = SymbolicTstops(sys; expression, eval_expression, eval_module)
        if tstops !== nothing
            kwargs1 = merge(kwargs1, (; tstops))
        end
    end

    return merge(kwargs1, kwargs)
end

function filter_kwargs(kwargs)
    kwargs = Dict(kwargs)
    for key in keys(kwargs)
        key in DiffEqBase.allowedkeywords || delete!(kwargs, key)
    end
    pairs(NamedTuple(kwargs))
end

struct SymbolicTstops{F}
    fn::F
end

function (st::SymbolicTstops)(p, tspan)
    buffer = reduce(vcat, st.fn(p, tspan...))
    if ArrayInterface.ismutable(buffer)
        return unique!(sort!(buffer))
    else
        return unique(sort(buffer))
    end
end

function SymbolicTstops(
        sys::AbstractSystem; expression = Val{false}, eval_expression = false,
        eval_module = @__MODULE__)
    tstops = symbolic_tstops(sys)
    isempty(tstops) && return nothing
    t0 = gensym(:t0)
    t1 = gensym(:t1)
    tstops = map(tstops) do val
        if is_array_of_symbolics(val) || val isa AbstractArray
            collect(val)
        else
            term(:, t0, unwrap(val), t1; type = AbstractArray{Real})
        end
    end
    rps = reorder_parameters(sys)
    tstops,
    _ = build_function_wrapper(sys, tstops,
        rps...,
        t0,
        t1;
        expression = Val{true},
        p_start = 1, p_end = length(rps), add_observed = false, force_SA = true)
    tstops = GeneratedFunctionWrapper{(1, 3, is_split(sys))}(
        expression, tstops, nothing; eval_expression, eval_module)

    if expression == Val{true}
        return :($SymbolicTstops($tstops))
    else
        return SymbolicTstops(tstops)
    end
end

"""
    $(TYPEDSIGNATURES)

Macro for writing problem/function constructors. Expects a function definition with type
parameters for `iip` and `specialize`. Generates fallbacks with
`specialize = SciMLBase.FullSpecialize` and `iip = true`.
"""
macro fallback_iip_specialize(ex)
    @assert Meta.isexpr(ex, :function)
    # fnname is ODEProblem{iip, spec}(args...) where {iip, spec}
    # body is function body
    fnname, body = ex.args
    @assert Meta.isexpr(fnname, :where)
    # fnname_call is ODEProblem{iip, spec}(args...)
    # where_args are `iip, spec`
    fnname_call, where_args... = fnname.args
    @assert length(where_args) == 2
    iiparg, specarg = where_args

    @assert Meta.isexpr(fnname_call, :call)
    # fnname_curly is ODEProblem{iip, spec}
    fnname_curly, args... = fnname_call.args
    # the function should have keyword arguments
    @assert Meta.isexpr(args[1], :parameters)

    # arguments to call with
    call_args = map(args) do arg
        # keyword args are in `Expr(:parameters)` so any `Expr(:kw)` here
        # are optional positional arguments. Analyze `:(f(a, b = 1; k = 1, l...))`
        # to understand
        Meta.isexpr(arg, :kw) && return arg.args[1]
        return arg
    end
    call_kwargs = map(call_args[1].args) do arg
        Meta.isexpr(arg, :...) && return arg
        @assert Meta.isexpr(arg, :kw)
        return Expr(:kw, arg.args[1], arg.args[1])
    end
    call_args[1] = Expr(:parameters, call_kwargs...)

    @assert Meta.isexpr(fnname_curly, :curly)
    # fnname_name is `ODEProblem`
    # curly_args is `iip, spec`
    fnname_name, curly_args... = fnname_curly.args
    @assert curly_args == where_args

    # callexpr_iip is `ODEProblem{iip, FullSpecialize}(call_args...)`
    callexpr_iip = Expr(
        :call, Expr(:curly, fnname_name, curly_args[1], SciMLBase.FullSpecialize), call_args...)
    # `ODEProblem{iip}`
    fnname_iip = Expr(:curly, fnname_name, curly_args[1])
    # `ODEProblem{iip}(args...)`
    fncall_iip = Expr(:call, fnname_iip, args...)
    # ODEProblem{iip}(args...) where {iip}
    fnwhere_iip = Expr(:where, fncall_iip, where_args[1])
    fn_iip = Expr(:function, fnwhere_iip, callexpr_iip)

    # `ODEProblem{true}(call_args...)`
    callexpr_base = Expr(:call, Expr(:curly, fnname_name, true), call_args...)
    # `ODEProblem(args...)`
    fncall_base = Expr(:call, fnname_name, args...)
    fn_base = Expr(:function, fncall_base, callexpr_base)

    # Handle case when this is a problem constructor and `u0map` is a `StaticArray`,
    # where `iip` should default to `false`.
    fn_sarr = nothing
    if occursin("Problem", string(fnname_name))
        # args should at least contain an argument for the `u0map`
        @assert length(args) > 2
        u0_arg = args[3]
        # should not have a type-annotation
        @assert !Meta.isexpr(u0_arg, :(::))
        if Meta.isexpr(u0_arg, :kw)
            argname, default = u0_arg.args
            u0_arg = Expr(:kw, Expr(:(::), argname, StaticArray), default)
        else
            u0_arg = Expr(:(::), u0_arg, StaticArray)
        end

        callexpr_sarr = Expr(:call, Expr(:curly, fnname_name, false), call_args...)
        fncall_sarr = Expr(:call, fnname_name, args[1], args[2], u0_arg, args[4:end]...)
        fn_sarr = Expr(:function, fncall_sarr, callexpr_sarr)
    end
    return quote
        $fn_base
        $fn_sarr
        $fn_iip
        Base.@__doc__ $ex
    end |> esc
end

"""
    $(TYPEDSIGNATURES)

Turn key-value pairs in `kws` into assignments and append them to `block.args`. `head` is
the head of the `Expr` used to create the assignment. `filter` is a function that takes the
key and returns whether or not to include it in the assignments.
"""
function namedtuple_to_assignments!(
        block, kws::NamedTuple; head = :(=), filter = Returns(true))
    for (k, v) in pairs(kws)
        filter(k) || continue
        push!(block.args, Expr(head, k, v))
    end
end

"""
    $(TYPEDSIGNATURES)

Build an expression that constructs SciMLFunction `T`. `args` is a `NamedTuple` mapping
names of positional arguments to `T` to their (expression) values. `kwargs` are parsed
as keyword arguments to the constructor.
"""
function build_scimlfn_expr(T, args::NamedTuple; kwargs...)
    kwargs = NamedTuple(kwargs)
    let_args = Expr(:block)
    namedtuple_to_assignments!(let_args, args)

    kwexpr = Expr(:parameters)
    # don't include initialization data in the generated expression
    filter = !isequal(:initialization_data)
    namedtuple_to_assignments!(let_args, kwargs; filter = filter)
    namedtuple_to_assignments!(kwexpr, kwargs; head = :kw, filter)
    let_body = Expr(:call, T, kwexpr, keys(args)...)
    return Expr(:let, let_args, let_body)
end

"""
    $(TYPEDSIGNATURES)

Build an expression that constructs SciMLProblem `T`. `args` is a `NamedTuple` mapping
names of positional arguments to `T` to their (expression) values. `kwargs` are parsed
as keyword arguments to the constructor.
"""
function build_scimlproblem_expr(T, args::NamedTuple; kwargs...)
    kwargs = NamedTuple(kwargs)
    let_args = Expr(:block)
    namedtuple_to_assignments!(let_args, args)

    kwexpr = Expr(:parameters)
    namedtuple_to_assignments!(let_args, kwargs)
    namedtuple_to_assignments!(kwexpr, kwargs; head = :kw)
    let_body = Expr(:call, remake, Expr(:call, T, kwexpr, keys(args)...))
    return Expr(:let, let_args, let_body)
end

"""
    $(TYPEDSIGNATURES)

Return an expression constructing SciMLFunction `T` with positional arguments `args`
and keywords `kwargs`.
"""
function maybe_codegen_scimlfn(::Type{Val{true}}, T, args::NamedTuple; kwargs...)
    build_scimlfn_expr(T, args; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Construct SciMLFunction `T` with positional arguments `args` and keywords `kwargs`.
"""
function maybe_codegen_scimlfn(::Type{Val{false}}, T, args::NamedTuple; kwargs...)
    T(args...; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Return an expression constructing SciMLProblem `T` with positional arguments `args`
and keywords `kwargs`.
"""
function maybe_codegen_scimlproblem(::Type{Val{true}}, T, args::NamedTuple; kwargs...)
    build_scimlproblem_expr(T, args; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Construct SciMLProblem `T` with positional arguments `args` and keywords `kwargs`.
"""
function maybe_codegen_scimlproblem(::Type{Val{false}}, T, args::NamedTuple; kwargs...)
    # Call `remake` so it runs initialization if it is trivial
    # Use `@invokelatest` to avoid world-age issues with `eval_expression = true`
    @invokelatest remake(T(args...; kwargs...))
end

"""
    $(TYPEDSIGNATURES)

Return the `u0` vector for the given system `sys` and variable-value mapping `varmap`. All
keyword arguments are forwarded to [`varmap_to_vars`](@ref).
"""
function get_u0(sys::AbstractSystem, varmap; kwargs...)
    dvs = unknowns(sys)
    ps = parameters(sys; initial_parameters = true)
    op = to_varmap(varmap, dvs)
    add_observed!(sys, op)
    add_parameter_dependencies!(sys, op)
    missing_dvs, _ = build_operating_point!(
        sys, op, Dict(), Dict(), defaults(sys), dvs, ps)

    isempty(missing_dvs) || throw(MissingVariablesError(collect(missing_dvs)))

    return varmap_to_vars(op, dvs; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Return the `p` object for the given system `sys` and variable-value mapping `varmap`. All
keyword arguments are forwarded to [`MTKParameters`](@ref) for split systems and
[`varmap_to_vars`](@ref) for non-split systems.
"""
function get_p(sys::AbstractSystem, varmap; split = is_split(sys), kwargs...)
    dvs = unknowns(sys)
    ps = parameters(sys; initial_parameters = true)
    op = to_varmap(varmap, dvs)
    add_observed!(sys, op)
    add_parameter_dependencies!(sys, op)
    _, missing_ps = build_operating_point!(
        sys, op, Dict(), Dict(), defaults(sys), dvs, ps)

    isempty(missing_ps) || throw(MissingParametersError(collect(missing_ps)))

    if split
        MTKParameters(sys, op; kwargs...)
    else
        varmap_to_vars(op, ps; kwargs...)
    end
end
