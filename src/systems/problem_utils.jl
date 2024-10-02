const AnyDict = Dict{Any, Any}

"""
    $(TYPEDSIGNATURES)

If called without arguments, return `Dict{Any, Any}`. Otherwise, interpret the input
as a symbolic map and turn it into a `Dict{Any, Any}`. Handles `SciMLBase.NullParameters`
and `nothing`.
"""
anydict() = AnyDict()
anydict(::SciMLBase.NullParameters) = AnyDict()
anydict(::Nothing) = AnyDict()
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
    for var in vars
        haskey(varmap, var) && continue
        ttvar = toterm(var)
        haskey(varmap, ttvar) && continue

        # array symbolics with a defined size may be present in the scalarized form
        if Symbolics.isarraysymbolic(var) && is_sized_array_symbolic(var)
            val = map(eachindex(var)) do idx
                # @something is lazy and saves from writing a massive if-elseif-else
                @something(get(varmap, var[idx], nothing),
                    get(varmap, ttvar[idx], nothing), get(fallbacks, var, nothing)[idx],
                    get(fallbacks, ttvar, nothing)[idx], get(fallbacks, var[idx], nothing),
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
    return anydict(unwrap(k) => unwrap(v) for (k, v) in anydict(vals))
end

"""
    $(TYPEDSIGNATURES)

Return the appropriate zero value for a symbolic variable representing a number or array of
numbers. Sized array symbolics return a zero-filled array of matching size. Unsized array
symbolics return an empty array of the appropriate `eltype`.
"""
function zero_var(x::Symbolic{T}) where {V <: Number, T <: Union{V, AbstractArray{V}}}
    if Symbolics.isarraysymbolic(x)
        if is_sized_array_symbolic(x)
            return zeros(T, size(x))
        else
            return T[]
        end
    else
        return zero(T)
    end
end

"""
    $(TYPEDSIGNATURES)

Add equations `eqs` to `varmap`. Assumes each element in `eqs` maps a single symbolic
variable to an expression representing its value. In case `varmap` already contains an
entry for `eq.lhs`, insert the reverse mapping if `eq.rhs` is not a number.
"""
function add_observed_equations!(varmap::AbstractDict, eqs)
    for eq in eqs
        if haskey(varmap, eq.lhs)
            eq.rhs isa Number && continue
            haskey(varmap, eq.rhs) && continue
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

"""
    $(TYPEDSIGNATURES)

Return an array of values where the `i`th element corresponds to the value of `vars[i]`
in `varmap`. Does not perform symbolic substitution in the values of `varmap`.

Keyword arguments:
- `tofloat`: Convert values to floating point numbers using `float`.
- `use_union`: Use a `Union`-typed array if the values have heterogeneous types.
- `container_type`: The type of container to use for the values.
- `toterm`: The `toterm` method to use for converting symbolics.
- `promotetoconcrete`: whether the promote to a concrete buffer (respecting
  `tofloat` and `use_union`). Defaults to `container_type <: AbstractArray`.
- `check`: Error if any variables in `vars` do not have a mapping in `varmap`. Uses
  [`missingvars`](@ref) to perform the check.
- `allow_symbolic` allows the returned array to contain symbolic values. If this is `true`,
  `promotetoconcrete` is set to `false`.
"""
function better_varmap_to_vars(varmap::AbstractDict, vars::Vector;
        tofloat = true, use_union = true, container_type = Array,
        toterm = default_toterm, promotetoconcrete = nothing, check = true, allow_symbolic = false)
    isempty(vars) && return nothing

    if check
        missing_vars = missingvars(varmap, vars; toterm)
        isempty(missing_vars) || throw(MissingVariablesError(missing_vars))
    end
    vals = map(x -> varmap[x], vars)

    if container_type <: Union{AbstractDict, Tuple, Nothing}
        container_type = Array
    end

    promotetoconcrete === nothing && (promotetoconcrete = container_type <: AbstractArray)
    if promotetoconcrete && !allow_symbolic
        vals = promote_to_concrete(vals; tofloat = tofloat, use_union = use_union)
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

Performs symbolic substitution on the values in `varmap`, using `varmap` itself as the
set of substitution rules. 
"""
function evaluate_varmap!(varmap::AbstractDict)
    for (k, v) in varmap
        varmap[k] = fixpoint_sub(v, varmap)
    end
end

struct GetUpdatedMTKParameters{G, S}
    # `getu` functor which gets parameters that are unknowns during initialization
    getpunknowns::G
    # `setu` functor which returns a modified MTKParameters using those parameters
    setpunknowns::S
end

function (f::GetUpdatedMTKParameters)(prob, initializesol)
    mtkp = copy(parameter_values(prob))
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

function get_temporary_value(p)
    stype = symtype(unwrap(p))
    return if stype == Real
        zero(Float64)
    elseif stype <: AbstractArray{Real}
        zeros(Float64, size(p))
    elseif stype <: Real
        zero(stype)
    elseif stype <: AbstractArray
        zeros(eltype(stype), size(p))
    else
        error("Nonnumeric parameter $p with symtype $stype cannot be solved for during initialization")
    end
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
- `check_units`: Enable or disable unit checks.
- `tofloat`, `use_union`: Passed to [`better_varmap_to_vars`](@ref) for building `u0` (and
  possibly `p`).
- `u0_constructor`: A function to apply to the `u0` value returned from `better_varmap_to_vars`
  to construct the final `u0` value.
- `du0map`: A map of derivatives to values. See `implicit_dae`.
- `check_length`: Whether to check the number of equations along with number of unknowns and
  length of `u0` vector for consistency. If `false`, do not check with equations. This is
  forwarded to `check_eqs_u0`
- `symbolic_u0` allows the returned `u0` to be an array of symbolics.

All other keyword arguments are passed as-is to `constructor`.
"""
function process_SciMLProblem(
        constructor, sys::AbstractSystem, u0map, pmap; build_initializeprob = true,
        implicit_dae = false, t = nothing, guesses = AnyDict(),
        warn_initialize_determined = true, initialization_eqs = [],
        eval_expression = false, eval_module = @__MODULE__, fully_determined = false,
        check_units = true, tofloat = true, use_union = false,
        u0_constructor = identity, du0map = nothing, check_length = true, symbolic_u0 = false, kwargs...)
    dvs = unknowns(sys)
    ps = parameters(sys)
    iv = has_iv(sys) ? get_iv(sys) : nothing
    eqs = equations(sys)

    check_array_equations_unknowns(eqs, dvs)

    u0Type = typeof(u0map)
    pType = typeof(pmap)
    _u0map = u0map
    u0map = to_varmap(u0map, dvs)
    _pmap = pmap
    pmap = to_varmap(pmap, ps)
    defs = add_toterms(defaults(sys))
    cmap, cs = get_cmap(sys)
    kwargs = NamedTuple(kwargs)

    op = add_toterms(u0map)
    missing_unknowns = add_fallbacks!(op, dvs, defs)
    for (k, v) in defs
        haskey(op, k) && continue
        op[k] = v
    end
    merge!(op, pmap)
    missing_pars = add_fallbacks!(op, ps, defs)
    for eq in cmap
        op[eq.lhs] = eq.rhs
    end
    if sys isa ODESystem
        guesses = merge(ModelingToolkit.guesses(sys), todict(guesses))
        has_observed_u0s = any(
            k -> has_observed_with_lhs(sys, k) || has_parameter_dependency_with_lhs(sys, k),
            keys(op))
        solvablepars = [p
                        for p in parameters(sys)
                        if is_parameter_solvable(p, pmap, defs, guesses)]
        if build_initializeprob &&
           (((implicit_dae || has_observed_u0s || !isempty(missing_unknowns) ||
              !isempty(solvablepars)) &&
             get_tearing_state(sys) !== nothing) ||
            !isempty(initialization_equations(sys))) && t !== nothing
            initializeprob = ModelingToolkit.InitializationProblem(
                sys, t, u0map, pmap; guesses, warn_initialize_determined,
                initialization_eqs, eval_expression, eval_module, fully_determined, check_units)
            initializeprobmap = getu(initializeprob, unknowns(sys))

            punknowns = [p
                         for p in all_variable_symbols(initializeprob)
                         if is_parameter(sys, p)]
            getpunknowns = getu(initializeprob, punknowns)
            setpunknowns = setp(sys, punknowns)
            initializeprobpmap = GetUpdatedMTKParameters(getpunknowns, setpunknowns)

            reqd_syms = parameter_symbols(initializeprob)
            update_initializeprob! = UpdateInitializeprob(
                getu(sys, reqd_syms), setu(initializeprob, reqd_syms))
            for p in punknowns
                p = unwrap(p)
                stype = symtype(p)
                op[p] = get_temporary_value(p)
                delete!(missing_pars, p)
            end

            for v in missing_unknowns
                op[v] = zero_var(v)
            end
            empty!(missing_unknowns)
            kwargs = merge(kwargs,
                (; initializeprob, initializeprobmap,
                    initializeprobpmap, update_initializeprob!))
        end
    end

    if t !== nothing && !(constructor <: Union{DDEFunction, SDDEFunction})
        op[iv] = t
    end

    add_observed!(sys, op)
    add_parameter_dependencies!(sys, op)

    evaluate_varmap!(op)

    u0 = better_varmap_to_vars(
        op, dvs; tofloat = true, use_union = false,
        container_type = u0Type, allow_symbolic = symbolic_u0)

    if u0 !== nothing
        u0 = u0_constructor(u0)
    end

    check_eqs_u0(eqs, dvs, u0; check_length, kwargs...)

    if is_split(sys)
        p = MTKParameters(sys, op)
    else
        p = better_varmap_to_vars(op, ps; tofloat, use_union, container_type = pType)
    end

    if implicit_dae && du0map !== nothing
        ddvs = map(Differential(iv), dvs)
        du0map = to_varmap(du0map, ddvs)
        merge!(op, du0map)

        du0 = varmap_to_vars(du0map, ddvs; toterm = identity,
            tofloat = true)
        kwargs = merge(kwargs, (; ddvs))
    else
        du0 = nothing
    end

    f = constructor(sys, dvs, ps, u0; p = p,
        eval_expression = eval_expression,
        eval_module = eval_module,
        kwargs...)
    implicit_dae ? (f, du0, u0, p) : (f, u0, p)
end
