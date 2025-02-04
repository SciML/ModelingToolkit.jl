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
    if !allow_symbolic
        for (sym, val) in zip(vars, vals)
            symbolic_type(val) == NotSymbolic() && continue
            throw(UnexpectedSymbolicValueInVarmap(sym, val))
        end
    end

    if container_type <: Union{AbstractDict, Tuple, Nothing, SciMLBase.NullParameters}
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
"""
function filter_missing_values!(varmap::AbstractDict)
    filter!(kvp -> kvp[2] !== nothing, varmap)
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
"""
function build_operating_point(
        u0map::AbstractDict, pmap::AbstractDict, defs::AbstractDict, cmap, dvs, ps)
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
    return op, missing_unknowns, missing_pars
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
        guesses, missing_unknowns; implicit_dae = false, kwargs...)
    guesses = merge(ModelingToolkit.guesses(sys), todict(guesses))
    has_observed_u0s = any(
        k -> has_observed_with_lhs(sys, k) || has_parameter_dependency_with_lhs(sys, k),
        keys(op))
    solvablepars = [p
                    for p in parameters(sys)
                    if is_parameter_solvable(p, pmap, defs, guesses)]
    has_dependent_unknowns = any(unknowns(sys)) do sym
        val = get(op, sym, nothing)
        val === nothing && return false
        return symbolic_type(val) != NotSymbolic() || is_array_of_symbolics(val)
    end
    if (((implicit_dae || has_observed_u0s || !isempty(missing_unknowns) ||
          !isempty(solvablepars) || has_dependent_unknowns) &&
         (!has_tearing_state(sys) || get_tearing_state(sys) !== nothing)) ||
        !isempty(initialization_equations(sys))) &&
       (!is_time_dependent(sys) || t !== nothing)
        initializeprob = ModelingToolkit.InitializationProblem(
            sys, t, u0map, pmap; guesses, kwargs...)

        if is_time_dependent(sys)
            all_init_syms = Set(all_symbols(initializeprob))
            solved_unknowns = filter(var -> var in all_init_syms, unknowns(sys))
            initializeprobmap = getu(initializeprob, solved_unknowns)
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
            p = unwrap(p)
            stype = symtype(p)
            op[p] = get_temporary_value(p)
            if iscall(p) && operation(p) === getindex
                arrp = arguments(p)[1]
                op[arrp] = collect(arrp)
            end
        end

        if is_time_dependent(sys)
            for v in missing_unknowns
                op[v] = zero_var(v)
            end
            empty!(missing_unknowns)
        end
        return (;
            initialization_data = SciMLBase.OverrideInitData(
                initializeprob, update_initializeprob!, initializeprobmap,
                initializeprobpmap))
    end
    return (;)
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
- `tofloat`, `use_union`: Passed to [`better_varmap_to_vars`](@ref) for building `u0` (and
  possibly `p`).
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

All other keyword arguments are passed as-is to `constructor`.
"""
function process_SciMLProblem(
        constructor, sys::AbstractSystem, u0map, pmap; build_initializeprob = true,
        implicit_dae = false, t = nothing, guesses = AnyDict(),
        warn_initialize_determined = true, initialization_eqs = [],
        eval_expression = false, eval_module = @__MODULE__, fully_determined = nothing,
        check_initialization_units = false, tofloat = true, use_union = false,
        u0_constructor = identity, du0map = nothing, check_length = true,
        symbolic_u0 = false, warn_cyclic_dependency = false,
        circular_dependency_max_cycle_length = length(all_symbols(sys)),
        circular_dependency_max_cycles = 10,
        substitution_limit = 100, use_scc = true, kwargs...)
    dvs = unknowns(sys)
    ps = parameters(sys)
    iv = has_iv(sys) ? get_iv(sys) : nothing
    eqs = equations(sys)

    check_array_equations_unknowns(eqs, dvs)

    u0Type = typeof(u0map)
    pType = typeof(pmap)
    _u0map = u0map
    u0map = to_varmap(u0map, dvs)
    symbols_to_symbolics!(sys, u0map)
    _pmap = pmap
    pmap = to_varmap(pmap, ps)
    symbols_to_symbolics!(sys, pmap)
    defs = add_toterms(recursive_unwrap(defaults(sys)))
    cmap, cs = get_cmap(sys)
    kwargs = NamedTuple(kwargs)

    op, missing_unknowns, missing_pars = build_operating_point(
        u0map, pmap, defs, cmap, dvs, ps)

    if build_initializeprob
        kws = maybe_build_initialization_problem(
            sys, op, u0map, pmap, t, defs, guesses, missing_unknowns;
            implicit_dae, warn_initialize_determined, initialization_eqs,
            eval_expression, eval_module, fully_determined,
            warn_cyclic_dependency, check_units = check_initialization_units,
            circular_dependency_max_cycle_length, circular_dependency_max_cycles, use_scc)

        kwargs = merge(kwargs, kws)
    end

    if t !== nothing && !(constructor <: Union{DDEFunction, SDDEFunction})
        op[iv] = t
    end

    add_observed!(sys, op)
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
        op, dvs; tofloat = true, use_union = false,
        container_type = u0Type, allow_symbolic = symbolic_u0)

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
        p = MTKParameters(sys, op)
    else
        p = better_varmap_to_vars(op, ps; tofloat, use_union, container_type = pType)
    end

    if implicit_dae && du0map !== nothing
        ddvs = map(Differential(iv), dvs)
        du0map = to_varmap(du0map, ddvs)
        merge!(op, du0map)
        du0 = varmap_to_vars(op, ddvs; toterm = identity,
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
        use_union = true,
        tofloat = true,
        symbolic_u0 = false)
    dvs = unknowns(sys)
    ps = parameters(sys)

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
        allobs = Set(getproperty.(observed(sys), :lhs))
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
        u0 = varmap_to_vars(u0map, dvs; defaults = defs, tofloat = true, use_union)
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
