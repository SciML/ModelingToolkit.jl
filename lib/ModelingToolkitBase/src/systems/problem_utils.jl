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
anydict(x::AbstractDict) = AnyDict(x)
function anydict(x)
    op = AnyDict()
    for (k, v) in x
        if haskey(op, k)
            throw(
                ArgumentError(
                    """
                    Found duplicate entries in symbolic map. Key $k is provided multiple \
                    times.
                    """
                )
            )
        end
        op[k] = v
    end
    return op
end

"""
    $(TYPEDSIGNATURES)

Check if `x` is a symbolic with known size. Assumes `SymbolicUtils.shape(unwrap(x))`
is a valid operation.
"""
symbolic_has_known_size(x) = !(SU.shape(unwrap(x)) isa SU.Unknown)

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
    return if is_split(sys)
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

Return the list of variables in `varlist` not present in `varmap`. Uses the same criteria
for missing array variables and `toterm` forms as [`add_fallbacks!`](@ref).
"""
function missingvars(
        varmap::AtomicArrayDict, varlist::Vector; toterm = default_toterm
    )
    missings = Set{SymbolicT}()
    for var in varlist
        var = unwrap(var)
        get_possibly_indexed(varmap, var, COMMON_NOTHING) === COMMON_NOTHING || continue
        ttsym = toterm(var)
        get_possibly_indexed(varmap, ttsym, COMMON_NOTHING) === COMMON_NOTHING || continue
        push!(missings, var)
    end
    return missings
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
    return symbolic_type(x) == ArraySymbolic() ? value(x) : recursive_unwrap.(x)
end

function recursive_unwrap(x::SparseMatrixCSC)
    I, J, V = findnz(x)
    V = recursive_unwrap(V)
    m, n = size(x)
    return sparse(I, J, V, m, n)
end

recursive_unwrap(x) = value(x)

function recursive_unwrap(x::AbstractDict)
    return anydict(unwrap(k) => recursive_unwrap(v) for (k, v) in x)
end

"""
    $(TYPEDSIGNATURES)

Add equations `eqs` to `varmap`. Assumes each element in `eqs` maps a single symbolic
variable to an expression representing its value. In case `varmap` already contains an
entry for `eq.lhs`, insert the reverse mapping if `eq.rhs` is not a number.
"""
function add_observed_equations!(varmap::AtomicArrayDict{SymbolicT}, eqs::Vector{Equation}, bound_ps::Union{Nothing, ROSymmapT} = nothing)
    for eq in eqs
        if get_possibly_indexed(varmap, eq.lhs, COMMON_NOTHING) !== COMMON_NOTHING
            SU.isconst(eq.rhs) && continue
            get_possibly_indexed(varmap, eq.rhs, COMMON_NOTHING) !== COMMON_NOTHING && continue
            bound_ps isa ROSymmapT && has_possibly_indexed_key(parent(bound_ps), eq.rhs) && continue
            Moshi.Match.@match eq.rhs begin
                BSImpl.Term(; f, args) && if f isa SymbolicT end => nothing
                BSImpl.Sym() => nothing
                _ => continue
            end
            write_possibly_indexed_array!(varmap, eq.rhs, eq.lhs, COMMON_NOTHING)
        else
            write_possibly_indexed_array!(varmap, eq.lhs, eq.rhs, COMMON_NOTHING)
        end
    end
    return
end

"""
    $(TYPEDSIGNATURES)

Add all equations in `observed(sys)` to `varmap` using [`add_observed_equations!`](@ref).
"""
function add_observed!(sys::AbstractSystem, varmap::AbstractDict)
    return add_observed_equations!(varmap, observed(sys))
end

struct UnexpectedSymbolicValueInVarmap <: Exception
    sym::Any
    val::Any
end

function Base.showerror(io::IO, err::UnexpectedSymbolicValueInVarmap)
    return println(
        io,
        """
        Found symbolic value $(err.val) for variable $(err.sym). You may be missing an \
        initial condition or have cyclic initial conditions. If this is intended, pass \
        `symbolic_u0 = true`. In case the initial conditions are not cyclic but \
        require more substitutions to resolve, increase `substitution_limit`. To report \
        cycles in initial conditions of unknowns/parameters, pass \
        `warn_cyclic_dependency = true`. If the cycles are still not reported, you \
        may need to pass a larger value for `circular_dependency_max_cycle_length` \
        or `circular_dependency_max_cycles`.
        """
    )
end

struct MissingGuessError <: Exception
    syms::Vector{Any}
    vals::Vector{Any}
end

function Base.showerror(io::IO, err::MissingGuessError)
    println(
        io,
        """
        Cyclic guesses detected in the system. Symbolic values were found for the following \
        variables/parameters in the map: \
        """
    )
    for (sym, val) in zip(err.syms, err.vals)
        println(io, "$sym  => $val")
    end
    return println(
        io,
        """
        In order to resolve this, please provide additional numeric guesses so that the \
        chain can be resolved to assign numeric values to each variable. Alternatively, the \
        `missing_guess_value` keyword can be used to set a fallback guess for all \
        variables. The keyword must be passed an instance of the `MissingGuessValue` sum-type.
        """
    )
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
    return println(io, join(e.vars, ", "))
end

"""
    $TYPEDEF

A Moshi.jl enum to allow choosing what happens with missing guess values when building a
numerical problem from a `System`.

# Variants

- `MissingGuessValue.Constant(val::Number)`: Missing guesses are set to the given value
  `val`.
- `MissingGuessValue.Random(rng::AbstractRNG)`: Missing guesses are set to `rand(rng)`.
- `MissingGuessValue.HashedRandom`: Missing guesses are set to a
  deterministically determined random-like value based on the hash of the variable name
- `MissingGuessValue.Error()`: Missing guess values cause an error.
"""
Moshi.Data.@data MissingGuessValue begin
    Constant(Number)
    Random(AbstractRNG)
    HashedRandom
    Error
end

# To be overloaded downstream by MTK
default_missing_guess_value() = default_missing_guess_value(nothing)
default_missing_guess_value(_) = MissingGuessValue.HashedRandom()

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
function varmap_to_vars(
        varmap::AbstractDict, vars::Vector; ir = nothing,
        tofloat = true, use_union = false, container_type = Array, buffer_eltype = Nothing,
        toterm = default_toterm, check = true, allow_symbolic = false,
        is_initializeprob = false, substitution_limit = 100, missing_values = MissingGuessValue.Error()
    )
    isempty(vars) && return nothing

    if !(varmap isa SymmapT)
        varmap = as_atomic_dict_with_defaults(Dict{SymbolicT, SymbolicT}(varmap), COMMON_NOTHING)
    end
    if toterm !== nothing
        add_toterms!(varmap; toterm)
    end
    if ir !== nothing
        evaluate_varmap!(
            ir, AtomicArrayDictSubstitutionWrapper(varmap), vars;
            limit = substitution_limit,
            allow_symbolic = allow_symbolic ||
                !Moshi.Data.isa_variant(missing_values, MissingGuessValue.Error)
        )
    else
        evaluate_varmap!(
            AtomicArrayDictSubstitutionWrapper(varmap), vars;
            limit = substitution_limit,
            allow_symbolic = allow_symbolic ||
                !Moshi.Data.isa_variant(missing_values, MissingGuessValue.Error)
        )
    end
    if check && !allow_symbolic
        missing_vars = missingvars(varmap, vars; toterm)
        for var in vars
            var = unwrap(var)
            val = get_possibly_indexed(varmap, var, COMMON_NOTHING)
            SU.isconst(val) || push!(missing_vars, var)
        end
        Moshi.Match.@match missing_values begin
            MissingGuessValue.Constant(val) => begin
                cval = BSImpl.Const{VartypeT}(val)
                for var in missing_vars
                    if Symbolics.isarraysymbolic(var)
                        varmap[var] = BSImpl.Const{VartypeT}(fill(val, size(var)))
                    else
                        write_possibly_indexed_array!(varmap, var, cval, COMMON_NOTHING)
                    end
                end
            end
            MissingGuessValue.Random(rng) => begin
                for var in missing_vars
                    if Symbolics.isarraysymbolic(var)
                        varmap[var] = rand(rng, size(var))
                    else
                        write_possibly_indexed_array!(varmap, var, Symbolics.SConst(rand(rng)), COMMON_NOTHING)
                    end
                end
            end
            MissingGuessValue.HashedRandom() => begin
                for var in missing_vars
                    if Symbolics.isarraysymbolic(var)
                        varmap[var] = [hash(var, hash(i)) for i in SU.stable_eachindex(var)] ./ 0x1p64
                    else
                        write_possibly_indexed_array!(varmap, var, Symbolics.SConst(hash(var) / 0x1p64), COMMON_NOTHING)
                    end
                end
            end
            MissingGuessValue.Error() => begin
                if !isempty(missing_vars)
                    if is_initializeprob
                        throw(MissingGuessError(collect(missing_vars), collect(missing_vars)))
                    else
                        throw(MissingVariablesError(missing_vars))
                    end
                end

            end
        end
    end
    vals = map(vars) do x
        x = unwrap(x)
        v = get_possibly_indexed(varmap, x, x)
        Moshi.Match.@match v begin
            BSImpl.Const(; val) => return val
            _ => begin
                Moshi.Match.@match x begin
                    BSImpl.Term(; f, args) && if f isa Initial && isequal(v, args[1]) end => begin
                        if Symbolics.isarraysymbolic(v)
                            return fill(false, size(v))
                        else
                            return false
                        end
                    end
                    _ => return v
                end
            end
        end
    end
    if !allow_symbolic
        missingsyms = Any[]
        missingvals = Any[]
        for (sym, val) in zip(vars, vals)
            val !== nothing && symbolic_type(val) == NotSymbolic() && continue
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
        return SymbolicUtils.Code.create_array(
            container_type, eltype(vals), Val{1}(),
            Val(length(vals)), vals...
        )
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
        varmap::AbstractDict, vars; max_cycle_length = length(varmap), max_cycles = 10
    )
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
            SU.search_variables!(buffer, v)
            for var in buffer
                haskey(var_to_idx, var) || continue
                add_edge!(graph, kidx, var_to_idx[var])
            end
        elseif v isa AbstractArray
            for val in v
                SU.search_variables!(buffer, val)
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

    return map(cycles) do cycle
        map(Base.Fix1(getindex, allvars), cycle)
    end
end

"""
    $(TYPEDSIGNATURES)

Performs symbolic substitution on the values in `varmap` for the keys in `vars`, using
`varmap` itself as the set of substitution rules. If an entry in `vars` is not a key
in `varmap`, it is ignored.
"""
function evaluate_varmap!(varmap::AbstractDict{SymbolicT, SymbolicT}, vars; limit = 100, allow_symbolic = false)
    for k in vars
        arr, _ = split_indexed_var(unwrap(k))
        v = get(varmap, arr, COMMON_NOTHING)
        v === COMMON_NOTHING && continue
        SU.isconst(v) && continue
        varmap[arr] = fixpoint_sub(v, varmap; maxiters = limit, fold = Val(true), warn_maxiters = !allow_symbolic)
    end
    return
end

function evaluate_varmap!(
        ir::IRStructure{SymReal}, varmap::AtomicArrayDictSubstitutionWrapper, vars;
        limit = 100, allow_symbolic = false
    )
    subber = Symbolics.FixpointSubstituter(
        SU.IRSubstituter{true}(ir, varmap; filterer = Symbolics.FPSubFilterer{Nothing}());
        maxiters = limit, warn_maxiters = !allow_symbolic
    )
    for k in vars
        v = get(varmap, k, COMMON_NOTHING)
        v === COMMON_NOTHING && continue
        SU.isconst(v) && continue
        varmap[k] = subber(v)
    end
    return
end

"""
    $(TYPEDSIGNATURES)

Remove keys in `varmap` whose values are `nothing`.

If `missing_values` is not `nothing`, it is assumed to be a collection and all removed
keys will be added to it.
"""
function filter_missing_values!(varmap::AbstractDict; missing_values = nothing)
    return filter!(varmap) do kvp
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
        symbolic_has_known_size(var) || continue
        haskey(varmap, var) || continue
        for i in eachindex(var)
            haskey(varmap, var[i]) && continue
            varmap[var[i]] = varmap[var][i]
        end
    end
    return
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
        error(lazy"Nonnumeric parameter $p with symtype $stype cannot be solved for during initialization")
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
    $TYPEDSIGNATURES

For every `Initial(x)` parameter in `sys`, add `Initial(x) => x` to `op` if it does not
already contain that key.
"""
function add_initials!(sys::AbstractSystem, op::SymmapT)
    for p in get_ps(sys)
        haskey(op, p) && continue
        Moshi.Match.@match p begin
            BSImpl.Term(; f, args) && if f isa Initial end => begin
                write_possibly_indexed_array!(
                    op, p, if Symbolics.isarraysymbolic(p)
                        BSImpl.Const{VartypeT}(fill(false, size(p)))
                    else
                        COMMON_FALSE
                    end, COMMON_FALSE
                )
            end
            _ => nothing
        end
    end
    return
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
    return ow.f(state_values(prob), parameter_values(prob), current_time(prob))
end

function (ow::ObservedWrapper{false})(prob)
    return ow.f(state_values(prob), parameter_values(prob))
end

"""
    $(TYPEDSIGNATURES)

Given an index provider `indp` and a vector of symbols `syms` return a type-stable getter
function.

Note that the getter ONLY works for problem-like objects, since it generates an observed
function. It does NOT work for solutions.
"""
Base.@nospecializeinfer function concrete_getu(
        indp, syms; wrap_as_any = false,
        eval_expression, eval_module, force_time_independent = false, kwargs...
    )
    @nospecialize
    obsfn = build_explicit_observed_function(
        indp, syms; wrap_delays = false, eval_expression, eval_module,
        force_time_independent, kwargs...
    )
    if wrap_as_any
        return ObservedWrapper{is_time_dependent(indp) && !force_time_independent, Any}(obsfn)
    end
    return ObservedWrapper{is_time_dependent(indp) && !force_time_independent}(obsfn)
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
    return pca.p_constructor(x)
end

function (pca::PConstructorApplicator)(x::AbstractArray{Bool})
    return pca.p_constructor(BitArray(x))
end

function (pca::PConstructorApplicator{typeof(identity)})(x::SubArray)
    return collect(x)
end

function (pca::PConstructorApplicator{typeof(identity)})(x::SubArray{Bool})
    return BitArray(x)
end

function (pca::PConstructorApplicator{typeof(identity)})(x::SubArray{<:AbstractArray})
    return collect(pca.(x))
end

function (pca::PConstructorApplicator)(x::AbstractArray{<:AbstractArray})
    return pca.p_constructor(pca.(x))
end

"""
    $TYPEDEF

Callable struct designed for use by `MTKParametersReconstructor`. Uses a fixed set of templates to
act as a very dynamic (and limited) observed function returning an array. See `__apply_copy_template`
for the supported templates.
"""
struct CopyParamsByTemplate{IsRoot, T, N, G}
    """
    List of templates.
    """
    template::T # TODO: This field is parametric because I thought we might want to specialize it in some cases.
    """
    Size of the returned buffer.
    """
    size::NTuple{N, Int}
    """
    Merged getter for all symbolic-fallback batches, or `nothing`.
    """
    fallback_getter::G
end

function CopyParamsByTemplate{IR}(
        temp::T, size::NTuple{N, Int}, fallback_getter::G = nothing
    ) where {IR, T, N, G}
    return CopyParamsByTemplate{IR, T, N, G}(temp, size, fallback_getter)
end

"""
    $TYPEDEF

Template entry indexing into the result of a `CopyParamsByTemplate`'s merged
`fallback_getter`. `range` is the contiguous slice of that result corresponding to this
fallback batch's symbols.
"""
struct FallbackSlice
    range::UnitRange{Int}
end

function __apply_copy_template(valp, template)
    p = parameter_values(valp)
    u = state_values(valp)
    if template isa ParameterIndex{SciMLStructures.Tunable, UnitRange{Int}}
        if p isa MTKParameters
            return p.tunable[template.idx]
        else
            return p[template.idx]
        end
    elseif template isa ParameterIndex{SciMLStructures.Initials, UnitRange{Int}}
        return p.initials[template.idx]
    elseif template isa ParameterIndex{SciMLStructures.Discrete, Tuple{Int, UnitRange{Int}}}
        return p.discrete[template.idx[1]][template.idx[2]]
    elseif template isa ParameterIndex{SciMLStructures.Constants, Tuple{Int, UnitRange{Int}}}
        return p.constant[template.idx[1]][template.idx[2]]
    elseif template isa ParameterIndex{Nonnumeric, Tuple{Int, UnitRange{Int}}}
        return p.nonnumeric[template.idx[1]][template.idx[2]]
    elseif template isa StaticBufferIndex{SciMLStructures.Discrete}
        return _static_buffer(p.discrete, template)[template.range]
    elseif template isa StaticBufferIndex{SciMLStructures.Constants}
        return _static_buffer(p.constant, template)[template.range]
    elseif template isa StaticBufferIndex{Nonnumeric}
        return _static_buffer(p.nonnumeric, template)[template.range]
    elseif template isa UnitRange{Int}
        return u[template]
    elseif template isa ObservedWrapper
        return template(valp)
    elseif template isa CopyParamsByTemplate
        return template(valp)
    elseif template isa IndepVarTemplate
        return current_time(valp)
    elseif template isa ParameterIndex{SciMLStructures.Constants, <:Tuple{Vararg{Int}}}
        i, j, rest... = template.idx
        return p.constant[i][j][rest...]
    elseif template isa ParameterIndex{SciMLStructures.Nonnumeric, <:Tuple{Vararg{Int}}}
        i, j, rest... = template.idx
        return p.nonnumeric[i][j][rest...]
    else
        # MethodError because this is a manual dispatch chain
        throw(MethodError(__apply_copy_template, (valp, template)))
    end
end

@inline function __apply_root_template(src, template, fb)
    if template isa FallbackSlice
        return fb[template.range]
    else
        return __apply_copy_template(src, template)
    end
end

function (cp::CopyParamsByTemplate{IsRoot})(src) where {IsRoot}
    return if IsRoot
        if cp.fallback_getter === nothing
            reshape(mapreduce(Base.Fix1(__apply_copy_template, src), vcat, cp.template), cp.size)
        else
            fb = cp.fallback_getter(src)
            reshape(
                mapreduce(t -> __apply_root_template(src, t, fb), vcat, cp.template), cp.size
            )
        end
    else
        buffers = map(Base.Fix1(__apply_copy_template, src), cp.template)
        if cp.template isa Tuple
            buffers = collect(buffers)
        end
        reshape(buffers, cp.size)
    end
end

struct IndepVarTemplate end
const IV_TEMPLATE = IndepVarTemplate()

"""
    $TYPEDEF

Template entry for `CopyParamsByTemplate` indexing into one of the inner buffers of a
multi-buffer `MTKParameters` portion (discrete/constants/nonnumeric). Unlike
`ParameterIndex{P, Tuple{Int, UnitRange{Int}}}`, the buffer index `I` is lifted into the
type domain so that indexing the heterogeneously-typed tuple of buffers constant-folds and
infers concretely. With a runtime buffer index the result is a small `Union` of the buffer
types, which Enzyme's type analysis rejects (`IllegalTypeAnalysisException`) when it flows
into `reshape` inside the `CopyParamsByTemplate` compile unit.
"""
struct StaticBufferIndex{P, I}
    range::UnitRange{Int}
end

function StaticBufferIndex{P}(idx::Tuple{Int, UnitRange{Int}}) where {P}
    return StaticBufferIndex{P, idx[1]}(idx[2])
end

@inline _static_buffer(bufs::Tuple, ::StaticBufferIndex{P, I}) where {P, I} = bufs[I]

Base.@nospecializeinfer function __specialize_templates(template::Vector{Any}, elem_types::Set{DataType})
    if length(template) <= 4
        return Tuple(template)
    elseif length(elem_types) <= 4
        return Vector{Union{collect(elem_types)...}}(template)
    else
        return template
    end
end

# Memo for the symbolic-fallback getters built by `CopyParamsByTemplate`. The init problem
# builds several `CopyParamsByTemplate`s over the same `initsys` (e.g. `GetUpdatedU0` and the
# `initializeprobmap`)
abstract type TemplateGetuCache end
const TemplateGetuCacheT = Dict{Vector{SymbolicT}, Any}

function should_invalidate_mutable_cache_entry(::Type{TemplateGetuCache}, @nospecialize(patch::NamedTuple))
    return false
end

function cached_template_getu(srcsys::AbstractSystem, batch::Vector{SymbolicT}; kws...)
    if !(srcsys isa System)
        return concrete_getu(srcsys, Symbolics.SConst(batch); wrap_as_any = true, kws...)
    end
    cache = check_mutable_cache(srcsys, TemplateGetuCache, TemplateGetuCacheT, nothing)
    if cache === nothing
        cache = TemplateGetuCacheT()
        store_to_mutable_cache!(srcsys, TemplateGetuCache, cache)
    end
    return get!(cache, batch) do
        concrete_getu(srcsys, Symbolics.SConst(batch); wrap_as_any = true, kws...)
    end
end

function CopyParamsByTemplate(srcsys::AbstractSystem, syms::AbstractArray{SymbolicT}; kws...)
    template = []
    elem_types = Set{DataType}()
    iv = get_iv(srcsys)
    irinfo = get_ir_info(srcsys)

    for sym in syms
        if iv isa SymbolicT && isequal(iv, sym)
            push!(template, IV_TEMPLATE)
            push!(elem_types, IndepVarTemplate)
            continue
        end
        symidx = parameter_index(srcsys, sym)
        if symidx === nothing
            symidx = parameter_index(srcsys, irinfo.obs_subber(sym))
        end
        if symidx === nothing
            symidx = variable_index(srcsys, sym)
            if symidx === nothing
                symidx = variable_index(srcsys, irinfo.obs_subber(sym))
            end
            if symidx === nothing
                if isempty(template)
                    push!(elem_types, Vector{SymbolicT})
                    push!(template, SymbolicT[sym])
                    continue
                end
                prev = template[end]
                if prev isa Vector{SymbolicT}
                    push!(prev, sym)
                else
                    push!(elem_types, Vector{SymbolicT})
                    push!(template, SymbolicT[sym])
                end
                continue
            end
            if isempty(template)
                push!(elem_types, UnitRange{Int})
                push!(template, symidx:symidx)
                continue
            end
            prev = template[end]
            if prev isa UnitRange{Int} && last(prev) + 1 == symidx
                template[end] = first(prev):symidx
            else
                push!(elem_types, UnitRange{Int})
                push!(template, symidx:symidx)
            end
            continue
        elseif symidx isa Int
            symidx = ParameterIndex(SciMLStructures.Tunable(), symidx)
        end
        portion = symidx.portion
        _bufidx = symidx.idx
        bufidx::UnitRange{Int} = if _bufidx isa AbstractVector{Int}
            @assert isequal(vec(_bufidx), first(_bufidx):last(_bufidx))
            subidx = nothing
            first(_bufidx):last(_bufidx)
        elseif _bufidx isa Int
            subidx = nothing
            _bufidx:_bufidx
        elseif _bufidx isa NTuple{2, Int}
            subidx = _bufidx[1]
            _bufidx[2]:_bufidx[2]
        elseif _bufidx isa Tuple{Vararg{Int}} # indexing into a non-tunable array parameter
            push!(template, symidx)
            continue
        else
            # Will error due to the typeassert on `bufidx`
            nothing
        end
        if isempty(template)
            pidx = if subidx === nothing
                ParameterIndex(symidx.portion, bufidx)
            else
                ParameterIndex(symidx.portion, (subidx, bufidx))
            end
            push!(template, pidx)
            push!(elem_types, typeof(pidx))
            continue
        end
        prev = template[end]
        if prev isa ParameterIndex && prev.portion === symidx.portion && (
                subidx === nothing && prev.idx isa UnitRange{Int} &&
                    last(prev.idx) + 1 == first(bufidx) ||
                    prev.idx isa Tuple{Int, UnitRange{Int}} && subidx == prev.idx[1] &&
                    last(prev.idx[2]) + 1 == first(bufidx)
            )
            if subidx === nothing
                template[end] = ParameterIndex(prev.portion, first(prev.idx):last(bufidx))
            else
                template[end] = ParameterIndex(prev.portion, (subidx, first(prev.idx[2]):last(bufidx)))
            end
        elseif subidx === nothing
            push!(template, ParameterIndex(symidx.portion, bufidx))
            push!(elem_types, typeof(template[end]))
        else
            push!(template, ParameterIndex(symidx.portion, (subidx, bufidx)))
            push!(elem_types, typeof(template[end]))
        end
    end

    fallback_idxs = Int[]
    for i in eachindex(template)
        template[i] isa Vector{SymbolicT} && push!(fallback_idxs, i)
    end
    fallback_getter = nothing
    if length(fallback_idxs) == 1
        i = fallback_idxs[1]
        template[i] = cached_template_getu(srcsys, template[i]::Vector{SymbolicT}; kws...)
        delete!(elem_types, Vector{SymbolicT})
        push!(elem_types, typeof(template[i]))
    elseif length(fallback_idxs) > 1
        all_fallback = SymbolicT[]
        for i in fallback_idxs
            append!(all_fallback, template[i]::Vector{SymbolicT})
        end
        fallback_getter = cached_template_getu(srcsys, all_fallback; kws...)
        offset = 0
        for i in fallback_idxs
            len = length(template[i]::Vector{SymbolicT})
            template[i] = FallbackSlice((offset + 1):(offset + len))
            offset += len
        end
        delete!(elem_types, Vector{SymbolicT})
        push!(elem_types, FallbackSlice)
    end

    # Lift buffer indices of multi-buffer portions (discrete/constants/nonnumeric) into
    # the type domain. This is done as a final pass so the contiguous-range merging above
    # can keep operating on plain `ParameterIndex`es.
    for i in eachindex(template)
        entry = template[i]
        # Only lift the `(bufidx, range)` form into the type domain.
        if entry isa ParameterIndex && entry.portion isa Union{SciMLStructures.Discrete, SciMLStructures.Constants, Nonnumeric} && entry.idx isa Tuple{Int, UnitRange{Int}}
            delete!(elem_types, typeof(entry))
            template[i] = StaticBufferIndex{typeof(entry.portion)}(entry.idx)
            push!(elem_types, typeof(template[i]))
        end
    end

    return CopyParamsByTemplate{true}(
        __specialize_templates(template, elem_types), size(syms), fallback_getter
    )
end

function CopyParamsByTemplate(srcsys::AbstractSystem, syms::AbstractArray; kws...)
    template = []
    elem_types = Set{DataType}()
    for sym in syms
        push!(template, CopyParamsByTemplate(srcsys, sym; kws...))
        push!(elem_types, typeof(template[end]))
    end
    return CopyParamsByTemplate{false}(__specialize_templates(template, elem_types), size(syms))
end

struct MTKParametersReconstructor{T, I, D, C, N}
    tunables_fn::T
    initials_fn::I
    discretes_fn::D
    consts_fn::C
    nonnumerics_fn::N
    diffcache_buffer_idx::Int
end

# TODO: make this infer when the nonnumerics are non-trivial
function (recon::MTKParametersReconstructor)(src, dst)
    src_ps = parameter_values(src)
    dst_ps = parameter_values(dst)
    oldcache = dst_ps.caches
    # I don't know why but this makes it infer properly
    if recon.tunables_fn isa ComposedFunction
        tunablevals = recon.tunables_fn.outer(recon.tunables_fn.inner(src))
    else
        tunablevals = recon.tunables_fn(src)
    end
    initialvals = recon.initials_fn(src)
    nonnumerics = recon.nonnumerics_fn(src)::typeof(dst_ps.nonnumeric)
    (; diffcache_buffer_idx) = recon
    if !iszero(diffcache_buffer_idx)
        @set! nonnumerics[diffcache_buffer_idx] = DiffCacheAllocatorAPIWrapper{ForwardDiff.valtype(eltype(initialvals))}.(nonnumerics[diffcache_buffer_idx])
    end
    return MTKParameters(
        tunablevals, initialvals, recon.discretes_fn(src),
        recon.consts_fn(src), nonnumerics, oldcache isa Tuple{} ? () : copy.(oldcache)
    )
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
function MTKParametersReconstructor(
        srcsys::AbstractSystem, dstsys::AbstractSystem;
        initials = false, unwrap_initials = false, p_constructor = identity,
        force_time_independent = false,
        kwargs...
    )
    _p_constructor = p_constructor
    p_constructor = PConstructorApplicator(p_constructor)
    # if we call `getu` on this (and it were able to handle empty tuples) we get the
    # fields of `MTKParameters` except caches.
    syms = reorder_parameters(
        dstsys, parameters(dstsys; initial_parameters = initials); flatten = false
    )
    # `dstsys` is an initialization system, do basically everything is a tunable
    # and tunables are a mix of different types in `srcsys`. No initials. Constants
    # are going to be constants in `srcsys`, as are `nonnumeric`.

    # `syms[1]` is always the tunables because `srcsys` will have initials.
    tunable_syms = syms[1]
    tunable_getter = if isempty(tunable_syms)
        Returns(SVector{0, Float64}())
    else
        p_constructor ∘ CopyParamsByTemplate(srcsys, tunable_syms; kwargs...)
    end
    initials_getter = if initials && !isempty(syms[2])
        initsyms = syms[2]::Vector{SymbolicT}
        allsyms = Set{SymbolicT}(variable_symbols(srcsys))
        if unwrap_initials
            for i in eachindex(initsyms)
                sym = initsyms[i]
                arr, isarr = split_indexed_var(sym)
                innersym = if isarr
                    sidx = get_stable_index(sym)
                    first(arguments(arr))[sidx]
                else
                    first(arguments(arr))
                end
                if innersym in allsyms
                    initsyms[i] = innersym
                end
            end
        end
        p_constructor ∘ CopyParamsByTemplate(srcsys, initsyms; kwargs...)
    else
        Returns(SVector{0, Float64}())
    end
    discs_getter = if isempty(syms[3])
        Returns(())
    else
        ic = get_index_cache(dstsys)
        blockarrsizes = Tuple(
            map(ic.discrete_buffer_sizes) do bufsizes
                p_constructor(map(x -> x.length, bufsizes))
            end
        )

        # discretes need to be blocked arrays
        # the `getu` returns a tuple of arrays corresponding to `p.discretes`
        # `Base.Fix1(...)` applies `p_constructor` to each of the arrays in the tuple
        # `Base.Fix2(...)` does `BlockedArray.(tuple_of_arrs, blockarrsizes)` returning a
        # tuple of `BlockedArray`s
        Base.Fix2(Broadcast.BroadcastFunction(BlockedArray), blockarrsizes) ∘
            Base.Fix1(broadcast, p_constructor) ∘ Tuple ∘
            # This `broadcast.(collect, ...)` avoids `ReshapedArray`/`SubArray`s from
            # appearing in the result.
            CopyParamsByTemplate(srcsys, broadcast.(collect, syms[3]); kwargs...)
    end
    const_getter = if isempty(syms[4])
        Returns(())
    else
        Base.Fix1(broadcast, p_constructor) ∘ Tuple ∘ CopyParamsByTemplate(srcsys, syms[4]; kwargs...)
    end
    diffcache_buffer_idx = 0
    nonnumeric_getter = if isempty(syms[5])
        Returns(())
    else
        ic = get_index_cache(dstsys)
        buftypes = Tuple(
            map(ic.nonnumeric_buffer_sizes) do bufsize
                Vector{bufsize.type}
            end
        )

        diffcache_params = SU.getmetadata(dstsys, DiffCacheParams, Dict{SymbolicT, Int}())::Dict{SymbolicT, Int}
        if !isempty(diffcache_params)
            representative = first(keys(diffcache_params))
            diffcache_buffer_idx, _ = ic.nonnumeric_idx[representative]
            @set! buftypes[diffcache_buffer_idx] = identity
            for (i, sym) in enumerate(syms[5][diffcache_buffer_idx])
            end
        end
        # nonnumerics retain the assigned buffer type without narrowing
        Base.Fix1(broadcast, _p_constructor) ∘
            Base.Fix1(Broadcast.BroadcastFunction(call), buftypes) ∘ Tuple ∘ CopyParamsByTemplate(srcsys, syms[5]; kwargs...)
    end

    return MTKParametersReconstructor(tunable_getter, initials_getter, discs_getter, const_getter, nonnumeric_getter, diffcache_buffer_idx)
end

function call(f, args...)
    return f(args...)
end

"""
    $(TYPEDSIGNATURES)

Construct a `ReconstructInitializeprob` which reconstructs the `u0` and `p` of `dstsys`
with values from `srcsys`.

Extra keyword arguments are forwarded to `build_function_wrapper`.
"""
function ReconstructInitializeprob(
        srcsys::AbstractSystem, dstsys::AbstractSystem; u0_constructor = identity, p_constructor = identity,
        eval_expression = false, eval_module = @__MODULE__, is_steadystateprob = false, kwargs...
    )
    @assert is_initializesystem(dstsys)
    ugetter = u0_constructor ∘
        concrete_getu(
        srcsys, unknowns(dstsys);
        eval_expression, eval_module, force_time_independent = is_steadystateprob,
        iip_config = (true, false),
        kwargs...
    )
    if is_split(dstsys)
        pgetter = MTKParametersReconstructor(
            srcsys, dstsys; p_constructor, eval_expression, eval_module,
            force_time_independent = is_steadystateprob, kwargs...
        )
    else
        syms = parameters(dstsys)
        pgetter = let inner = concrete_getu(
                srcsys, syms; eval_expression, eval_module,
                force_time_independent = is_steadystateprob, kwargs...
            ),
                p_constructor = p_constructor

            function _getter2(valp, initprob)
                return p_constructor(inner(valp))
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
    if T != eltype(u0) && T != Union{} && T !== Any
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
        if eltype(buf) != T && !(buf isa SVector{0})
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
        sys::AbstractSystem, initsys::AbstractSystem; p_constructor = identity, eval_expression, eval_module,
        kwargs...
    )
    @assert is_initializesystem(initsys)
    if is_split(sys)
        return let getter = MTKParametersReconstructor(
                initsys, sys; initials = true, unwrap_initials = true, p_constructor,
                eval_expression, eval_module, kwargs...
            )
            function initprobpmap_split(prob, initsol)
                return getter(initsol, prob)
            end
        end
    else
        return let getter = concrete_getu(
                initsys, parameters(sys; initial_parameters = true);
                eval_expression, eval_module, kwargs...
            ), p_constructor = p_constructor

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

Any changes to this method should also be made to the one in ChainRulesCoreExt.
"""
function update_initializeprob!(initprob, prob)
    pgetter = get_scimlfn(prob).initialization_data.metadata.oop_reconstruct_u0_p.pgetter
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
    op::SymmapT
    """
    The `guesses` used to construct the initialization.
    """
    guesses::SymmapT
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
    """
    The value of the `missing_guess_value` keyword indicating how to handle missing guesses.
    """
    missing_guess_value::MissingGuessValue.Type
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

function GetUpdatedU0(sys::AbstractSystem, initsys::AbstractSystem, op::AbstractDict; kwargs...)
    dvs = unknowns(sys)
    eqs = equations(sys)
    guessvars = trues(length(dvs))
    for (i, var) in enumerate(dvs)
        varval = get(op, var, COMMON_NOTHING)
        guessvars[i] = varval === COMMON_NOTHING || !SU.isconst(varval)
    end
    get_guessvars = iszero(count(guessvars)) ? nothing : CopyParamsByTemplate(initsys, dvs[guessvars]; kwargs...)
    get_initial_unknowns = getu(sys, Initial.(dvs))
    return GetUpdatedU0(guessvars, get_guessvars, get_initial_unknowns)
end

function (guu::GetUpdatedU0)(prob, initprob)
    buffer = guu.get_initial_unknowns(prob)
    if guu.get_guessvars !== nothing
        algebuf = view(buffer, guu.guessvars)
        copyto!(algebuf, guu.get_guessvars(initprob))
    end
    return buffer
end

struct SetInitialUnknowns{S}
    setter!::S
    idxs_in_initials::Vector{Int}
end

function SetInitialUnknowns(sys::AbstractSystem)
    initpars = Initial.(unknowns(sys))
    idxs_in_initials = Int[]
    sizehint!(idxs_in_initials, length(initpars))
    if is_split(sys)
        for par in initpars
            idx = parameter_index(sys, par)::ParameterIndex{SciMLStructures.Initials, Int}
            push!(idxs_in_initials, idx.idx)
        end
    else
        for par in initpars
            idx = parameter_index(sys, par)::Int
            push!(idxs_in_initials, idx)
        end
    end
    return SetInitialUnknowns(setu(sys, initpars), idxs_in_initials)
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
    PromoteToTunableEltype(observed, floatT)

Wraps an `initializeprob` observed function so its output array is promoted to an
eltype compatible with the current tunable parameters. Addresses the case where
the observed function is generated from fully constant RHS (e.g. `initialization_eqs
= [s ~ 0]`): the resulting `create_array(Array, nothing, …, 0, 0)` would otherwise
produce `Vector{Int64}`, which — when downstream `remake` reinstalls it as `u0` —
silently defeats ForwardDiff/Tracker/Measurements promotion of `u0`.

`floatT` is the static floor (the same `floatT` the rest of the construction pipeline
commits to, derived from the user's varmap). It guarantees integer RHS gets lifted
to a float without overriding the user's chosen precision (e.g. `Float32`). The
dynamic tunable eltype is read fresh from `parameter_values(nlsol)` on every call,
so a later `remake` with `ForwardDiff.Dual` parameters still wins via `promote_type`.
"""
struct PromoteToTunableEltype{F, floatT}
    observed::F
end

PromoteToTunableEltype(observed, ::Type{T}) where {T} =
    PromoteToTunableEltype{typeof(observed), T}(observed)

function (p::PromoteToTunableEltype{F, floatT})(nlsol) where {F, floatT}
    raw = p.observed(nlsol)
    raw isa AbstractArray || return raw
    isempty(raw) && return raw
    T = promote_type(eltype(raw), _tunable_eltype(parameter_values(nlsol)), floatT)
    return T === eltype(raw) ? raw : convert(AbstractArray{T}, raw)
end

_tunable_eltype(p::MTKParameters) = isempty(p.tunable) ? Bool : eltype(p.tunable)
function _tunable_eltype(p)
    if SciMLStructures.isscimlstructure(p)
        tun = SciMLStructures.canonicalize(SciMLStructures.Tunable(), p)[1]
        return isempty(tun) ? Bool : eltype(tun)
    end
    return Bool
end

"""
    $TYPEDSIGNATURES

Utility function present in `initializeprobmap` for iip problems. If the problem is `iip`,
then `u0` cannot be a struct-of-arrays (e.g. `Tracker.TrackedVector`). It needs to be
turned into an array-of-structs (`Vector{Tracker.TrackedReal{..}}`). Methods are added to
this function in extensions.
"""
__iip_u0_ad_wrapper(x) = x

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
        sys::AbstractSystem, iip::Bool, op::SymmapT, t, guesses;
        time_dependent_init = is_time_dependent(sys), u0_constructor = identity,
        p_constructor = identity, floatT = Float64, initialization_eqs = [],
        use_scc = true, eval_expression = false, eval_module = @__MODULE__,
        missing_guess_value = default_missing_guess_value(),
        # Intercept `expression` because we don't support it here yet
        implicit_dae = false, is_steadystateprob = false, expression = Val{false}, kwargs...
    )
    guesses = merge(ModelingToolkitBase.guesses(sys), todict(guesses))

    if t === nothing && is_time_dependent(sys)
        t = zero(floatT)
    end

    orig_op = copy(op)
    initializeprob = ModelingToolkitBase.InitializationProblem{iip}(
        sys, t, op; guesses, time_dependent_init, initialization_eqs, fast_path = true,
        use_scc, u0_constructor, p_constructor, eval_expression, eval_module,
        missing_guess_value, is_steadystateprob, kwargs...
    )
    initsys = initializeprob.f.sys::System
    needs_remake = false
    _u0 = state_values(initializeprob)
    if _u0 !== nothing
        if ArrayInterface.ismutable(_u0)
            __u0 = floatT.(_u0)
        else
            __u0 = similar_type(_u0, floatT)(_u0)
        end
        if eltype(__u0) != eltype(_u0)
            _u0 = __u0
            needs_remake = true
        end
    end
    initp = parameter_values(initializeprob)
    if is_split(sys)
        buffer, repack, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), initp)
        _initp = repack(floatT.(buffer))
        if !(initp.initials isa StaticVector{0})
            buffer, repack, _ = SciMLStructures.canonicalize(SciMLStructures.Initials(), _initp)
            _initp = repack(floatT.(buffer))
        end
        if eltype(_initp.tunable) != eltype(initp.tunable) || eltype(_initp.initials) != eltype(initp.initials)
            initp = _initp
            needs_remake = true
        end
    elseif initp isa AbstractArray
        if ArrayInterface.ismutable(initp)
            initp′ = similar(initp, floatT)
            if eltype(initp′) != eltype(initp)
                copyto!(initp′, initp)
                initp = initp′
                needs_remake = true
            end
        else
            initp′ = similar_type(initp, floatT)(initp)
            if eltype(initp′) != eltype(initp)
                initp = initp′
                needs_remake = true
            end
        end
    end
    if needs_remake
        initializeprob = remake(initializeprob; u0 = _u0, p = initp)
    end

    get_initial_unknowns = if time_dependent_init
        GetUpdatedU0(sys, initsys, op; eval_expression, eval_module, kwargs...)
    else
        nothing
    end
    meta = InitializationMetadata(
        orig_op,
        as_atomic_dict_with_defaults(Dict{SymbolicT, SymbolicT}(guesses), COMMON_NOTHING),
        Vector{Equation}(initialization_eqs),
        use_scc, time_dependent_init,
        ReconstructInitializeprob(
            sys, initsys; u0_constructor,
            p_constructor, eval_expression, eval_module, is_steadystateprob, kwargs...
        ),
        get_initial_unknowns, SetInitialUnknowns(sys), missing_guess_value
    )

    if time_dependent_init
        all_init_syms = Set(all_symbols(initializeprob))
        solved_unknowns = filter(var -> var in all_init_syms, unknowns(sys))
        if isempty(solved_unknowns)
            initializeprobmap = nothing
        else
            initializeprobmap = u0_constructor ∘ PromoteToTunableEltype(CopyParamsByTemplate(initializeprob.f.sys, solved_unknowns; eval_expression, eval_module, kwargs...), floatT)
            if iip
                initializeprobmap = __iip_u0_ad_wrapper ∘ initializeprobmap
            end
        end
    else
        initializeprobmap = nothing
    end

    punknowns = [
        p
            for p in all_variable_symbols(initializeprob)
            if is_parameter(sys, p)
    ]
    if initializeprobmap === nothing && isempty(punknowns)
        initializeprobpmap = nothing
    else
        initializeprobpmap = construct_initializeprobpmap(
            sys, initsys; p_constructor, eval_expression, eval_module, kwargs...
        )
    end

    # we still want the `initialization_data` because it helps with `remake`
    if initializeprobmap === nothing && initializeprobpmap === nothing
        update_initializeprob! = nothing
    else
        update_initializeprob! = ModelingToolkitBase.update_initializeprob!
    end

    missingvars = Set{SymbolicT}()
    temp_op = copy(op)
    for (k, v) in op
        v === COMMON_MISSING || continue
        push!(missingvars, k)
        delete!(temp_op, k)
    end
    binds = bindings(sys)
    if time_dependent_init
        for v in unknowns(sys)
            has_possibly_indexed_key(parent(binds), v) && continue
            val = get_possibly_indexed(op, v, COMMON_NOTHING)
            if !SU.isconst(val) || val === COMMON_NOTHING
                push!(missingvars, v)
            end
        end
        if implicit_dae
            for v in unknowns(sys)
                v = Differential(get_iv(sys))(v)
                ttv = default_toterm(v)
                if get_possibly_indexed(op, v, COMMON_NOTHING) === COMMON_NOTHING &&
                        get_possibly_indexed(op, ttv, COMMON_NOTHING) === COMMON_NOTHING &&
                        # FIXME: Derivatives of algebraic variables aren't present
                        (is_variable(initsys, ttv) || has_observed_with_lhs(initsys, ttv))
                    push!(missingvars, ttv)
                end
            end
        end
    end
    for v in get_all_discretes_fast(sys)
        has_possibly_indexed_key(parent(binds), v) && continue
        has_possibly_indexed_key(op, v) || push!(missingvars, v)
    end
    for (k, v) in binds
        v === COMMON_MISSING && !has_possibly_indexed_key(op, k) && push!(missingvars, k)
    end
    for p in as_atomic_array_set(parameters(sys))
        haskey(binds, p) && continue
        haskey(op, p) || push!(missingvars, p)
    end
    missingvars = collect(missingvars)

    for (i, v) in enumerate(unknowns(initsys))
        write_possibly_indexed_array!(temp_op, v, SConst(_u0[i]), COMMON_NOTHING)
    end
    add_observed!(initsys, temp_op)
    left_merge!(temp_op, ModelingToolkitBase.guesses(sys))
    subber = Symbolics.FixpointSubstituter{true}(AADSubWrapper(temp_op))
    for p in missingvars
        write_possibly_indexed_array!(op, p, subber(p), COMMON_NOTHING)
    end

    return (;
        initialization_data = SciMLBase.OverrideInitData(
            initializeprob, update_initializeprob!, initializeprobmap,
            initializeprobpmap; metadata = meta, is_update_oop = Val(true)
        ),
    )
end

function rm_union(::Type{T}) where {T}
    types = filter(!=(Nothing), Base.uniontypes(T))
    isempty(types) && return T
    return Core.apply_type(Union, types...)
end

"""
    $(TYPEDSIGNATURES)

Calculate the floating point type to use from the given `varmap` by looking at variables
with a constant value.
"""
function float_type_from_varmap(varmap, floatT = Bool)
    for (k, v) in varmap
        is_variable_floatingpoint(k) || continue
        SU.isconst(v) || symbolic_type(v) isa NotSymbolic || continue
        is_array_of_symbolics(v) && continue
        v = unwrap_const(v)
        if v isa AbstractArray
            # Remove union in case some elements of the array are `nothing`
            floatT = promote_type(floatT, rm_union(eltype(unwrap_const(v))))
        elseif v isa Number
            floatT = promote_type(floatT, typeof(unwrap_const(v)))
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
            u0ElType
        )
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
        elT = if symbolic_u0 && any(x -> x === nothing || symbolic_type(x) != NotSymbolic(), vals)
            nothing
        else
            floatT
        end
        return SymbolicUtils.Code.create_array(u0Type, elT, Val(1), Val(length(vals)), vals...)
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
        return SymbolicUtils.Code.create_array(
            pType, floatT, Val(ndims(vals)), Val(size(vals)), vals...
        )
    end
end

abstract type ProblemConstructionHook end

function operating_point_preprocess(sys::AbstractSystem, op; name = "operating_point")
    if op !== nothing && !(eltype(op) <: Pair) && !isempty(op)
        throw(
            ArgumentError(
                """
                The $name passed to the problem constructor must be a symbolic map.
                """
            )
        )
    end
    op = recursive_unwrap(anydict(op))
    symbols_to_symbolics!(sys, op)
    return op
end

function build_operating_point(sys::AbstractSystem, op; fast_path = false)
    if !fast_path
        op = operating_point_preprocess(sys, op)
    end
    # Replace `nothing`s with sentinels so that `left_merge!` thinks they're values
    # and doesn't override them. This is because explicit `nothing` values in `op`
    # should be considered as overrides for initial conditions in `ics`.
    map!(x -> @something(x, CommonSentinel()), values(op))
    op = as_atomic_dict_with_defaults(Dict{SymbolicT, SymbolicT}(op), COMMON_NOTHING)
    ics = add_toterms(initial_conditions(sys); replace = is_discrete_system(sys))
    left_merge!(op, ics)
    map!(values(op)) do v
        v === COMMON_SENTINEL && return COMMON_NOTHING
        Symbolics.isarraysymbolic(v) || return v
        any(Base.Fix2(===, COMMON_SENTINEL) ∘ Base.Fix1(getindex, v), SU.stable_eachindex(v)) || return v

        new_v = map(SU.stable_eachindex(v)) do i
            v[i] === COMMON_SENTINEL ? COMMON_NOTHING : v[i]
        end
        return SU.Const{VartypeT}(new_v)
    end
    filter!(Base.Fix2(!==, COMMON_NOTHING) ∘ last, op)
    return op
end

struct MissingNecessaryInitialConditionsError <: Exception
    missing_ics::NecessaryInitialConditionsT
end

const MISSING_NECESSARY_ICS_ERR_PRELUDE = """
Missing necessary initial conditions for some variables. The missing variables are listed \
below, along with the reason why the initial condition is required:
"""

function Base.showerror(io::IO, err::MissingNecessaryInitialConditionsError)
    println(io, MISSING_NECESSARY_ICS_ERR_PRELUDE)
    for (k, v) in err.missing_ics
        printstyled(io, k; bold = true)
        println(io, ": ", v)
    end
    return
end

function check_necessary_initial_conditions(sys::AbstractSystem, op::SymmapT)
    ics = get_necessary_initial_conditions(sys)
    isempty(ics) && return
    missing_ics = filter(!Base.Fix1(has_possibly_indexed_key, op) ∘ first, ics)
    isempty(missing_ics) || throw(MissingNecessaryInitialConditionsError(missing_ics))
    return
end

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
Base.@nospecializeinfer function process_SciMLProblem(
        @nospecialize(constructor), sys::AbstractSystem, @nospecialize(op);
        u0_eltype = nothing, u0_constructor = identity, p_constructor = identity,
        symbolic_u0 = false, kwargs...
    )
    u0Type = pType = typeof(op)
    op = operating_point_preprocess(sys, op)
    floatT = calculate_float_type(op, u0Type)
    u0_eltype = something(u0_eltype, floatT)
    u0_constructor = get_u0_constructor(u0_constructor, u0Type, floatT, symbolic_u0)
    p_constructor = get_p_constructor(p_constructor, pType, floatT)

    __process_SciMLProblem(constructor, sys, op; floatT, u0Type, u0_eltype, u0_constructor, p_constructor, symbolic_u0, kwargs...)
end

function __process_SciMLProblem(
        @nospecialize(constructor), sys::AbstractSystem, op::AnyDict;
        floatT, u0Type, u0_eltype,
        build_initializeprob = supports_initialization(sys),
        implicit_dae = false, t = nothing, guesses = AnyDict(),
        warn_initialize_determined = true, initialization_eqs = [],
        eval_expression = false, eval_module = @__MODULE__, fully_determined = nothing,
        check_initialization_units = false, tofloat = true,
        u0_constructor = identity, p_constructor = identity,
        check_length = true, symbolic_u0 = false, warn_cyclic_dependency = false,
        circular_dependency_max_cycle_length = length(all_symbols(sys)),
        circular_dependency_max_cycles = 10, initsys_mtkcompile_kwargs = (;),
        substitution_limit = 100, use_scc = true, time_dependent_init = is_time_dependent(sys),
        algebraic_only = false, missing_guess_value = default_missing_guess_value(),
        allow_incomplete = false, is_initializeprob = false, is_steadystateprob = false,
        return_operating_point = false,
        compiler_options::CompilerOptions = CompilerOptions(),
        init_compiler_options::CompilerOptions = CompilerOptions(),
        kwargs...
    )
    dvs = unknowns(sys)
    ps = parameters(sys; initial_parameters = true)
    iv = has_iv(sys) ? get_iv(sys) : nothing
    eqs = equations(sys)

    check_array_equations_unknowns(eqs, dvs)

    op = build_operating_point(sys, op; fast_path = true)

    check_inputmap_keys(sys, op)

    op = getmetadata(sys, ProblemConstructionHook, identity)(op)::SymmapT
    check_necessary_initial_conditions(sys, op)

    kwargs = NamedTuple(kwargs)

    add_initials!(sys, op)

    _sys = reverse_all_default_reversible_transformations(sys)
    obs = observed(_sys)

    guesses = operating_point_preprocess(sys, guesses; name = "guesses")

    if !is_time_dependent(sys) || is_initializesystem(sys)
        add_observed_equations!(op, obs, bindings(sys))
    end

    if build_initializeprob
        kws = maybe_build_initialization_problem(
            sys, constructor <: SciMLBase.AbstractSciMLFunction{true},
            op, t, guesses; initsys_mtkcompile_kwargs,
            warn_initialize_determined, initialization_eqs,
            eval_expression, eval_module, fully_determined,
            warn_cyclic_dependency, check_units = check_initialization_units,
            circular_dependency_max_cycle_length, circular_dependency_max_cycles, use_scc,
            algebraic_only, allow_incomplete, u0_constructor, p_constructor, floatT,
            time_dependent_init, missing_guess_value, is_steadystateprob, implicit_dae,
            compiler_options = init_compiler_options,
            kwargs...
        )

        kwargs = merge(kwargs, kws)
    end

    if t !== nothing && !(constructor <: Union{DDEFunction, SDDEFunction})
        op[iv] = t
    end

    binds = bindings(sys)
    # If we aren't building an initialization problem, we aren't respecting non-parameter
    # bindings anyway.
    if build_initializeprob
        no_override_merge_except_missing!(op, binds)
    else
        for p in bound_parameters(sys)
            haskey(op, p) && throw(ArgumentError("Cannot provide initial value for bound parameter $p."))
            op[p] = binds[p]
        end
        left_merge!(op, binds)
    end
    add_observed_equations!(op, obs)

    if warn_cyclic_dependency
        cycles = check_substitution_cycles(
            op, dvs; max_cycle_length = circular_dependency_max_cycle_length,
            max_cycles = circular_dependency_max_cycles
        )
        if !isempty(cycles)
            buffer = IOBuffer()
            for cycle in cycles
                println(buffer, cycle)
            end
            msg = String(take!(buffer))
            @warn "Cycles in unknowns:\n$msg"
        end
    end

    ir = get_irstructure(sys)
    if is_initializeprob
        u0 = varmap_to_vars(
            op, dvs; ir, buffer_eltype = u0_eltype, container_type = u0Type,
            allow_symbolic = symbolic_u0, is_initializeprob, substitution_limit,
            missing_values = missing_guess_value
        )
    else
        u0 = varmap_to_vars(
            op, dvs; ir, buffer_eltype = u0_eltype, container_type = u0Type,
            allow_symbolic = symbolic_u0, is_initializeprob, substitution_limit
        )
    end
    if u0 !== nothing
        u0 = u0_constructor(u0)
    end

    check_eqs_u0(eqs, dvs, u0; check_length, kwargs...)

    if warn_cyclic_dependency
        cycles = check_substitution_cycles(
            op, ps; max_cycle_length = circular_dependency_max_cycle_length,
            max_cycles = circular_dependency_max_cycles
        )
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
        p = MTKParameters(sys, op; floatT = floatT, p_constructor, fast_path = true)
    else
        p = p_constructor(varmap_to_vars(op, ps; tofloat, container_type = u0Type))
    end

    if implicit_dae
        ddvs = map(default_toterm ∘ Differential(iv), dvs)
        du0 = varmap_to_vars(
            op, ddvs; toterm = default_toterm,
            tofloat
        )
        kwargs = merge(kwargs, (; ddvs))
    else
        du0 = nothing
    end

    if constructor <: NonlinearFunction && length(dvs) != length(eqs)
        kwargs = merge(
            kwargs,
            (;
                resid_prototype = u0_constructor(
                    calculate_resid_prototype(
                        length(eqs), u0, p
                    )
                ),
            )
        )
    end

    f = constructor(
        sys; u0 = u0, p = p, t = t,
        eval_expression = eval_expression,
        eval_module = eval_module,
        compiler_options,
        kwargs...
    )
    if return_operating_point
        return implicit_dae ? (f, du0, u0, p, op) : (f, u0, p, op)
    else
        return implicit_dae ? (f, du0, u0, p) : (f, u0, p)
    end
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

    return if !isempty(badvarkeys)
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
    return println(io, join(e.vars, ", "))
end

function SciMLBase.detect_cycles(sys::AbstractSystem, varmap::Dict{Any, Any}, vars)
    varmap = AnyDict(unwrap(k) => unwrap(v) for (k, v) in varmap)
    vars = map(unwrap, vars)
    cycles = check_substitution_cycles(varmap, vars)
    return !isempty(cycles)
end

function process_kwargs(
        sys::System; expression = Val{false}, callback = nothing,
        eval_expression = false, eval_module = @__MODULE__,
        _skip_events = false, _skip_tstops = false, tspan = nothing, kwargs...
    )
    kwargs = filter_kwargs(kwargs)
    kwargs1 = (;)

    if is_time_dependent(sys)
        if expression == Val{false} && !_skip_events
            cbs = process_events(
                sys; callback, eval_expression, eval_module, tspan, kwargs...
            )
            if cbs !== nothing
                kwargs1 = merge(kwargs1, (callback = cbs,))
            end
        end

        if !_skip_tstops
            tstops = SymbolicTstops(
                sys, GeneratedFunctionOptions(; expression, eval_expression, eval_module)
            )
            if tstops !== nothing
                kwargs1 = merge(kwargs1, (; tstops))
            end
        end
    end

    return merge(kwargs1, kwargs)
end

function filter_kwargs(kwargs)
    kwargs = Dict(kwargs)
    for key in keys(kwargs)
        key in DiffEqBase.allowedkeywords || delete!(kwargs, key)
    end
    return pairs(NamedTuple(kwargs))
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

function SymbolicTstops(sys::AbstractSystem, opts::GeneratedFunctionOptions)
    expression = expression_val(opts)
    (; eval_expression, eval_module) = opts
    tstops = symbolic_tstops(sys)
    isempty(tstops) && return nothing
    t0 = gensym(:t0)
    t1 = gensym(:t1)
    tstops = map(tstops) do val
        if is_array_of_symbolics(val) || val isa AbstractArray
            collect(val)
        else
            term(:, term(+, t0, unwrap(val)), unwrap(val), t1; type = Vector{Real})
        end
    end
    rps = reorder_parameters(sys)
    tstops,
        _ = build_function_wrapper(
        sys, Symbolics.SConst(tstops),
        [rps; Any[t0, t1]],
        BuildFunctionWrapperOptions(;
            p_start = 1, p_end = length(rps),
            codegen_function_options = Symbolics.CodegenFunctionOptions(;
                expression = Val{true}, force_SA = true
            )
        )
    )
    tstops = GeneratedFunctionWrapper{(1, 3, is_split(sys))}(
        expression, tstops, nothing; eval_expression, eval_module
    )

    if expression == Val{true}
        return :($SymbolicTstops($tstops))
    else
        return SymbolicTstops(tstops)
    end
end

# Backward-compatibility keyword method. The positional `opts::GeneratedFunctionOptions`
# method above is the primary; this wrapper preserves the historical keyword API. It must
# be defined after `struct SymbolicTstops` so it registers as a method on that binding.
function SymbolicTstops(
        sys::AbstractSystem; expression = Val{false}, eval_expression = false,
        eval_module = @__MODULE__
    )
    return SymbolicTstops(
        sys, GeneratedFunctionOptions(; expression, eval_expression, eval_module)
    )
end

"""
    Both

Sentinel type used as the `iip` type parameter in problem constructors when the caller
did not explicitly specify in-place vs. out-of-place behavior. The actual value is
resolved at construction time: `iip = false` when the operating-point `op` is a
`StaticArray`, and `iip = true` otherwise.
"""
struct Both end

"""
    resolve_iip(iip, op)

Resolve the `iip` type parameter for a problem constructor. When `iip` is `Both`,
returns `false` if `op` is a `StaticArray` and `true` otherwise. For any other value
of `iip`, returns `iip` unchanged.
"""
resolve_iip(iip, @nospecialize(op)) = iip
resolve_iip(::Type{Both}, @nospecialize(op)) = !(op isa StaticArray)

"""
    $(TYPEDSIGNATURES)

Macro for writing problem/function constructors. Expects a function definition with type
parameters for `iip` and `specialize`. Generates fallbacks with
`specialize = SciMLBase.AutoSpecialize` and `iip = Both` (resolved at construction time).
"""
# Unwrap `@nospecialize(arg)` to get the underlying argument expression.
# Returns the argument unchanged if not wrapped in @nospecialize.
function _unwrap_nospecialize(arg)
    if Meta.isexpr(arg, :macrocall) && length(arg.args) >= 3 &&
            arg.args[1] in (Symbol("@nospecialize"), GlobalRef(Base, Symbol("@nospecialize")))
        return arg.args[3]
    end
    return arg
end

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

    # Create signature args with @nospecialize stripped (for fallback function signatures)
    sig_args = map(args) do arg
        unwrapped = _unwrap_nospecialize(arg)
        # Handle :parameters specially - unwrap each kwarg inside
        if Meta.isexpr(unwrapped, :parameters)
            new_params = map(unwrapped.args) do kwarg
                kw = _unwrap_nospecialize(kwarg)
                # Convert :(=) to :kw if needed
                if Meta.isexpr(kw, :(=))
                    Expr(:kw, kw.args...)
                else
                    kw
                end
            end
            return Expr(:parameters, new_params...)
        end
        return unwrapped
    end

    # arguments to call with (for forwarding calls)
    call_args = map(sig_args) do arg
        # keyword args are in `Expr(:parameters)` so any `Expr(:kw)` here
        # are optional positional arguments. Analyze `:(f(a, b = 1; k = 1, l...))`
        # to understand
        Meta.isexpr(arg, :kw) && return arg.args[1]
        return arg
    end
    call_kwargs = map(call_args[1].args) do arg
        Meta.isexpr(arg, :...) && return arg
        @assert Meta.isexpr(arg, :kw) "Expected keyword argument, got $(arg)"
        return Expr(:kw, arg.args[1], arg.args[1])
    end
    call_args[1] = Expr(:parameters, call_kwargs...)

    @assert Meta.isexpr(fnname_curly, :curly)
    # fnname_name is `ODEProblem`
    # curly_args is `iip, spec`
    fnname_name, curly_args... = fnname_curly.args
    @assert curly_args == where_args

    # callexpr_iip is `ODEProblem{iip, AutoSpecialize}(call_args...)`
    callexpr_iip = Expr(
        :call, Expr(:curly, fnname_name, curly_args[1], SciMLBase.AutoSpecialize), call_args...
    )
    # `ODEProblem{iip}`
    fnname_iip = Expr(:curly, fnname_name, curly_args[1])
    # `ODEProblem{iip}(sig_args...)` - use sig_args (no @nospecialize) for fallback signature
    fncall_iip = Expr(:call, fnname_iip, sig_args...)
    # ODEProblem{iip}(sig_args...) where {iip}
    fnwhere_iip = Expr(:where, fncall_iip, where_args[1])
    fn_iip = Expr(:function, fnwhere_iip, callexpr_iip)

    # Problem constructors default to `Both` (iip resolved at construction time).
    # Function constructors keep `true` as the iip default.
    is_problem = occursin("Problem", string(fnname_name))
    default_iip = is_problem ? Both : true
    callexpr_base = Expr(:call, Expr(:curly, fnname_name, default_iip), call_args...)
    # `ODEProblem(sig_args...)` - use sig_args for fallback signature
    fncall_base = Expr(:call, fnname_name, sig_args...)
    fn_base = Expr(:function, fncall_base, callexpr_base)

    # The StaticArray-specific fallback is no longer needed: `Both` defers the
    # iip decision to the body of the fully-parameterized method.
    fn_sarr = nothing
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
        block, kws::NamedTuple; head = :(=), filter = Returns(true)
    )
    for (k, v) in pairs(kws)
        filter(k) || continue
        push!(block.args, Expr(head, k, v))
    end
    return
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
    return build_scimlfn_expr(T, args; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Construct SciMLFunction `T` with positional arguments `args` and keywords `kwargs`.
"""
function maybe_codegen_scimlfn(::Type{Val{false}}, ::Type{T}, args::NamedTuple; kwargs...) where {T}
    @nospecialize args kwargs
    return T(args...; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Return an expression constructing SciMLProblem `T` with positional arguments `args`
and keywords `kwargs`.
"""
function maybe_codegen_scimlproblem(::Type{Val{true}}, T, args::NamedTuple; kwargs...)
    return build_scimlproblem_expr(T, args; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Construct SciMLProblem `T` with positional arguments `args` and keywords `kwargs`.
"""
function maybe_codegen_scimlproblem(::Type{Val{false}}, ::Type{T}, args::NamedTuple; kwargs...) where {T}
    @nospecialize args kwargs
    # Call `remake` so it runs initialization if it is trivial
    # Use `@invokelatest` to avoid world-age issues with `eval_expression = true`
    return @invokelatest remake(T(args...; kwargs...))
end

"""
    $(TYPEDSIGNATURES)

Return the `u0` vector for the given system `sys` and variable-value mapping `varmap`. All
keyword arguments are forwarded to [`varmap_to_vars`](@ref).
"""
function get_u0(sys::AbstractSystem, varmap; kwargs...)
    op = build_operating_point(sys, varmap)
    binds = bindings(sys)
    no_override_merge_except_missing!(op, binds)
    obs = observed(reverse_all_default_reversible_transformations(sys))
    add_observed_equations!(op, obs)

    return varmap_to_vars(op, unknowns(sys); kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Return the `p` object for the given system `sys` and variable-value mapping `varmap`. All
keyword arguments are forwarded to [`MTKParameters`](@ref) for split systems and
[`varmap_to_vars`](@ref) for non-split systems.
"""
function get_p(sys::AbstractSystem, varmap; split = is_split(sys), kwargs...)
    op = build_operating_point(sys, varmap)
    binds = bindings(sys)
    no_override_merge_except_missing!(op, binds)
    add_initials!(sys, op)
    obs = observed(reverse_all_default_reversible_transformations(sys))
    add_observed_equations!(op, obs)

    return if split
        MTKParameters(sys, op; fast_path = true, kwargs...)
    else
        varmap_to_vars(op, parameters(sys; initial_parameters = true); kwargs...)
    end
end
