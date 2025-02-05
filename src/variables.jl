struct VariableUnit end
struct VariableConnectType end
struct VariableNoiseType end
struct VariableInput end
struct VariableOutput end
struct VariableIrreducible end
struct VariableStatePriority end
struct VariableMisc end
Symbolics.option_to_metadata_type(::Val{:unit}) = VariableUnit
Symbolics.option_to_metadata_type(::Val{:connect}) = VariableConnectType
Symbolics.option_to_metadata_type(::Val{:input}) = VariableInput
Symbolics.option_to_metadata_type(::Val{:output}) = VariableOutput
Symbolics.option_to_metadata_type(::Val{:irreducible}) = VariableIrreducible
Symbolics.option_to_metadata_type(::Val{:state_priority}) = VariableStatePriority
Symbolics.option_to_metadata_type(::Val{:misc}) = VariableMisc

"""
    dump_variable_metadata(var)

Return all the metadata associated with symbolic variable `var` as a `NamedTuple`.

```@example
using ModelingToolkit

@parameters p::Int [description = "My description", bounds = (0.5, 1.5)]
ModelingToolkit.dump_variable_metadata(p)
```
"""
function dump_variable_metadata(var)
    uvar = unwrap(var)
    variable_source, name = Symbolics.getmetadata(
        uvar, VariableSource, (:unknown, :unknown))
    type = symtype(uvar)
    if type <: AbstractArray
        shape = Symbolics.shape(var)
        if shape == ()
            shape = nothing
        end
    else
        shape = nothing
    end
    unit = getunit(uvar)
    connect = getconnect(uvar)
    input = isinput(uvar) || nothing
    output = isoutput(uvar) || nothing
    irreducible = isirreducible(var)
    state_priority = Symbolics.getmetadata(uvar, VariableStatePriority, nothing)
    misc = getmisc(uvar)
    bounds = hasbounds(uvar) ? getbounds(uvar) : nothing
    desc = getdescription(var)
    if desc == ""
        desc = nothing
    end
    default = hasdefault(uvar) ? getdefault(uvar) : nothing
    guess = getguess(uvar)
    disturbance = isdisturbance(uvar) || nothing
    tunable = istunable(uvar, isparameter(uvar))
    dist = getdist(uvar)
    variable_type = getvariabletype(uvar)

    meta = (
        var = var,
        variable_source,
        name,
        variable_type,
        shape,
        unit,
        connect,
        input,
        output,
        irreducible,
        state_priority,
        misc,
        bounds,
        desc,
        guess,
        disturbance,
        tunable,
        dist,
        type,
        default
    )

    return NamedTuple(k => v for (k, v) in pairs(meta) if v !== nothing)
end

### Connect
abstract type AbstractConnectType end
struct Equality <: AbstractConnectType end # Equality connection
struct Flow <: AbstractConnectType end     # sum to 0
struct Stream <: AbstractConnectType end   # special stream connector

"""
    getconnect(x)

Get the connect type of x. See also [`hasconnect`](@ref).
"""
getconnect(x) = getconnect(unwrap(x))
getconnect(x::Symbolic) = Symbolics.getmetadata(x, VariableConnectType, nothing)
"""
    hasconnect(x)

Determine whether variable `x` has a connect type. See also [`getconnect`](@ref).
"""
hasconnect(x) = getconnect(x) !== nothing
function setconnect(x, t::Type{T}) where {T <: AbstractConnectType}
    setmetadata(x, VariableConnectType, t)
end

### Input, Output, Irreducible 
isvarkind(m, x::Union{Num, Symbolics.Arr}) = isvarkind(m, value(x))
function isvarkind(m, x)
    iskind = getmetadata(x, m, nothing)
    iskind !== nothing && return iskind
    x = getparent(x, x)
    getmetadata(x, m, false)
end

setinput(x, v::Bool) = setmetadata(x, VariableInput, v)
setoutput(x, v::Bool) = setmetadata(x, VariableOutput, v)
setio(x, i::Bool, o::Bool) = setoutput(setinput(x, i), o)

isinput(x) = isvarkind(VariableInput, x)
isoutput(x) = isvarkind(VariableOutput, x)

# Before the solvability check, we already have handled IO variables, so
# irreducibility is independent from IO.
isirreducible(x) = isvarkind(VariableIrreducible, x)
setirreducible(x, v::Bool) = setmetadata(x, VariableIrreducible, v)
state_priority(x) = convert(Float64, getmetadata(x, VariableStatePriority, 0.0))::Float64

function default_toterm(x)
    if iscall(x) && (op = operation(x)) isa Operator
        if !(op isa Differential)
            if op isa Shift && op.steps < 0
                return x
            end
            x = normalize_to_differential(op)(arguments(x)...)
        end
        Symbolics.diff2term(x)
    else
        x
    end
end

"""
$(SIGNATURES)

Takes a list of pairs of `variables=>values` and an ordered list of variables
and creates the array of values in the correct order with default values when
applicable.
"""
function varmap_to_vars(varmap, varlist; defaults = Dict(), check = true,
        toterm = default_toterm, promotetoconcrete = nothing,
        tofloat = true, use_union = true)
    varlist = collect(map(unwrap, varlist))

    # Edge cases where one of the arguments is effectively empty.
    is_incomplete_initialization = varmap isa DiffEqBase.NullParameters ||
                                   varmap === nothing
    if is_incomplete_initialization || isempty(varmap)
        if isempty(defaults)
            if !is_incomplete_initialization && check
                isempty(varlist) || throw(MissingVariablesError(varlist))
            end
            return nothing
        else
            varmap = Dict()
        end
    end

    # We respect the input type if it's a static array
    # otherwise canonicalize to a normal array
    # container_type = T <: Union{Dict,Tuple} ? Array : T
    if varmap isa StaticArray
        container_type = typeof(varmap)
    else
        container_type = Array
    end

    vals = if eltype(varmap) <: Pair # `varmap` is a dict or an array of pairs
        varmap = todict(varmap)
        _varmap_to_vars(varmap, varlist; defaults = defaults, check = check,
            toterm = toterm)
    else # plain array-like initialization
        varmap
    end

    promotetoconcrete === nothing && (promotetoconcrete = container_type <: AbstractArray)
    if promotetoconcrete
        vals = promote_to_concrete(vals; tofloat = tofloat, use_union = use_union)
    end

    if isempty(vals)
        return nothing
    elseif container_type <: Tuple
        (vals...,)
    else
        SymbolicUtils.Code.create_array(container_type, eltype(vals), Val{1}(),
            Val(length(vals)), vals...)
    end
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
    println(io, e.vars)
end

function _varmap_to_vars(varmap::Dict, varlist; defaults = Dict(), check = false,
        toterm = Symbolics.diff2term, initialization_phase = false)
    varmap = canonicalize_varmap(varmap; toterm)
    defaults = canonicalize_varmap(defaults; toterm)
    varmap = merge(defaults, varmap)
    values = Dict()

    T = Union{}
    for var in varlist
        var = unwrap(var)
        val = unwrap(fixpoint_sub(var, varmap; operator = Symbolics.Operator))
        if !isequal(val, var)
            values[var] = val
        end
    end
    missingvars = setdiff(varlist, collect(keys(values)))
    check && (isempty(missingvars) || throw(MissingVariablesError(missingvars)))
    return [values[unwrap(var)] for var in varlist]
end

function varmap_with_toterm(varmap; toterm = Symbolics.diff2term)
    return merge(todict(varmap), Dict(toterm(unwrap(k)) => v for (k, v) in varmap))
end

function canonicalize_varmap(varmap; toterm = Symbolics.diff2term)
    new_varmap = Dict()
    for (k, v) in varmap
        k = unwrap(k)
        v = unwrap(v)
        new_varmap[k] = v
        new_varmap[toterm(k)] = v
        if Symbolics.isarraysymbolic(k) && Symbolics.shape(k) !== Symbolics.Unknown()
            for i in eachindex(k)
                new_varmap[k[i]] = v[i]
                new_varmap[toterm(k[i])] = v[i]
            end
        end
    end
    return new_varmap
end

@noinline function throw_missingvars(vars)
    throw(ArgumentError("$vars are missing from the variable map."))
end

struct IsHistory end
ishistory(x) = ishistory(unwrap(x))
ishistory(x::Symbolic) = getmetadata(x, IsHistory, false)
hist(x, t) = wrap(hist(unwrap(x), t))
function hist(x::Symbolic, t)
    setmetadata(
        toparam(maketerm(typeof(x), operation(x), [unwrap(t)], metadata(x))),
        IsHistory, true)
end

## Bounds ======================================================================
struct VariableBounds end
Symbolics.option_to_metadata_type(::Val{:bounds}) = VariableBounds

"""
    getbounds(x)

Get the bounds associated with symbolic variable `x`.
Create parameters with bounds like this

```
@parameters p [bounds=(-1, 1)]
```
"""
function getbounds(x::Union{Num, Symbolics.Arr, SymbolicUtils.Symbolic})
    x = unwrap(x)
    p = Symbolics.getparent(x, nothing)
    if p === nothing
        bounds = Symbolics.getmetadata(x, VariableBounds, (-Inf, Inf))
        if symbolic_type(x) == ArraySymbolic() && Symbolics.shape(x) != Symbolics.Unknown()
            bounds = map(bounds) do b
                b isa AbstractArray && return b
                return fill(b, size(x))
            end
        end
    else
        # if we reached here, `x` is the result of calling `getindex`
        bounds = @something Symbolics.getmetadata(x, VariableBounds, nothing) getbounds(p)
        idxs = arguments(x)[2:end]
        bounds = map(bounds) do b
            if b isa AbstractArray
                if Symbolics.shape(p) != Symbolics.Unknown() && size(p) != size(b)
                    throw(DimensionMismatch("Expected array variable $p with shape $(size(p)) to have bounds of identical size. Found $bounds of size $(size(bounds))."))
                end
                return b[idxs...]
            elseif symbolic_type(x) == ArraySymbolic()
                return fill(b, size(x))
            else
                return b
            end
        end
    end
    return bounds
end

"""
    hasbounds(x)

Determine whether symbolic variable `x` has bounds associated with it.
See also [`getbounds`](@ref).
"""
function hasbounds(x)
    b = getbounds(x)
    any(isfinite.(b[1]) .|| isfinite.(b[2]))
end

## Disturbance =================================================================
struct VariableDisturbance end
Symbolics.option_to_metadata_type(::Val{:disturbance}) = VariableDisturbance

isdisturbance(x::Num) = isdisturbance(Symbolics.unwrap(x))

"""
    isdisturbance(x)

Determine whether symbolic variable `x` is marked as a disturbance input.
"""
function isdisturbance(x)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, VariableDisturbance, false)
end

function disturbances(sys)
    [filter(isdisturbance, unknowns(sys)); filter(isdisturbance, parameters(sys))]
end

## Tunable =====================================================================
struct VariableTunable end
Symbolics.option_to_metadata_type(::Val{:tunable}) = VariableTunable

istunable(x::Num, args...) = istunable(Symbolics.unwrap(x), args...)

"""
    istunable(x, default = true)

Determine whether symbolic variable `x` is marked as a tunable for an automatic tuning algorithm.

`default` indicates whether variables without `tunable` metadata are to be considered tunable or not.

Create a tunable parameter by

```
@parameters u [tunable=true]
```

See also [`tunable_parameters`](@ref), [`getbounds`](@ref)
"""
function istunable(x, default = true)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, VariableTunable, default)
end

## Dist ========================================================================
struct VariableDistribution end
Symbolics.option_to_metadata_type(::Val{:dist}) = VariableDistribution
getdist(x::Num) = getdist(Symbolics.unwrap(x))

"""
    getdist(x)

Get the probability distribution associated with symbolic variable `x`. If no distribution
is associated with `x`, `nothing` is returned.
Create parameters with associated distributions like this

```julia
using Distributions
d = Normal(0, 1)
@parameters u [dist = d]
hasdist(u) # true
getdist(u) # retrieve distribution
```
"""
function getdist(x)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, VariableDistribution, nothing)
end

"""
    hasdist(x)

Determine whether symbolic variable `x` has a probability distribution associated with it.
"""
function hasdist(x)
    b = getdist(x)
    b !== nothing
end

## System interface

"""
    tunable_parameters(sys, p = parameters(sys); default=true)

Get all parameters of `sys` that are marked as `tunable`.

Keyword argument `default` indicates whether variables without `tunable` metadata are to be considered tunable or not.

Create a tunable parameter by

```
@parameters u [tunable=true]
```

For systems created with `split = true` (the default) and `default = true` passed to this function, the order
of parameters returned is the order in which they are stored in the tunables portion of `MTKParameters`. Note
that array variables will not be scalarized. To obtain the flattened representation of the tunables portion,
call `Symbolics.scalarize(tunable_parameters(sys))` and concatenate the resulting arrays.

See also [`getbounds`](@ref), [`istunable`](@ref), [`MTKParameters`](@ref), [`complete`](@ref)
"""
function tunable_parameters(sys, p = parameters(sys); default = true)
    filter(x -> istunable(x, default), p)
end

"""
    getbounds(sys::ModelingToolkit.AbstractSystem, p = parameters(sys))

Returns a dict with pairs `p => (lb, ub)` mapping parameters of `sys` to lower and upper bounds.
Create parameters with bounds like this

```
@parameters p [bounds=(-1, 1)]
```

To obtain unknown variable bounds, call `getbounds(sys, unknowns(sys))`
"""
function getbounds(sys::ModelingToolkit.AbstractSystem, p = parameters(sys))
    Dict(p .=> getbounds.(p))
end

"""
    lb, ub = getbounds(p::AbstractVector)

Return vectors of lower and upper bounds of parameter vector `p`.
Create parameters with bounds like this

```
@parameters p [bounds=(-1, 1)]
```

See also [`tunable_parameters`](@ref), [`hasbounds`](@ref)
"""
function getbounds(p::AbstractVector)
    bounds = getbounds.(p)
    lb = first.(bounds)
    ub = last.(bounds)
    (; lb, ub)
end

## Description =================================================================
struct VariableDescription end
Symbolics.option_to_metadata_type(::Val{:description}) = VariableDescription

getdescription(x::Num) = getdescription(Symbolics.unwrap(x))
getdescription(x::Symbolics.Arr) = getdescription(Symbolics.unwrap(x))
"""
    getdescription(x)

Return any description attached to variables `x`. If no description is attached, an empty string is returned.
"""
function getdescription(x)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, VariableDescription, "")
end

function hasdescription(x)
    getdescription(x) != ""
end

## Brownian
"""
    tobrownian(s::Sym)

Maps the brownianiable to an unknown.
"""
tobrownian(s::Symbolic) = setmetadata(s, MTKVariableTypeCtx, BROWNIAN)
tobrownian(s::Num) = Num(tobrownian(value(s)))
isbrownian(s) = getvariabletype(s) === BROWNIAN

"""
$(SIGNATURES)

Define one or more Brownian variables.
"""
macro brownian(xs...)
    all(
        x -> x isa Symbol || Meta.isexpr(x, :call) && x.args[1] == :$ || Meta.isexpr(x, :$),
        xs) ||
        error("@brownian only takes scalar expressions!")
    Symbolics._parse_vars(:brownian,
        Real,
        xs,
        tobrownian) |> esc
end

## Guess ======================================================================
struct VariableGuess end
Symbolics.option_to_metadata_type(::Val{:guess}) = VariableGuess
getguess(x::Union{Num, Symbolics.Arr}) = getguess(Symbolics.unwrap(x))

"""
    getguess(x)

Get the guess for the initial value associated with symbolic variable `x`.
Create variables with a guess like this

```
@variables x [guess=1]
```
"""
function getguess(x)
    Symbolics.getmetadata(x, VariableGuess, nothing)
end

"""
    setguess(x, v)

Set the guess for the initial value associated with symbolic variable `x` to `v`.
See also [`hasguess`](@ref).
"""
function setguess(x, v)
    Symbolics.setmetadata(x, VariableGuess, v)
end

"""
    hasguess(x)

Determine whether symbolic variable `x` has a guess associated with it.
See also [`getguess`](@ref).
"""
function hasguess(x)
    getguess(x) !== nothing
end

function get_default_or_guess(x)
    if hasdefault(x) && !((def = getdefault(x)) isa Equation)
        return def
    else
        return getguess(x)
    end
end

## Miscellaneous metadata ======================================================================
"""
    getmisc(x)

Fetch any miscellaneous data associated with symbolic variable `x`.
See also [`hasmisc(x)`](@ref).
"""
getmisc(x) = getmisc(unwrap(x))
getmisc(x::Symbolic) = Symbolics.getmetadata(x, VariableMisc, nothing)
"""
    hasmisc(x)

Determine whether a symbolic variable `x` has misc
metadata associated with it. 

See also [`getmisc(x)`](@ref).
"""
hasmisc(x) = getmisc(x) !== nothing
setmisc(x, miscdata) = setmetadata(x, VariableMisc, miscdata)

## Units ======================================================================
"""
    getunit(x)

Fetch the unit associated with variable `x`. This function is a metadata getter for an individual variable, while `get_unit` is used for unit inference on more complicated sdymbolic expressions.
"""
getunit(x) = getunit(unwrap(x))
getunit(x::Symbolic) = Symbolics.getmetadata(x, VariableUnit, nothing)
"""
    hasunit(x)

Check if the variable `x` has a unit.
"""
hasunit(x) = getunit(x) !== nothing
