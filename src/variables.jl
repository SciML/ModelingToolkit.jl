struct VariableUnit end
struct VariableConnectType end
struct VariableNoiseType end
struct VariableInput end
struct VariableOutput end
struct VariableIrreducible end
struct VariableStatePriority end
Symbolics.option_to_metadata_type(::Val{:unit}) = VariableUnit
Symbolics.option_to_metadata_type(::Val{:connect}) = VariableConnectType
Symbolics.option_to_metadata_type(::Val{:noise}) = VariableNoiseType
Symbolics.option_to_metadata_type(::Val{:input}) = VariableInput
Symbolics.option_to_metadata_type(::Val{:output}) = VariableOutput
Symbolics.option_to_metadata_type(::Val{:irreducible}) = VariableIrreducible
Symbolics.option_to_metadata_type(::Val{:state_priority}) = VariableStatePriority

abstract type AbstractConnectType end
struct Equality <: AbstractConnectType end # Equality connection
struct Flow <: AbstractConnectType end     # sum to 0
struct Stream <: AbstractConnectType end   # special stream connector

isvarkind(m, x::Num) = isvarkind(m, value(x))
function isvarkind(m, x)
    iskind = getmetadata(x, m, nothing)
    iskind !== nothing && return iskind
    x = getparent(x, x)
    getmetadata(x, m, false)
end

setinput(x, v) = setmetadata(x, VariableInput, v)
setoutput(x, v) = setmetadata(x, VariableOutput, v)
setio(x, i, o) = setoutput(setinput(x, i), o)
isinput(x) = isvarkind(VariableInput, x)
isoutput(x) = isvarkind(VariableOutput, x)
# Before the solvability check, we already have handled IO varibales, so
# irreducibility is independent from IO.
isirreducible(x) = isvarkind(VariableIrreducible, x)
state_priority(x) = convert(Float64, getmetadata(x, VariableStatePriority, 0.0))::Float64

const DEFAULT_EMPTY_DICT = Dict()

"""
$(SIGNATURES)

Takes a list of pairs of `variables=>values` and an ordered list of variables
and creates the array of values in the correct order with default values when
applicable.
"""
function varmap_to_vars(varmap, varlist; defaults = Dict(), check = true,
                        toterm = Symbolics.diff2term, promotetoconcrete = nothing,
                        tofloat = true, use_union = false)
    varlist = collect(map(unwrap, varlist))

    # Edge cases where one of the arguments is effectively empty.
    is_incomplete_initialization = varmap isa DiffEqBase.NullParameters ||
                                   varmap === nothing
    if is_incomplete_initialization || isempty(varmap)
        if isempty(defaults)
            if !is_incomplete_initialization && check
                isempty(varlist) || throw_missingvars(varlist)
            end
            return nothing
        else
            varmap = Dict()
        end
    end

    T = typeof(varmap)
    # We respect the input type
    container_type = T <: Dict ? Array : T

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

function _varmap_to_vars(varmap::Dict, varlist; defaults = Dict(), check = false,
                         toterm = Symbolics.diff2term)
    varmap = merge(defaults, varmap) # prefers the `varmap`
    varmap = Dict(toterm(value(k)) => value(varmap[k]) for k in keys(varmap))
    # resolve symbolic parameter expressions
    for (p, v) in pairs(varmap)
        varmap[p] = fixpoint_sub(v, varmap)
    end

    missingvars = setdiff(varlist, collect(keys(varmap)))
    check && (isempty(missingvars) || throw_missingvars(missingvars))

    out = [varmap[var] for var in varlist]
end

@noinline function throw_missingvars(vars)
    throw(ArgumentError("$vars are missing from the variable map."))
end

"""
$(SIGNATURES)

Intercept the call to `handle_varmap` and convert it to an ordered list if the user has
  ModelingToolkit loaded, and the problem has a symbolic origin.
"""
function SciMLBase.handle_varmap(varmap, sys::AbstractSystem; field = :states, kwargs...)
    out = varmap_to_vars(varmap, getfield(sys, field); kwargs...)
    return out
end

struct IsHistory end
ishistory(x) = ishistory(unwrap(x))
ishistory(x::Symbolic) = getmetadata(x, IsHistory, false)
hist(x, t) = wrap(hist(unwrap(x), t))
function hist(x::Symbolic, t)
    setmetadata(toparam(similarterm(x, operation(x), [unwrap(t)], metadata = metadata(x))),
                IsHistory, true)
end

## Bounds ======================================================================
struct VariableBounds end
Symbolics.option_to_metadata_type(::Val{:bounds}) = VariableBounds
getbounds(x::Num) = getbounds(Symbolics.unwrap(x))

"""
    getbounds(x)

Get the bounds associated with symbolc variable `x`.
Create parameters with bounds like this
```
@parameters p [bounds=(-1, 1)]
```
"""
function getbounds(x)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, VariableBounds, (-Inf, Inf))
end

"""
    hasbounds(x)

Determine whether or not symbolic variable `x` has bounds associated with it.
See also [`getbounds`](@ref).
"""
function hasbounds(x)
    b = getbounds(x)
    isfinite(b[1]) || isfinite(b[2])
end

## Disturbance =================================================================
struct VariableDisturbance end
Symbolics.option_to_metadata_type(::Val{:disturbance}) = VariableDisturbance

isdisturbance(x::Num) = isdisturbance(Symbolics.unwrap(x))

"""
    isdisturbance(x)

Determine whether or not symbolic variable `x` is marked as a disturbance input.
"""
function isdisturbance(x)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, VariableDisturbance, false)
end

function disturbances(sys)
    [filter(isdisturbance, states(sys)); filter(isdisturbance, parameters(sys))]
end

## Tunable =====================================================================
struct VariableTunable end
Symbolics.option_to_metadata_type(::Val{:tunable}) = VariableTunable

istunable(x::Num, args...) = istunable(Symbolics.unwrap(x), args...)

"""
    istunable(x, default = false)

Determine whether or not symbolic variable `x` is marked as a tunable for an automatic tuning algorithm.

`default` indicates whether variables without `tunable` metadata are to be considered tunable or not.

Create a tunable parameter by
```
@parameters u [tunable=true]
```
See also [`tunable_parameters`](@ref), [`getbounds`](@ref)
"""
function istunable(x, default = false)
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

Get the probability distribution associated with symbolc variable `x`. If no distribution
is associated with `x`, `nothing` is returned.
Create parameters with associated distributions like this
```julia
using Distributions
d = Normal(0, 1)
@parameters u [dist=d]
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

Determine whether or not symbolic variable `x` has a probability distribution associated with it.
"""
function hasdist(x)
    b = getdist(x)
    b !== nothing
end

## System interface

"""
    tunable_parameters(sys, p = parameters(sys); default=false)

Get all parameters of `sys` that are marked as `tunable`.

Keyword argument `default` indicates whether variables without `tunable` metadata are to be considered tunable or not.

Create a tunable parameter by
```
@parameters u [tunable=true]
```
See also [`getbounds`](@ref), [`istunable`](@ref)
"""
function tunable_parameters(sys, p = parameters(sys); default = false)
    filter(x -> istunable(x, default), p)
end

"""
    getbounds(sys::ModelingToolkit.AbstractSystem, p = parameters(sys))

Returns a dict with pairs `p => (lb, ub)` mapping parameters of `sys` to lower and upper bounds.
Create parameters with bounds like this
```
@parameters p [bounds=(-1, 1)]
```

To obtain state bounds, call `getbounds(sys, states(sys))`
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

## binary variables =================================================================
struct VariableBinary end
Symbolics.option_to_metadata_type(::Val{:binary}) = VariableBinary

isbinaryvar(x::Num) = isbinaryvar(Symbolics.unwrap(x))

"""
    isbinaryvar(x)

Determine if a variable is binary.
"""
function isbinaryvar(x)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    return Symbolics.getmetadata(x, VariableBinary, false)
end

## integer variables =================================================================
struct VariableInteger end
Symbolics.option_to_metadata_type(::Val{:integer}) = VariableInteger

isintegervar(x::Num) = isintegervar(Symbolics.unwrap(x))

"""
    isintegervar(x)

Determine if a variable is integer.
"""
function isintegervar(x)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    return Symbolics.getmetadata(x, VariableInteger, false)
end
