"""
    $TYPEDEF

Symbolic metadata key for storing the unit associated with a symbolic variable.
"""
struct VariableUnit end
"""
    $TYPEDEF

Symbolic metadata key for storing the type of connector that a variable is.
"""
struct VariableConnectType end
"""
    $TYPEDEF

Symbolic metadata key for storing whether a symbolic variable is an input or not.
"""
struct VariableInput end
"""
    $TYPEDEF

Symbolic metadata key for storing whether a symbolic variable is an output or not.
"""
struct VariableOutput end
"""
  $TYPEDEF

Symbolic metadata key for storing whether a symbolic variable is irreducible or not.
"""
struct VariableIrreducible end
"""
    $TYPEDEF

Symbolic metadata key for storing the priority a variable has for being selected during
state selection.
"""
struct VariableStatePriority end
"""
    $TYPEDEF

Symbolic metadata key for storing miscellaneous information about a symbolic variable.
"""
struct VariableMisc end
# Metadata for renamed shift variables xₜ₋₁
struct VariableUnshifted end
struct VariableShift end
Symbolics.option_to_metadata_type(::Val{:unit}) = VariableUnit
Symbolics.option_to_metadata_type(::Val{:connect}) = VariableConnectType
Symbolics.option_to_metadata_type(::Val{:input}) = VariableInput
Symbolics.option_to_metadata_type(::Val{:output}) = VariableOutput
Symbolics.option_to_metadata_type(::Val{:irreducible}) = VariableIrreducible
Symbolics.option_to_metadata_type(::Val{:state_priority}) = VariableStatePriority
Symbolics.option_to_metadata_type(::Val{:misc}) = VariableMisc
Symbolics.option_to_metadata_type(::Val{:unshifted}) = VariableUnshifted
Symbolics.option_to_metadata_type(::Val{:shift}) = VariableShift

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
    variable_source,
    name = Symbolics.getmetadata(
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
"""
    $(TYPEDEF)

Flag which is meant to be passed to the `connect` metadata of a variable to affect how it
behaves when the connector it is in is part of a `connect` equation. `Equality` is the
default value and such variables when connected are made equal. For example, electric
potential is equated at a junction.

For more information, refer to the [Connection semantics](@ref connect_semantics) section
of the docs.

See also: [`connect`](@ref), [`@connector`](@ref), [`Flow`](@ref),
[`Stream`](@ref).
"""
struct Equality <: AbstractConnectType end # Equality connection
"""
    $(TYPEDEF)

Flag which is meant to be passed to the `connect` metadata of a variable to affect how it
behaves when the connector it is in is part of a `connect` equation. `Flow` denotes that
the sum of marked variable in all connectors in the connection set must sum to zero. For
example, electric current sums to zero at a junction (assuming appropriate signs are used
for current flowing in and out of the function).

For more information, refer to the [Connection semantics](@ref connect_semantics) section
of the docs.

See also: [`connect`](@ref), [`@connector`](@ref), [`Equality`](@ref),
[`Stream`](@ref).
"""
struct Flow <: AbstractConnectType end     # sum to 0
"""
    $(TYPEDEF)

Flag which is meant to be passed to the `connect` metadata of a variable to affect how it
behaves when the connector it is in is part of a `connect` equation. `Stream` denotes that
the variable is part of a special stream connector.

For more information, refer to the [Connection semantics](@ref connect_semantics) section
of the docs.

See also: [`connect`](@ref), [`@connector`](@ref), [`Equality`](@ref),
[`Flow`](@ref).
"""
struct Stream <: AbstractConnectType end   # special stream connector

"""
    getconnect(x)

Get the connect type of x. See also [`hasconnect`](@ref).
"""
getconnect(x::Num) = getconnect(unwrap(x))
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

"""
    $(TYPEDSIGNATURES)

Set the `input` metadata of variable `x` to `v`.
"""
setinput(x, v::Bool) = setmetadata(x, VariableInput, v)
"""
    $(TYPEDSIGNATURES)

Set the `output` metadata of variable `x` to `v`.
"""
setoutput(x, v::Bool) = setmetadata(x, VariableOutput, v)
setio(x, i::Bool, o::Bool) = setoutput(setinput(x, i), o)

"""
    $(TYPEDSIGNATURES)

Check if variable `x` is marked as an input.
"""
isinput(x) = isvarkind(VariableInput, x)
"""
    $(TYPEDSIGNATURES)

Check if variable `x` is marked as an output.
"""
isoutput(x) = isvarkind(VariableOutput, x)

# Before the solvability check, we already have handled IO variables, so
# irreducibility is independent from IO.
"""
    $(TYPEDSIGNATURES)

Check if `x` is marked as irreducible. This prevents it from being eliminated as an
observed variable in `mtkcompile`.
"""
isirreducible(x) = isvarkind(VariableIrreducible, x)
setirreducible(x, v::Bool) = setmetadata(x, VariableIrreducible, v)
state_priority(x::Union{Num, Symbolics.Arr}) = state_priority(unwrap(x))
"""
    $(TYPEDSIGNATURES)

Return the `state_priority` metadata of variable `x`. This influences its priority to be
chosen as a state in `mtkcompile`.
"""
state_priority(x) = convert(Float64, getmetadata(x, VariableStatePriority, 0.0))::Float64

normalize_to_differential(x) = x

function default_toterm(x)
    if iscall(x) && (op = operation(x)) isa Operator
        if !(op isa Differential)
            if op isa Shift && op.steps < 0
                return shift2term(x)
            end
            x = normalize_to_differential(op)(arguments(x)...)
        end
        Symbolics.diff2term(x)
    else
        x
    end
end

## Bounds ======================================================================
"""
    $TYPEDEF

Symbolic metadata key for specifying the bounds of a symbolic variable.
"""
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

function setbounds(x::Num, bounds)
    (lb, ub) = bounds
    setmetadata(x, VariableBounds, (lb, ub))
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

setdisturbance(x, v) = setmetadata(x, VariableDisturbance, v)

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
    tunable_parameters(sys, p = parameters(sys; initial_parameters = true); default=true)

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
function tunable_parameters(
        sys, p = parameters(sys; initial_parameters = true); default = true)
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
"""
    $TYPEDEF

Symbolic metadata key for storing the description of a symbolic variable.
"""
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

"""
    $(TYPEDSIGNATURES)

Check if variable `x` has a non-empty attached description.
"""
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
macro brownians(xs...)
    all(
        x -> x isa Symbol || Meta.isexpr(x, :call) && x.args[1] == :$ || Meta.isexpr(x, :$),
        xs) ||
        error("@brownians only takes scalar expressions!")
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
getmisc(x::Num) = getmisc(unwrap(x))
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
getunit(x::Num) = getunit(unwrap(x))
getunit(x::Symbolic) = Symbolics.getmetadata(x, VariableUnit, nothing)
"""
    hasunit(x)

Check if the variable `x` has a unit.
"""
hasunit(x) = getunit(x) !== nothing

getunshifted(x::Num) = getunshifted(unwrap(x))
getunshifted(x::Symbolic) = Symbolics.getmetadata(x, VariableUnshifted, nothing)

getshift(x::Num) = getshift(unwrap(x))
getshift(x::Symbolic) = Symbolics.getmetadata(x, VariableShift, 0)

###################
### Evaluate at ###
###################
"""
    EvalAt(t)

An operator that evaluates time-dependent variables at a specific absolute time point `t`.

# Fields
- `t::Union{Symbolic, Number}`: The absolute time at which to evaluate the variable.

# Description
`EvalAt` is used to evaluate time-dependent variables at a specific time point. This is particularly 
useful in optimization problems where you need to specify constraints or costs at particular moments 
in time, or delay differential equations for setting a delay time.

The operator works by replacing the time argument of time-dependent variables with the specified 
time `t`. For variables that don't depend on time, `EvalAt` returns them unchanged.

# Behavior
- For time-dependent variables like `x(t)`, `EvalAt(τ)(x)` returns `x(τ)` 
- For time-independent parameters, `EvalAt` returns them unchanged
- For derivatives, `EvalAt` evaluates the derivative at the specified time
- For arrays of variables, `EvalAt` is applied element-wise


# Examples
```julia
using ModelingToolkit

@variables t x(t) y(t)
@parameters p

# Evaluate x at time t=1.0
EvalAt(1.0)(x)  # Returns x(1.0)

# Works with parameters (returns unchanged)
EvalAt(1.0)(p)  # Returns p

# Works with derivatives
D = Differential(t)
EvalAt(1.0)(D(x))  # Returns D(x) evaluated at t=1.0

# Use in optimization constraints
@optimization_model model begin
    @constraints begin
        EvalAt(0.5)(x) ~ 2.0  # x must equal 2.0 at t=0.5
    end
end
```

# Errors
- Throws an error when applied to variables with more than one argument (e.g., `z(u, t)`)

See also: [`Differential`](@ref)
"""
struct EvalAt <: Symbolics.Operator
    t::Union{Symbolic, Number}
end

function (A::EvalAt)(x::Symbolic)
    if symbolic_type(x) == NotSymbolic() || !iscall(x)
        if x isa Symbolics.CallWithMetadata
            return x(A.t)
        else
            return x
        end
    end

    if iscall(x) && operation(x) == getindex
        arr = arguments(x)[1]
        term(getindex, A(arr), arguments(x)[2:end]...)
    elseif operation(x) isa Differential
        x = default_toterm(x)
        A(x)
    else
        length(arguments(x)) !== 1 &&
            error("Variable $x has too many arguments. EvalAt can only be applied to one-argument variables.")
        (symbolic_type(only(arguments(x))) !== ScalarSymbolic()) && return x
        return operation(x)(A.t)
    end
end

function (A::EvalAt)(x::Union{Num, Symbolics.Arr})
    wrap(A(unwrap(x)))
end
SymbolicUtils.isbinop(::EvalAt) = false

Base.nameof(::EvalAt) = :EvalAt
Base.show(io::IO, A::EvalAt) = print(io, "EvalAt(", A.t, ")")
Base.:(==)(A1::EvalAt, A2::EvalAt) = isequal(A1.t, A2.t)
Base.hash(A::EvalAt, u::UInt) = hash(A.t, u)
