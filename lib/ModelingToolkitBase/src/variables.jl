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
struct VariableTimeDomain end

Symbolics.option_to_metadata_type(::Val{:unit}) = VariableUnit
Symbolics.option_to_metadata_type(::Val{:connect}) = VariableConnectType
Symbolics.option_to_metadata_type(::Val{:input}) = VariableInput
Symbolics.option_to_metadata_type(::Val{:output}) = VariableOutput
Symbolics.option_to_metadata_type(::Val{:irreducible}) = VariableIrreducible
Symbolics.option_to_metadata_type(::Val{:state_priority}) = VariableStatePriority
Symbolics.option_to_metadata_type(::Val{:misc}) = VariableMisc
Symbolics.option_to_metadata_type(::Val{:unshifted}) = VariableUnshifted
Symbolics.option_to_metadata_type(::Val{:shift}) = VariableShift
Symbolics.option_to_metadata_type(::Val{:timedomain}) = VariableTimeDomain

"""
    dump_variable_metadata(var)

Return all the metadata associated with symbolic variable `var` as a `NamedTuple`.

```@example
using ModelingToolkitBase

@parameters p::Int [description = "My description", bounds = (0.5, 1.5)]
ModelingToolkitBase.dump_variable_metadata(p)
```
"""
function dump_variable_metadata(var)
    uvar = unwrap(var)
    variable_source,
        name = Symbolics.getmetadata(
        uvar, VariableSource, (:unknown, :unknown)
    )
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
        default,
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
getconnect(x::SymbolicT) = Symbolics.getmetadata(x, VariableConnectType, nothing)
"""
    hasconnect(x)

Determine whether variable `x` has a connect type. See also [`getconnect`](@ref).
"""
hasconnect(x) = getconnect(x) !== nothing
function setconnect(x, t::Type{T}) where {T <: AbstractConnectType}
    return setmetadata(x, VariableConnectType, t)
end

### Input, Output, Irreducible
isvarkind(m, x, def = false) = safe_getmetadata(m, x, def)
safe_getmetadata(m, x::Union{Num, Symbolics.Arr}, def) = safe_getmetadata(m, value(x), def)
function safe_getmetadata(m::DataType, x::SymbolicT, default)
    hasmetadata(x, m) && return getmetadata(x, m)
    iscall(x) && operation(x) === getindex && return safe_getmetadata(m, arguments(x)[1], default)
    return default
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
isinput(x) = isvarkind(VariableInput, x)::Bool
"""
    $(TYPEDSIGNATURES)

Check if variable `x` is marked as an output.
"""
isoutput(x) = isvarkind(VariableOutput, x)::Bool

# Before the solvability check, we already have handled IO variables, so
# irreducibility is independent from IO.
"""
    $(TYPEDSIGNATURES)

Check if `x` is marked as irreducible. This prevents it from being eliminated as an
observed variable in `mtkcompile`.
"""
isirreducible(x) = isvarkind(VariableIrreducible, x)::Bool
setirreducible(x, v::Bool) = setmetadata(x, VariableIrreducible, v)
state_priority(x::Union{Num, Symbolics.Arr}) = state_priority(unwrap(x))
"""
    $(TYPEDSIGNATURES)

Return the `state_priority` metadata of variable `x`. This influences its priority to be
chosen as a state in `mtkcompile`.
"""
state_priority(x) = convert(Float64, getmetadata(x, VariableStatePriority, 0.0))::Float64

function normalize_to_differential(@nospecialize(op))
    if op isa Shift && op.t isa SymbolicT
        return Differential(op.t)^op.steps
    else
        return op
    end
end

default_toterm(x) = x
function default_toterm(x::SymbolicT)
    return Moshi.Match.@match x begin
        BSImpl.Term(; f, args, shape, type, metadata) && if f isa Operator end => begin
            if f isa Shift && f.steps < 0
                return shift2term(x)
            elseif f isa Differential
                return Symbolics.diff2term(x)
            else
                newf = normalize_to_differential(f)
                f === newf && return x
                x = BSImpl.Term{VartypeT}(newf, args; type, shape, metadata)
                return Symbolics.diff2term(x)
            end
        end
        _ => return x
    end
end

"""
Rename a Shift variable with negative shift, Shift(t, k)(x(t)) to xₜ₋ₖ(t).
"""
function shift2term(var::SymbolicT)
    return Moshi.Match.@match var begin
        BSImpl.Term(f, args) && if f isa Shift end => begin
            op = f
            arg = args[1]
            Moshi.Match.@match arg begin
                BSImpl.Term(; f, args, type, shape, metadata) && if f === getindex end => begin
                    newargs = copy(parent(args))
                    newargs[1] = shift2term(op(newargs[1]))
                    unshifted_args = copy(newargs)
                    unshifted_args[1] = ModelingToolkitBase.getunshifted(newargs[1])::SymbolicT
                    unshifted = BSImpl.Term{VartypeT}(getindex, unshifted_args; type, shape, metadata)
                    if metadata === nothing
                        metadata = Base.ImmutableDict{DataType, Any}(VariableUnshifted, unshifted)
                    elseif metadata isa Base.ImmutableDict{DataType, Any}
                        metadata = Base.ImmutableDict(metadata, VariableUnshifted, unshifted)
                    end
                    return BSImpl.Term{VartypeT}(getindex, newargs; type, shape, metadata)
                end
                _ => nothing
            end
            unshifted = ModelingToolkitBase.getunshifted(arg)
            is_lowered = unshifted !== nothing
            backshift = op.steps + ModelingToolkitBase.getshift(arg)
            iszero(backshift) && return unshifted::SymbolicT
            io = IOBuffer()
            O = (is_lowered ? unshifted : arg)::SymbolicT
            write(io, getname(O))
            # Char(0x209c) = ₜ
            write(io, Char(0x209c))
            # Char(0x208b) = ₋ (subscripted minus)
            # Char(0x208a) = ₊ (subscripted plus)
            pm = backshift > 0 ? Char(0x208a) : Char(0x208b)
            write(io, pm)
            _backshift = backshift
            backshift = abs(backshift)
            N = ndigits(backshift)
            den = 10^(N - 1)
            for _ in 1:N
                # subscripted number, e.g. ₁
                write(io, Char(0x2080 + div(backshift, den) % 10))
                den = div(den, 10)
            end
            newname = Symbol(take!(io))
            newvar = Symbolics.rename(arg, newname)
            newvar = setmetadata(newvar, ModelingToolkitBase.VariableUnshifted, O)
            newvar = setmetadata(newvar, ModelingToolkitBase.VariableShift, _backshift)
            return newvar
        end
        _ => return var
    end
end

simplify_shifts(eq::Equation) = simplify_shifts(eq.lhs) ~ simplify_shifts(eq.rhs)

function _simplify_shifts(var::SymbolicT)
    return Moshi.Match.@match var begin
        BSImpl.Term(; f, args) && if f isa Shift && f.steps == 0 end => return args[1]
        BSImpl.Term(; f = op1, args) && if op1 isa Shift end => begin
            vv1 = args[1]
            Moshi.Match.@match vv1 begin
                BSImpl.Term(; f = op2, args = a2) && if op2 isa Shift end => begin
                    vv2 = a2[1]
                    s1 = op1.steps
                    s2 = op2.steps
                    t1 = op1.t
                    t2 = op2.t
                    return simplify_shifts(ModelingToolkitBase.Shift(t1 === nothing ? t2 : t1, s1 + s2)(vv2))
                end
                _ => return var
            end
        end
        _ => var
    end
end

"""
Simplify multiple shifts: Shift(t, k1)(Shift(t, k2)(x)) becomes Shift(t, k1+k2)(x).
"""
function simplify_shifts(var::SymbolicT)
    ModelingToolkitBase.hasshift(var) || return var
    return SU.Rewriters.Postwalk(_simplify_shifts)(var)::SymbolicT
end

distribute_shift(eq::Equation) = distribute_shift(eq.lhs) ~ distribute_shift(eq.rhs)
distribute_shift(var::Union{Num, Arr}) = distribute_shift(unwrap(var))
"""
Distribute a shift applied to a whole expression or equation. 
Shift(t, 1)(x + y) will become Shift(t, 1)(x) + Shift(t, 1)(y).
Only shifts variables whose independent variable is the same t that appears in the Shift (i.e. constants, time-independent parameters, etc. do not get shifted).
"""
function distribute_shift(var::SymbolicT)
    return Moshi.Match.@match var begin
        BSImpl.Term(; f, args) && if f isa Shift end => begin
            shiftexpr = _distribute_shift(args[1], f)
            return simplify_shifts(shiftexpr)
        end
        _ => return var
    end
end

"""
    $TYPEDSIGNATURES

Whether `distribute_shift` should distribute shifts into the given operation.
"""
distribute_shift_into_operator(_) = true

function _distribute_shift(expr::SymbolicT, shift)
    if iscall(expr)
        op = operation(expr)
        distribute_shift_into_operator(op)::Bool || return expr
        args = arguments(expr)

        if ModelingToolkitBase.isvariable(expr) && operation(expr) !== getindex &&
                !ModelingToolkitBase.iscalledparameter(expr)
            (length(args) == 1 && isequal(shift.t, only(args))) ? (return shift(expr)) :
                (return expr)
        elseif op isa Shift
            return shift(expr)
        else
            return maketerm(
                typeof(expr), operation(expr), Base.Fix2(_distribute_shift, shift).(args),
                unwrap(expr).metadata
            )
        end
    else
        return expr
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
getbounds(x::Union{Num, Symbolics.Arr}) = getbounds(unwrap(x))
function getbounds(x::SymbolicT)
    if hasmetadata(x, VariableBounds)
        return getmetadata(x, VariableBounds)
    end
    arrx, isarr = split_indexed_var(x)
    if !isarr || !hasmetadata(arrx, VariableBounds)
        if SU.is_array_shape(SU.shape(x))
            return (fill(-Inf, size(x)), fill(Inf, size(x)))
        end
        return (-Inf, Inf)
    end
    bounds = getmetadata(arrx, VariableBounds, nothing)::NTuple{2, Any}
    idxs = @views unwrap_const.(arguments(x)[2:end])
    return map(bounds) do b
        @assert !symbolic_has_known_size(arrx) || SU.shape(arrx) == SU.shape(b)
        return b[idxs...]
    end
end

"""
    hasbounds(x)

Determine whether symbolic variable `x` has bounds associated with it.
See also [`getbounds`](@ref).
"""
function hasbounds(x)
    b = getbounds(x)::NTuple{2, Any}
    return any(isfinite.(b[1]) .|| isfinite.(b[2]))::Bool
end

function setbounds(x::Num, bounds)
    (lb, ub) = bounds
    return setmetadata(x, VariableBounds, (lb, ub))
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
    return isvarkind(VariableDisturbance, x)::Bool
end

setdisturbance(x, v) = setmetadata(x, VariableDisturbance, v)

function disturbances(sys)
    return [filter(isdisturbance, unknowns(sys)); filter(isdisturbance, parameters(sys))]
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
    return isvarkind(VariableTunable, x, default)::Bool
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
    return safe_getmetadata(VariableDistribution, x, nothing)
end

"""
    hasdist(x)

Determine whether symbolic variable `x` has a probability distribution associated with it.
"""
function hasdist(x)
    b = getdist(x)
    return b !== nothing
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
        sys, p = parameters(sys; initial_parameters = true); default = true
    )
    return filter(x -> istunable(x, default), p)
end

"""
    getbounds(sys::ModelingToolkitBase.AbstractSystem, p = parameters(sys))

Returns a dict with pairs `p => (lb, ub)` mapping parameters of `sys` to lower and upper bounds.
Create parameters with bounds like this

```
@parameters p [bounds=(-1, 1)]
```

To obtain unknown variable bounds, call `getbounds(sys, unknowns(sys))`
"""
function getbounds(sys::ModelingToolkitBase.AbstractSystem, p = parameters(sys))
    return Dict(p .=> getbounds.(p))
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
    return (; lb, ub)
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
    return safe_getmetadata(VariableDescription, x, "")
end

"""
    $(TYPEDSIGNATURES)

Check if variable `x` has a non-empty attached description.
"""
function hasdescription(x)
    return getdescription(x) != ""
end

## Brownian
"""
    tobrownian(s::Sym)

Maps the brownianiable to an unknown.
"""
tobrownian(s::SymbolicT) = setmetadata(s, MTKVariableTypeCtx, BROWNIAN)
tobrownian(s::Num) = Num(tobrownian(value(s)))
isbrownian(s) = getvariabletype(s) === BROWNIAN

"""
$(SIGNATURES)

Define one or more Brownian variables.
"""
macro brownians(xs...)
    all(
        x -> x isa Symbol || Meta.isexpr(x, :call) && x.args[1] == :$ || Meta.isexpr(x, :$),
        xs
    ) ||
        error("@brownians only takes scalar expressions!")
    return Symbolics.parse_vars(
        :brownian,
        Real,
        xs,
        tobrownian
    )
end

## Poissonian ==================================================================
"""
    topoissonian(s::Sym, rate)

Maps the variable to a poissonian with the given rate expression stored in metadata.
"""
function topoissonian(s::SymbolicT, rate)
    s = setmetadata(s, MTKVariableTypeCtx, POISSONIAN)
    s = setmetadata(s, PoissonianRateCtx, rate)
    return s
end
topoissonian(s::Num, rate) = Num(topoissonian(value(s), rate))
ispoissonian(s) = getvariabletype(s) === POISSONIAN

"""
$(SIGNATURES)

Define one or more Poissonian variables with their rate expressions.

Each poissonian represents the differential of a Poisson counting process with
the specified rate. Unlike `@brownians`, a rate expression is required for each
poissonian.

# Examples
```julia
@poissonians dN(λ)              # Single declaration with constant rate
@poissonians dN₁(λ₁) dN₂(λ₂)    # Multiple inline declarations
@poissonians begin              # Block syntax
    dN₁(λ₁)
    dN₂(β*S*I)
end
```
"""
macro poissonians(exprs...)
    return esc(_poissonians(exprs...))
end

function _poissonians(exprs...)
    # Handle block syntax: @poissonians begin ... end
    if length(exprs) == 1 && exprs[1] isa Expr && exprs[1].head == :block
        # Filter out LineNumberNodes and process each expression
        inner_exprs = filter(x -> !(x isa LineNumberNode), exprs[1].args)
        return _poissonians(inner_exprs...)
    end

    assignments = Expr[]
    names = Symbol[]
    for expr in exprs
        # Must be a call expression: dN(rate)
        if !(expr isa Expr && expr.head == :call)
            error("@poissonians requires a rate expression: use @poissonians dN(rate)")
        end

        name = expr.args[1]
        if length(expr.args) < 2
            error("@poissonians requires a rate expression: use @poissonians $name(rate)")
        end
        rate = expr.args[2]

        if !(name isa Symbol)
            error("@poissonians variable name must be a symbol, got: $name")
        end

        push!(names, name)
        # Create the symbolic variable using Symbolics.variable and set poissonian metadata
        # Symbolics.variable creates a proper Sym{VartypeT} with VariableSource metadata
        push!(assignments, quote
            $name = $topoissonian($(Symbolics.variable)($(QuoteNode(name))), $rate)
        end)
    end

    # Return the variables as a tuple (or single if only one)
    if length(names) == 1
        return Expr(:block, assignments..., names[1])
    else
        return Expr(:block, assignments..., Expr(:tuple, names...))
    end
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
    return Symbolics.getmetadata_maybe_indexed(x, VariableGuess, nothing)
end

"""
    setguess(x, v)

Set the guess for the initial value associated with symbolic variable `x` to `v`.
See also [`hasguess`](@ref).
"""
function setguess(x, v)
    return Symbolics.setmetadata(x, VariableGuess, v)
end

"""
    hasguess(x)

Determine whether symbolic variable `x` has a guess associated with it.
See also [`getguess`](@ref).
"""
function hasguess(x)
    return getguess(x) !== nothing
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
getmisc(x::SymbolicT) = Symbolics.getmetadata(x, VariableMisc, nothing)
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
getunit(x::SymbolicT) = Symbolics.getmetadata(x, VariableUnit, nothing)
"""
    hasunit(x)

Check if the variable `x` has a unit.
"""
hasunit(x) = getunit(x) !== nothing

getunshifted(x::Num) = getunshifted(unwrap(x))
getunshifted(x::SymbolicT) = Symbolics.getmetadata(x, VariableUnshifted, nothing)::Union{SymbolicT, Nothing}

getshift(x::Num) = getshift(unwrap(x))
getshift(x::SymbolicT) = Symbolics.getmetadata(x, VariableShift, 0)::Int

###################
### Evaluate at ###
###################
"""
    EvalAt(t)

An operator that evaluates time-dependent variables at a specific absolute time point `t`.

# Fields
- `t::Union{SymbolicT, Number}`: The absolute time at which to evaluate the variable.

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
using ModelingToolkitBase

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
    t::Union{SymbolicT, Number}
end

function (A::EvalAt)(x::SymbolicT)
    iscall(x) || return x
    if operation(x) === getindex
        arr = arguments(x)[1]
        term(getindex, A(arr), arguments(x)[2:end]...)
    elseif operation(x) isa Differential
        x = default_toterm(x)
        A(x)
    else
        length(arguments(x)) !== 1 &&
            error("Variable $x has too many arguments. EvalAt can only be applied to one-argument variables.")
        SU.isconst(only(arguments(x))) && return x
        return operation(x)(A.t)
    end
end

function (A::EvalAt)(x::Union{Num, Symbolics.Arr})
    return wrap(A(unwrap(x)))
end

function (A::EvalAt)(x::CallAndWrap)
    return x(A.t)
end

SymbolicUtils.isbinop(::EvalAt) = false

Base.nameof(::EvalAt) = :EvalAt
Base.show(io::IO, A::EvalAt) = print(io, "EvalAt(", A.t, ")")
Base.:(==)(A1::EvalAt, A2::EvalAt) = isequal(A1.t, A2.t)
Base.hash(A::EvalAt, u::UInt) = hash(A.t, u)
