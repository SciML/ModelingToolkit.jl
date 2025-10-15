""""""
struct VariableUnit end
""""""
struct VariableConnectType end
""""""
struct VariableInput end
""""""
struct VariableOutput end
""""""
struct VariableIrreducible end
""""""
struct VariableStatePriority end
""""""
struct VariableMisc end
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
""""""
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
abstract type AbstractConnectType end
""""""
struct Equality <: AbstractConnectType end
""""""
struct Flow <: AbstractConnectType end
""""""
struct Stream <: AbstractConnectType end
""""""
getconnect(x::Num) = getconnect(unwrap(x))
getconnect(x::SymbolicT) = Symbolics.getmetadata(x, VariableConnectType, nothing)
""""""
hasconnect(x) = getconnect(x) !== nothing
function setconnect(x, t::Type{T}) where {T <: AbstractConnectType}
    setmetadata(x, VariableConnectType, t)
end
isvarkind(m, x, def = false) = safe_getmetadata(m, x, def)
safe_getmetadata(m, x::Union{Num, Symbolics.Arr}, def) = safe_getmetadata(m, value(x), def)
function safe_getmetadata(m::DataType, x::SymbolicT, default)
    hasmetadata(x, m) && return getmetadata(x, m)
    iscall(x) && operation(x) === getindex && return safe_getmetadata(m, arguments(x)[1], default)
    return default
end
""""""
setinput(x, v::Bool) = setmetadata(x, VariableInput, v)
""""""
setoutput(x, v::Bool) = setmetadata(x, VariableOutput, v)
setio(x, i::Bool, o::Bool) = setoutput(setinput(x, i), o)
""""""
isinput(x) = isvarkind(VariableInput, x)::Bool
""""""
isoutput(x) = isvarkind(VariableOutput, x)::Bool
""""""
isirreducible(x) = isvarkind(VariableIrreducible, x)::Bool
setirreducible(x, v::Bool) = setmetadata(x, VariableIrreducible, v)
state_priority(x::Union{Num, Symbolics.Arr}) = state_priority(unwrap(x))
""""""
state_priority(x) = convert(Float64, getmetadata(x, VariableStatePriority, 0.0))::Float64
function normalize_to_differential(@nospecialize(op))
    if op isa Shift && op.t isa SymbolicT
        return Differential(op.t) ^ op.steps
    else
        return op
    end
end
default_toterm(x) = x
function default_toterm(x::SymbolicT)
    Moshi.Match.@match x begin
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
""""""
struct VariableBounds end
Symbolics.option_to_metadata_type(::Val{:bounds}) = VariableBounds
""""""
getbounds(x::Union{Num, Symbolics.Arr}) = getbounds(unwrap(x))
function getbounds(x::SymbolicT)
    if operation(p) === getindex
        p = arguments(p)[1]
        bounds = Symbolics.getmetadata(x, VariableBounds, (-Inf, Inf))
        if symbolic_type(x) == ArraySymbolic() && symbolic_has_known_size(x)
            bounds = map(bounds) do b
                b isa AbstractArray && return b
                return fill(b, size(x))
            end
        end
    else
        bounds = @something Symbolics.getmetadata(x, VariableBounds, nothing) getbounds(p)
        idxs = arguments(x)[2:end]
        bounds = map(bounds) do b
            if b isa AbstractArray
                if symbolic_has_known_size(p) && size(p) != size(b)
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
""""""
function hasbounds(x)
    b = getbounds(x)::NTuple{2}
    any(isfinite.(b[1]) .|| isfinite.(b[2]))::Bool
end
function setbounds(x::Num, bounds)
    (lb, ub) = bounds
    setmetadata(x, VariableBounds, (lb, ub))
end
struct VariableDisturbance end
Symbolics.option_to_metadata_type(::Val{:disturbance}) = VariableDisturbance
isdisturbance(x::Num) = isdisturbance(Symbolics.unwrap(x))
""""""
function isdisturbance(x)
    isvarkind(VariableDisturbance, x)::Bool
end
setdisturbance(x, v) = setmetadata(x, VariableDisturbance, v)
function disturbances(sys)
    [filter(isdisturbance, unknowns(sys)); filter(isdisturbance, parameters(sys))]
end
struct VariableTunable end
Symbolics.option_to_metadata_type(::Val{:tunable}) = VariableTunable
istunable(x::Num, args...) = istunable(Symbolics.unwrap(x), args...)
""""""
function istunable(x, default = true)
    isvarkind(VariableTunable, x, default)::Bool
end
struct VariableDistribution end
Symbolics.option_to_metadata_type(::Val{:dist}) = VariableDistribution
getdist(x::Num) = getdist(Symbolics.unwrap(x))
""""""
function getdist(x)
    safe_getmetadata(VariableDistribution, x, nothing)
end
""""""
function hasdist(x)
    b = getdist(x)
    b !== nothing
end
""""""
function tunable_parameters(
        sys, p = parameters(sys; initial_parameters = true); default = true)
    filter(x -> istunable(x, default), p)
end
""""""
function getbounds(sys::ModelingToolkit.AbstractSystem, p = parameters(sys))
    Dict(p .=> getbounds.(p))
end
""""""
function getbounds(p::AbstractVector)
    bounds = getbounds.(p)
    lb = first.(bounds)
    ub = last.(bounds)
    (; lb, ub)
end
""""""
struct VariableDescription end
Symbolics.option_to_metadata_type(::Val{:description}) = VariableDescription
getdescription(x::Num) = getdescription(Symbolics.unwrap(x))
getdescription(x::Symbolics.Arr) = getdescription(Symbolics.unwrap(x))
""""""
function getdescription(x)
    safe_getmetadata(VariableDescription, x, "")
end
""""""
function hasdescription(x)
    getdescription(x) != ""
end
""""""
tobrownian(s::SymbolicT) = setmetadata(s, MTKVariableTypeCtx, BROWNIAN)
tobrownian(s::Num) = Num(tobrownian(value(s)))
isbrownian(s) = getvariabletype(s) === BROWNIAN
""""""
macro brownians(xs...)
    all(
        x -> x isa Symbol || Meta.isexpr(x, :call) && x.args[1] == :$ || Meta.isexpr(x, :$),
        xs) ||
        error("@brownians only takes scalar expressions!")
    Symbolics.parse_vars(:brownian,
        Real,
        xs,
        tobrownian)
end
struct VariableGuess end
Symbolics.option_to_metadata_type(::Val{:guess}) = VariableGuess
getguess(x::Union{Num, Symbolics.Arr}) = getguess(Symbolics.unwrap(x))
""""""
function getguess(x)
    Symbolics.getmetadata(x, VariableGuess, nothing)
end
""""""
function setguess(x, v)
    Symbolics.setmetadata(x, VariableGuess, v)
end
""""""
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
""""""
getmisc(x::Num) = getmisc(unwrap(x))
getmisc(x::SymbolicT) = Symbolics.getmetadata(x, VariableMisc, nothing)
""""""
hasmisc(x) = getmisc(x) !== nothing
setmisc(x, miscdata) = setmetadata(x, VariableMisc, miscdata)
""""""
getunit(x::Num) = getunit(unwrap(x))
getunit(x::SymbolicT) = Symbolics.getmetadata(x, VariableUnit, nothing)
""""""
hasunit(x) = getunit(x) !== nothing
getunshifted(x::Num) = getunshifted(unwrap(x))
getunshifted(x::SymbolicT) = Symbolics.getmetadata(x, VariableUnshifted, nothing)::Union{SymbolicT, Nothing}
getshift(x::Num) = getshift(unwrap(x))
getshift(x::SymbolicT) = Symbolics.getmetadata(x, VariableShift, 0)::Int
""""""
struct EvalAt <: Symbolics.Operator
    t::Union{SymbolicT, Number}
end
function (A::EvalAt)(x::SymbolicT)
    if symbolic_type(x) == NotSymbolic() || !iscall(x)
        if x isa CallAndWrap
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
