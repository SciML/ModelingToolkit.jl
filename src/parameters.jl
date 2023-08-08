import SymbolicUtils: symtype, term, hasmetadata, issym
@enum VariableType VARIABLE PARAMETER BROWNIAN
struct MTKVariableTypeCtx end

getvariabletype(x, def = VARIABLE) = getmetadata(unwrap(x), MTKVariableTypeCtx, def)

function isparameter(x)
    x = unwrap(x)

    if x isa Symbolic && (varT = getvariabletype(x, nothing)) !== nothing
        return varT === PARAMETER
        #TODO: Delete this branch
    elseif x isa Symbolic && Symbolics.getparent(x, false) !== false
        p = Symbolics.getparent(x)
        isparameter(p) ||
            (hasmetadata(p, Symbolics.VariableSource) &&
             getmetadata(p, Symbolics.VariableSource)[1] == :parameters)
    elseif istree(x) && operation(x) isa Symbolic
        varT === PARAMETER || isparameter(operation(x))
    elseif istree(x) && operation(x) == (getindex)
        isparameter(arguments(x)[1])
    elseif x isa Symbolic
        varT === PARAMETER
    else
        false
    end
end

"""
    toparam(s)

Maps the variable to a parameter.
"""
function toparam(s)
    if s isa Symbolics.Arr
        Symbolics.wrap(toparam(Symbolics.unwrap(s)))
    elseif s isa AbstractArray
        map(toparam, s)
    else
        setmetadata(s, MTKVariableTypeCtx, PARAMETER)
    end
end
toparam(s::Num) = wrap(toparam(value(s)))

"""
    tovar(s)

Maps the variable to a state.
"""
tovar(s::Symbolic) = setmetadata(s, MTKVariableTypeCtx, VARIABLE)
tovar(s::Num) = Num(tovar(value(s)))

"""
$(SIGNATURES)

Define one or more known parameters.
"""
macro parameters(xs...)
    Symbolics._parse_vars(:parameters,
        Real,
        xs,
        toparam) |> esc
end

struct Parameter{T <: Real}
    data::Vector{T}
    ref::T
    circular_buffer::Bool
end

Parameter(data::Vector{T}, ref::T) where {T <: Real} = Parameter(data, ref, true)
Parameter(x::Parameter) = x
function Parameter(x::T; tofloat = true) where {T <: Real}
    if tofloat
        x = float(x)
        P = typeof(x)
    else
        P = T
    end

    return Parameter(P[], x)
end

function Base.isequal(x::Parameter, y::Parameter)
    b0 = length(x.data) == length(y.data)
    if b0
        b1 = all(x.data .== y.data)
        b2 = x.ref == y.ref
        return b1 & b2
    else
        return false
    end
end

Base.:*(x::Number, y::Parameter) = x * y.ref
Base.:*(y::Parameter, x::Number) = Base.:*(x, y)
Base.:*(x::Parameter, y::Parameter) = x.ref * y.ref

Base.:/(x::Number, y::Parameter) = x / y.ref
Base.:/(y::Parameter, x::Number) = y.ref / x
Base.:/(x::Parameter, y::Parameter) = x.ref / y.ref

Base.:+(x::Number, y::Parameter) = x + y.ref
Base.:+(y::Parameter, x::Number) = Base.:+(x, y)
Base.:+(x::Parameter, y::Parameter) = x.ref + y.ref

Base.:-(y::Parameter) = -y.ref
Base.:-(x::Number, y::Parameter) = x - y.ref
Base.:-(y::Parameter, x::Number) = y.ref - x
Base.:-(x::Parameter, y::Parameter) = x.ref - y.ref

Base.:^(x::Number, y::Parameter) = Base.:^(x, y.ref)
Base.:^(y::Parameter, x::Number) = Base.:^(y.ref, x)
Base.:^(x::Parameter, y::Parameter) = Base.:^(x.ref, y.ref)

Base.isless(x::Parameter, y::Number) = Base.isless(x.ref, y)
Base.isless(y::Number, x::Parameter) = Base.isless(y, x.ref)

Base.copy(x::Parameter{T}) where {T} = Parameter{T}(copy(x.data), x.ref)

IfElse.ifelse(c::Bool, x::Parameter, y::Parameter) = ifelse(c, x.ref, y.ref)
IfElse.ifelse(c::Bool, x::Parameter, y::Number) = ifelse(c, x.ref, y)
IfElse.ifelse(c::Bool, x::Number, y::Parameter) = ifelse(c, x, y.ref)
Base.max(x::Number, y::Parameter) = max(x, y.ref)
Base.max(x::Parameter, y::Number) = max(x.ref, y)
Base.max(x::Parameter, y::Parameter) = max(x.ref, y.ref)

Base.min(x::Number, y::Parameter) = min(x, y.ref)
Base.min(x::Parameter, y::Number) = min(x.ref, y)
Base.min(x::Parameter, y::Parameter) = min(x.ref, y.ref)

function Base.show(io::IO, m::MIME"text/plain", p::Parameter)
    if !isempty(p.data)
        print(io, p.data)
    else
        print(io, p.ref)
    end
end

Base.convert(::Type{T}, x::Parameter{T}) where {T <: Real} = x.ref
function Base.convert(::Type{<:Parameter{T}}, x::Number) where {T <: Real}
    Parameter{T}(T[], x, true)
end
