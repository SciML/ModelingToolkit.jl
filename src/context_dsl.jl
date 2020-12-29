import SymbolicUtils: symtype
struct Parameter{T} end

isparameter(x) = false
isparameter(::Sym{<:Parameter}) = true
isparameter(::Sym{<:FnType{<:Any, <:Parameter}}) = true

SymbolicUtils.symtype(s::Symbolic{Parameter{T}}) where T = T
SymbolicUtils.similarterm(t::Term{T}, f, args) where {T<:Parameter} = Term{T}(f, args)

Base.convert(::Type{Num}, x::Symbolic{Parameter{T}}) where {T<:Number} = Num(x)

"""
$(SIGNATURES)

Define one or more known variables.
"""
macro parameters(xs...)
    esc(_parse_vars(:parameters, Parameter{Real}, xs))
end
