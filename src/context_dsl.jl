import SymbolicUtils: symtype
struct Parameter{T} end

isparameter(::Sym) = false
isparameter(::Sym{<:Parameter}) = true
isparameter(::Sym{<:FnType{<:Any, <:Parameter}}) = true

SymbolicUtils.symtype(s::Sym{Parameter{T}}) where T = T
SymbolicUtils.symtype(s::Sym{Parameter{<:FnType{<:Any, T}}}) where T = T

"""
$(SIGNATURES)

Define one or more known variables.
"""
macro parameters(xs...)
    esc(_parse_vars(:parameters, Parameter{Number}, xs))
end
