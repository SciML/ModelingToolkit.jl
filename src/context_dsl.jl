import SymbolicUtils: symtype, term
struct Parameter{T} end

isparameter(x) = false
isparameter(::Sym{<:Parameter}) = true
isparameter(::Sym{<:FnType{<:Any, <:Parameter}}) = true

SymbolicUtils.@number_methods(Sym{Parameter{Real}},
                              term(f, a),
                              term(f, a, b), skipbasics)

SymbolicUtils.symtype(s::Symbolic{Parameter{T}}) where T = T
SymbolicUtils.similarterm(t::Term{<:Parameter}, f, args) = Term(f, args)

Base.convert(::Type{Num}, x::Symbolic{Parameter{T}}) where {T<:Number} = Num(x)

"""
$(SIGNATURES)

Define one or more known variables.
"""
macro parameters(xs...)
    esc(_parse_vars(:parameters, Parameter{Real}, xs))
end
