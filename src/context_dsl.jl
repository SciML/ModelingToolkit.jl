struct Parameter{T} end
isparameter(::Variable) = false
isparameter(::Variable{<:Parameter}) = true

"""
$(SIGNATURES)

Define one or more known variables.
"""
macro parameters(xs...)
    esc(_parse_vars(:parameters, Parameter{Real}, xs))
end
