module SciCompDSL

abstract type AbstractSystem end

# Parameterize by T so that way it can be Vector{Expression} which is defined after
struct Operation{T}
    op::Function
    args::Vector{T}
end
include("variables.jl")
const Expression = Union{Variable,Operation}

export Operation, Expression

include("operators.jl")
include("equations.jl")
include("function_registration.jl")

end # module
