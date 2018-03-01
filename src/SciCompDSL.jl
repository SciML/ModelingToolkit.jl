module SciCompDSL

abstract type AbstractSystem end
abstract type AbstractOperator end

include("operations.jl")
include("variables.jl")

const Expression = Union{Variable,Operation,AbstractOperator}

export Operation, Expression, AbstractOperator

include("operators.jl")
include("equations.jl")
include("function_registration.jl")
export @register
end # module
