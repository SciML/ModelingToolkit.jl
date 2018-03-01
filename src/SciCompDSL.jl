module SciCompDSL

using DiffEqBase
import MacroTools: splitdef, combinedef
import IterTools: product

abstract type AbstractSystem end
abstract type AbstractOperator end

include("operations.jl")
include("variables.jl")

const Expression = Union{Variable,Operation,AbstractOperator}

export Operation, Expression, AbstractOperator

include("operators.jl")
include("systems.jl")
include("function_registration.jl")
include("simplify.jl")
include("utils.jl")
export @register
end # module
