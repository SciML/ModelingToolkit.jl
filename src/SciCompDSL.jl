module SciCompDSL

using DiffEqBase
import MacroTools: splitdef, combinedef
import IterTools: product

abstract type AbstractOperation end
abstract type AbstractSystem end
abstract type AbstractOperator end

include("variables.jl")

const Expression = Union{Variable,AbstractOperation,AbstractOperator}

include("operations.jl")
include("operators.jl")
include("systems.jl")
include("function_registration.jl")
include("simplify.jl")
include("utils.jl")

export Operation, Expression, AbstractOperator
export @register
end # module
