module SciCompDSL

abstract type AbstractVariable <: Real end
abstract type AbstractStatement end
abstract type AbstractSystem end
const Expression = Union{AbstractVariable,AbstractStatement}

struct Operation <: AbstractStatement
    op::Function
    args::Vector{Expression}
end
struct Statement <: AbstractStatement
    lhs::Expression
    rhs::Expression
end

export AbstractVariable, AbstractStatement, AbstractSystem, Expression,
       Operation, Statement

include("variables.jl")
include("derivatives.jl")
include("equations.jl")
include("function_registration.jl")

end # module
