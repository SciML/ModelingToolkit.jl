struct Variable <: AbstractVariable
    name::Symbol
end
struct Parameter <: AbstractVariable
    name::Symbol
end
struct Constant{T<:Number} <: AbstractVariable
    value::T
end
struct DependentVariable <: AbstractVariable
    name::Symbol
end
struct IndependentVariable <: AbstractVariable
    name::Symbol
end

export Variable,Parameter,Constant,DependentVariable,IndependentVariable
