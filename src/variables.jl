struct Variable <: Real
    name::Symbol
    subtype::Symbol
    value
    value_type::DataType
end

Variable(name,subtype) = Variable(name,subtype,nothing,Void)
Variable(name) = Variable(name,:None)
Parameter(name) = Variable(name,:Parameter)
Constant(value) = Variable(:None,:Constant,value,typeof(value))
DependentVariable(name) = Variable(name,:DependentVariable)
IndependentVariable(name) = Variable(name,:IndependentVariable)

export Variable,Parameter,Constant,DependentVariable,IndependentVariable
