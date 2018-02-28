struct Variable <: Real
    name::Symbol
    subtype::Symbol
    value
    value_type::DataType
end

Variable(name,subtype,value,value_type = typeof(value)) = Variable(name,subtype,value,value_type)
Variable(name,subtype) = Variable(name,subtype,nothing,Void)
Variable(name,args...) = Variable(name,:None,args...)
Parameter(name,args...) = Variable(name,:Parameter,args...)
Constant(value) = Variable(:None,:Constant,value,typeof(value))
DependentVariable(name,args...) = Variable(name,:DependentVariable,args...)
IndependentVariable(name,args...) = Variable(name,:IndependentVariable,args...)

export Variable,Parameter,Constant,DependentVariable,IndependentVariable
