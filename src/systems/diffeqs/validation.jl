Base.:*(x::Expression,y::Unitful.AbstractQuantity) = x * ModelingToolkit.Constant(y)

instantiate(x::ModelingToolkit.Constant) = x.value
instantiate(x::ModelingToolkit.Variable{Real}) = 1.0
instantiate(x::ModelingToolkit.Variable) = oneunit(1*ModelingToolkit.vartype(x))
function instantiate(x::ModelingToolkit.Operation)
    if x.op isa Variable
        return instantiate(x.op)
    elseif x.op isa Differential
        instantiate(x.args[1])/instantiate(x.args[1].args[1])
    else
        x.op(instantiate.(x.args)...)
    end
end

function validate(eq::ModelingToolkit.Equation)
    try
        return typeof(instantiate(eq.lhs)) == typeof(instantiate(eq.rhs))
    catch
        return false
    end
end

function validate(sys::AbstractODESystem)
    all(validate.(equations(sys)))
end
