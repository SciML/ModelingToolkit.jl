Base.:*(x::Union{Num,Symbolic},y::Unitful.AbstractQuantity) = x * y

instantiate(x::Sym{Real}) = 1.0
instantiate(x::Symbolic) = oneunit(1*ModelingToolkit.vartype(x))
function instantiate(x::Num)
    x = value(x)
    if operation(x) isa Sym
        return instantiate(operation(x))
    elseif operation(x) isa Differential
        instantiate(arguments(x)[1])/instantiate(arguments(x)[1].args[1])
    else
        operation(x)(instantiate.(arguments(x))...)
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
