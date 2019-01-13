export Equation


mutable struct Equation
    lhs::Expression
    rhs::Expression
end
Base.broadcastable(eq::Equation) = Ref(eq)

Base.:~(lhs::Expression, rhs::Expression) = Equation(lhs, rhs)
Base.:~(lhs::Expression, rhs::Number    ) = Equation(lhs, rhs)
Base.:~(lhs::Number    , rhs::Expression) = Equation(lhs, rhs)


_is_derivative(x::Variable) = x.diff !== nothing
_is_dependent(x::Variable) = x.subtype === :Unknown && !isempty(x.dependents)
_is_parameter(ivs) = x -> x.subtype === :Parameter && x ∉ ivs
_subtype(subtype::Symbol) = x -> x.subtype === subtype

function extract_elements(eqs, predicates)
    result = [Variable[] for p ∈ predicates]
    vars = foldl(vars!, eqs; init=Set{Variable}())

    for var ∈ vars
        for (i, p) ∈ enumerate(predicates)
            p(var) && (push!(result[i], var); break)
        end
    end

    return result
end

function vars!(vars, op)
    args = isa(op, Equation) ? Expression[op.lhs, op.rhs] : op.args

    for arg ∈ args
        if isa(arg, Operation)
            vars!(vars, arg)
        elseif isa(arg, Variable)
            push!(vars, arg)
            for dep ∈ arg.dependents
                push!(vars, dep)
            end
        end
    end

    return vars
end
