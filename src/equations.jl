export Equation


mutable struct Equation
    lhs::Expression
    rhs::Expression
end
Base.broadcastable(eq::Equation) = Ref(eq)

Base.:~(lhs::Expression, rhs::Expression) = Equation(lhs, rhs)
Base.:~(lhs::Expression, rhs::Number    ) = Equation(lhs, rhs)
Base.:~(lhs::Number    , rhs::Expression) = Equation(lhs, rhs)


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

get_args(O::Operation) = O.args
get_args(eq::Equation) = Expression[eq.lhs, eq.rhs]
function vars!(vars, op)
    for arg ∈ get_args(op)
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
