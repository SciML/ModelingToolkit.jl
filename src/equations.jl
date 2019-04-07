export Equation


struct Equation
    lhs::Expression
    rhs::Expression
end
Base.:(==)(a::Equation, b::Equation) = isequal((a.lhs, a.rhs), (b.lhs, b.rhs))

Base.:~(lhs::Expression, rhs::Expression) = Equation(lhs, rhs)
Base.:~(lhs::Expression, rhs::Number    ) = Equation(lhs, rhs)
Base.:~(lhs::Number    , rhs::Expression) = Equation(lhs, rhs)

_is_parameter(iv) = (O::Operation) -> O.op.known && !isequal(O, iv)
_is_known(O::Operation) = O.op.known
_is_unknown(O::Operation) = !O.op.known

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
vars(exprs) = foldl(vars!, exprs; init = Set{Variable}())
function vars!(vars, O)
    isa(O, Operation) || return vars
    for arg ∈ O.args
        if isa(arg, Operation)
            isa(arg.op, Variable) && push!(vars, arg.op)
            vars!(vars, arg)
        end
    end

    return vars
end
