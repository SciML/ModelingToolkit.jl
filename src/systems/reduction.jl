export alias_elimination

function flatten(sys::ODESystem)
    ODESystem(equations(sys),
              independent_variable(sys),
              observed=observed(sys))
end


using SymbolicUtils: Rewriters

function fixpoint_sub(x, dict)
    y = substitute(x, dict)
    while !isequal(x, y)
        y = x
        x = substitute(y, dict)
    end

    return x
end

function substitute_aliases(diffeqs, outputs)
    lhss(diffeqs) .~ fixpoint_sub.(rhss(diffeqs), (Dict(lhss(outputs) .=> rhss(outputs)),))
end

function make_lhs_0(eq)
    if eq.lhs isa Constant && iszero(eq.lhs)
        return eq
    else
        0 ~ eq.lhs - eq.rhs
    end
end

function alias_elimination(sys::ODESystem)
    eqs = vcat(equations(sys),
               make_lhs_0.(observed(sys)))

    new_stateops = map(eqs) do eq
        if eq.lhs isa Operation && eq.lhs.op isa Differential
            get_variables(eq.lhs)
        else
            []
        end
    end |> Iterators.flatten |> collect |> unique

    all_vars = map(eqs) do eq
        filter(x->!isparameter(x.op), get_variables(eq.rhs))
    end |> Iterators.flatten |> collect |> unique

    newstates = convert.(Variable, new_stateops)

    alg_idxs = findall(x->x.lhs isa Constant && iszero(x.lhs), eqs)

    eliminate = setdiff(convert.(Variable, all_vars), newstates)

    outputs = solve_for(eqs[alg_idxs], map(x->x(sys.iv()), eliminate))

    diffeqs = eqs[setdiff(1:length(eqs), alg_idxs)]

    diffeqs′ = substitute_aliases(diffeqs, outputs)

    ODESystem(diffeqs′, sys.iv(), new_stateops, parameters(sys), observed=outputs)
end

