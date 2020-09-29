export alias_elimination

function flatten(sys::ODESystem)
    if isempty(sys.systems)
        return sys
    else
        return ODESystem(equations(sys),
                         independent_variable(sys),
                         observed=observed(sys))
    end
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

function substitute_aliases(diffeqs, dict)
    lhss(diffeqs) .~ fixpoint_sub.(rhss(diffeqs), (dict,))
end

function make_lhs_0(eq)
    if eq.lhs isa Constant && iszero(eq.lhs)
        return eq
    else
        0 ~ eq.lhs - eq.rhs
    end
end

function alias_elimination(sys::ODESystem)
    eqs = vcat(equations(sys), observed(sys))

    # make all algebraic equations have 0 on LHS
    eqs = map(eqs) do eq
        if eq.lhs isa Operation && eq.lhs.op isa Differential
            eq
        else
            make_lhs_0(eq)
        end
    end

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

    vars = map(x->x(sys.iv()), eliminate)

    outputs = solve_for(eqs[alg_idxs], vars)

    diffeqs = eqs[setdiff(1:length(eqs), alg_idxs)]

    diffeqs′ = substitute_aliases(diffeqs, Dict(vars .=> outputs))

    ODESystem(diffeqs′, sys.iv(), new_stateops, parameters(sys), observed=vars .~ outputs)
end

