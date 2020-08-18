export alias_elimination

function substitute_aliases(diffeqs, outputs)
    lhss(diffeqs) .~ substitute.(rhss(diffeqs), (Dict(lhss(outputs) .=> rhss(outputs)),))
end

function alias_elimination(sys::ODESystem)
    eqs = equations(sys)

    new_stateops = map(eqs) do eq
        if eq.lhs isa Operation && eq.lhs.op isa Differential
            get_variables(eq.lhs)
        else
            []
        end
    end |> Iterators.flatten |> collect |> unique
    newstates = convert.(Variable, new_stateops)

    alg_idxs = findall(x->x.lhs isa Constant && iszero(x.lhs), eqs)

    eliminate = setdiff(states(sys), newstates)
    outputs = solve_for(eqs[alg_idxs], map(x->x(sys.iv()), eliminate))

    diffeqs = eqs[setdiff(1:length(eqs), alg_idxs)]

    diffeqs′ = substitute_aliases(diffeqs, outputs)

    ODESystem(diffeqs′, sys.iv(), new_stateops, parameters(sys), observed=outputs)
end

