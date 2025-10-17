function masked_cumsum!(A::Vector)
    acc = zero(eltype(A))
    for i in eachindex(A)
        iszero(A[i]) && continue
        A[i] = (acc += A[i])
    end
end
function free_equations(graph, vars_scc, var_eq_matching, varfilter::F) where {F}
    ne = nsrcs(graph)
    seen_eqs = falses(ne)
    for vars in vars_scc, var in vars
        varfilter(var) || continue
        ieq = var_eq_matching[var]
        if ieq isa Int
            seen_eqs[ieq] = true
        end
    end
    findall(!, seen_eqs)
end
