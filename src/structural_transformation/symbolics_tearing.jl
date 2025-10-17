function substitute_derivatives_algevars!(
        ts::TearingState, neweqs::Vector{Equation}, var_eq_matching::Matching, dummy_sub::Dict{SymbolicT, SymbolicT}, iv::Union{Nothing, SymbolicT}, D::Union{Nothing, Differential})
    for var in 1:length(fullvars)
        dv = var_to_diff[var]
        dv === nothing && continue
        if var_eq_matching[var] !== SelectedState()
            dd = fullvars[dv]
            v_t = setio(diff2term_with_unit(dd, iv), false, false)
            diff_to_var[dv] = nothing
        end
    end
    for eq in 𝑑neighbors(solvable_graph, dv)
        mi = get(linear_eqs, eq, 0)
        iszero(mi) && continue
        if length(nzs) == 2 &&
           (abs(nzs[1]) == 1 && nzs[1] == -nzs[2]) &&
           (v_t = rvs[1] == dv ? rvs[2] : rvs[1];
               diff_to_var[v_t] === nothing)
            @assert dv in rvs
            return eq, v_t
        end
    end
    return nothing
end
""""""
struct EquationGenerator{S}
    state::S
    """
    List of linearly solved (observed) equations.
    `eq_ordering` for `solved_eqs`.
    """
    solved_vars::Vector{Int}
end
function EquationGenerator(state, total_sub, D, idep)
    EquationGenerator(
        state, total_sub, D, idep, Equation[], Int[], Int[], Equation[], Int[])
end
""""""
function is_solvable(eg::EquationGenerator, ieq, iv)
    solvable_graph = eg.state.structure.solvable_graph
    return ieq isa Int && iv isa Int && BipartiteEdge(ieq, iv) in solvable_graph
end
""""""
function codegen_equation!(eg::EquationGenerator,
        eq::Equation, ieq::Int, iv::Union{Int, Unassigned}; simplify = false)
    if issolvable && isdervar && (!isdisc || !is_highest_diff)
        var = fullvars[iv]
        isnothing(D) && throw(UnexpectedDifferentialError(equations(sys)[ieq]))
        neweq = make_solved_equation(var, eq, total_sub; simplify)
        if neweq !== nothing
            if isdisc
                neweq = backshift_expr(neweq, idep)
            end
            push!(solved_eqs, neweq)
        end
        push!(neweqs′, neweq)
        push!(eq_ordering, ieq)
        push!(var_ordering, 0)
    end
end
""""""
struct UnexpectedDifferentialError
    eq::Equation
end
function tearing_reassemble(state::TearingState, var_eq_matching::Matching,
        array_hack = true, fully_determined = true)
    extra_eqs_vars = get_extra_eqs_vars(
        state, var_eq_matching, full_var_eq_matching, fully_determined)
    sys = update_simplified_system!(state, neweqs, solved_eqs, dummy_sub, var_sccs,
        extra_unknowns, iv, D; array_hack)
end
function dummy_derivative(sys, state = TearingState(sys); simplify = false,
        mm = nothing, array_hack = true, fully_determined = true, kwargs...)
    tearing_reassemble(state, var_eq_matching, full_var_eq_matching, var_sccs;
        simplify, mm, array_hack, fully_determined)
end
