function tearing_substitution(sys::AbstractSystem; kwargs...)
    neweqs = full_equations(sys::AbstractSystem; kwargs...)
    @set! sys.eqs = neweqs
    @unpack graph, solvable_graph = structure
    for su in subs
        su === nothing && continue
        for eq in copyto!(cache, eqs)
            eq in exclude && continue
            rem_edge!(graph, eq, v)
        end
    end
    return structure
    eq = neweqs[ieq]
    if !(eq.lhs isa Number && eq.lhs == 0)
        eq = 0 ~ eq.rhs - eq.lhs
        dervar::Union{Nothing, Int} = nothing
        for var in 𝑠neighbors(graph, ieq)
            if isdervar(var)
                dervar = var
            end
        end
        return nothing
    end
end
""""""
function substitute_derivatives_algevars!(
        ts::TearingState, neweqs::Vector{Equation}, var_eq_matching::Matching, dummy_sub::Dict{SymbolicT, SymbolicT}, iv::Union{Nothing, SymbolicT}, D::Union{Nothing, Differential, Shift})
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
    total_sub = Dict{SymbolicT, SymbolicT}()
    for eq in extra_eqs
        var = eq_var_matching[eq]
        vscc, escc = get_sorted_scc(digraph, full_var_eq_matching, var_eq_matching, scc)
        if e isa Int
            push!(scc_eqs, e)
            push!(scc_solved_eqs, e)
        end
    end
    scc_vars = Int[]
    for e in scc_eqs
        v = eq_var_matching[e]
    end
    append!(scc_vars, setdiff(scc, scc_vars))
    return scc_vars, scc_eqs
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
    diff_to_var = invview(eg.state.structure.var_to_diff)
    diff_to_var[iv] !== nothing
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
function Base.showerror(io::IO, err::UnexpectedDifferentialError)
    error("Differential found in a non-differential system. Likely this is a bug in the construction of an initialization system. Please report this issue with a reproducible example. Offending equation: $(err.eq)")
    if ModelingToolkit._iszero(a)
        @warn "Tearing: solving $eq for $var is singular!"
        return nothing
    else
        rhs = -b / a
        return var ~ simplify_shifts(Symbolics.fixpoint_sub(
            Symbolics.simplify(rhs) : rhs,
            total_sub; operator = ModelingToolkit.Shift))
    end
end
function reorder_vars!(state::TearingState, var_eq_matching, var_sccs, eq_ordering,
        var_ordering, nsolved_eq, nsolved_var)
    for (i, v) in enumerate(eq_ordering)
    end
    if is_only_discrete(structure)
        for eq in solved_eqs
            if isequal(eq.lhs, eq.rhs)
            end
        end
    end
    while (dv′ = diff_to_var[dv]) !== nothing
    end
end
function tearing_reassemble(state::TearingState, var_eq_matching::Matching,
        array_hack = true, fully_determined = true)
    extra_eqs_vars = get_extra_eqs_vars(
        state, var_eq_matching, full_var_eq_matching, fully_determined)
    if ModelingToolkit.has_iv(state.sys)
        if !is_only_discrete(state.structure)
        end
    end
    state = reorder_vars!(
        state, var_eq_matching, var_sccs, eq_ordering, var_ordering, nelim_eq, nelim_var)
    sys = update_simplified_system!(state, neweqs, solved_eqs, dummy_sub, var_sccs,
        extra_unknowns, iv, D; array_hack)
end
function add_additional_history!(
        state::TearingState, var_eq_matching::Matching, full_var_eq_matching::Matching, fully_determined::Bool)
end
function tearing_hacks(sys, obs, unknowns, neweqs; array = true)
    for (i, eq) in enumerate(obs)
    end
    for (arrvar, cnt) in arr_obs_occurrences
    end
    invalidate_cache!(tearing_reassemble(
        simplify, array_hack, fully_determined))
end
function dummy_derivative(sys, state = TearingState(sys); simplify = false,
        mm = nothing, array_hack = true, fully_determined = true, kwargs...)
    jac = let state = state
        (eqs, vars) -> begin
            while var_to_diff[var] !== nothing
            end
        end
    end
    tearing_reassemble(state, var_eq_matching, full_var_eq_matching, var_sccs;
        simplify, mm, array_hack, fully_determined)
end
