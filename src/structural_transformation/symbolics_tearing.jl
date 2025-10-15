using OffsetArrays: Origin
function substitution_graph(graph, slist, dlist, var_eq_matching)
    ns = length(slist)
    nd = length(dlist)
    ns == nd || error("internal error")
    newgraph = BipartiteGraph(ns, nd)
    erename = uneven_invmap(nsrcs(graph), slist)
    vrename = uneven_invmap(ndsts(graph), dlist)
    for e in 𝑠vertices(graph)
        ie = erename[e]
        ie == 0 && continue
        for v in 𝑠neighbors(graph, e)
            iv = vrename[v]
            iv == 0 && continue
            add_edge!(newgraph, ie, iv)
        end
    end
    newmatching = Matching(ns)
    for (v, e) in enumerate(var_eq_matching)
        (e === unassigned || e === SelectedState()) && continue
        iv = vrename[v]
        ie = erename[e]
        iv == 0 && continue
        ie == 0 && error("internal error")
        newmatching[iv] = ie
    end
    return DiCMOBiGraph{true}(newgraph, complete(newmatching))
end
function var_derivative_graph!(s::SystemStructure, v::Int)
    sg = g = add_vertex!(s.graph, DST)
    var_diff = add_vertex!(s.var_to_diff)
    add_edge!(s.var_to_diff, v, var_diff)
    s.solvable_graph === nothing || (sg = add_vertex!(s.solvable_graph, DST))
    @assert sg == g == var_diff
    return var_diff
end
function var_derivative!(ts::TearingState, v::Int)
    s = ts.structure
    var_diff = var_derivative_graph!(s, v)
    sys = ts.sys
    D = Differential(get_iv(sys))
    push!(ts.fullvars, D(ts.fullvars[v]))
    return var_diff
end
function eq_derivative_graph!(s::SystemStructure, eq::Int)
    add_vertex!(s.graph, SRC)
    s.solvable_graph === nothing || add_vertex!(s.solvable_graph, SRC)
    eq_diff = add_vertex!(s.eq_to_diff)
    add_edge!(s.eq_to_diff, eq, eq_diff)
    return eq_diff
end
function eq_derivative!(ts::TearingState, ieq::Int; kwargs...)
    s = ts.structure
    eq_diff = eq_derivative_graph!(s, ieq)
    sys = ts.sys
    eq = equations(ts)[ieq]
    eq = 0 ~ substitute(
        ModelingToolkit.derivative(
            eq.rhs - eq.lhs, get_iv(sys); throw_no_derivative = true), ts.param_derivative_map)
    vs = Set{SymbolicT}()
    SU.search_variables!(vs, eq.rhs)
    for v in vs
        v in ts.no_deriv_params || continue
        _original_eq = equations(ts)[ieq]
        error("""
        Encountered derivative of discrete variable `$(only(arguments(v)))` when \
        differentiating equation `$(_original_eq)`. This may indicate a model error or a \
        missing equation of the form `$v ~ ...` that defines this derivative.
        """)
    end
    push!(equations(ts), eq)
    eq_diff = length(equations(ts))
    for var in 𝑠neighbors(s.graph, ieq)
        add_edge!(s.graph, eq_diff, var)
        add_edge!(s.graph, eq_diff, s.var_to_diff[var])
    end
    s.solvable_graph === nothing ||
        find_eq_solvables!(
            ts, eq_diff; may_be_zero = true, allow_symbolic = false, kwargs...)
    return eq_diff
end
function tearing_substitution(sys::AbstractSystem; kwargs...)
    neweqs = full_equations(sys::AbstractSystem; kwargs...)
    @set! sys.eqs = neweqs
    @set! sys.schedule = nothing
end
function solve_equation(eq, var, simplify)
    rhs = value(symbolic_linear_solve(eq, var; simplify = simplify, check = false))
    SU.query(in(var), rhs) && throw(EquationSolveErrors(eq, var, rhs))
    var ~ rhs
end
function substitute_vars!(structure, subs, cache = Int[], callback! = nothing;
        exclude = ())
    @unpack graph, solvable_graph = structure
    for su in subs
        su === nothing && continue
        v, v′ = su
        eqs = 𝑑neighbors(graph, v)
        resize!(cache, length(eqs))
        for eq in copyto!(cache, eqs)
            eq in exclude && continue
            rem_edge!(graph, eq, v)
            add_edge!(graph, eq, v′)
            if BipartiteEdge(eq, v) in solvable_graph
                rem_edge!(solvable_graph, eq, v)
                add_edge!(solvable_graph, eq, v′)
            end
            callback! !== nothing && callback!(eq, su)
        end
    end
    return structure
end
function to_mass_matrix_form(neweqs, ieq, graph, fullvars, isdervar::F,
        var_to_diff) where {F}
    eq = neweqs[ieq]
    if !(eq.lhs isa Number && eq.lhs == 0)
        eq = 0 ~ eq.rhs - eq.lhs
    end
    rhs = eq.rhs
    if rhs isa SymbolicT
        dervar::Union{Nothing, Int} = nothing
        for var in 𝑠neighbors(graph, ieq)
            if isdervar(var)
                if dervar !== nothing
                    error("$eq has more than one differentiated variable!")
                end
                dervar = var
            end
        end
        dervar === nothing && return (0 ~ rhs), dervar
        new_lhs = var = fullvars[dervar]
        a, b, islinear = linear_expansion(rhs, var)
        if !islinear
            return (0 ~ rhs), nothing
        end
        new_rhs = -b / a
        return (new_lhs ~ new_rhs), invview(var_to_diff)[dervar]
    else
        if abs(rhs) > 100eps(float(rhs))
            @warn "The equation $eq is not consistent. It simplified to 0 == $rhs."
        end
        return nothing
    end
end
""""""
function substitute_derivatives_algevars!(
        ts::TearingState, neweqs::Vector{Equation}, var_eq_matching::Matching, dummy_sub::Dict{SymbolicT, SymbolicT}, iv::Union{Nothing, SymbolicT}, D::Union{Nothing, Differential, Shift})
    @unpack fullvars, sys, structure = ts
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    diff_to_var = invview(var_to_diff)
    for var in 1:length(fullvars)
        dv = var_to_diff[var]
        dv === nothing && continue
        if var_eq_matching[var] !== SelectedState()
            dd = fullvars[dv]
            v_t = setio(diff2term_with_unit(dd, iv), false, false)
            for eq in 𝑑neighbors(graph, dv)
                dummy_sub[dd] = v_t
                neweqs[eq] = substitute(neweqs[eq], dd => v_t)
            end
            fullvars[dv] = v_t
            dx = dv
            x_t = v_t
            while (ddx = var_to_diff[dx]) !== nothing
                dx_t = D(x_t)
                for eq in 𝑑neighbors(graph, ddx)
                    neweqs[eq] = substitute(neweqs[eq], fullvars[ddx] => dx_t)
                end
                fullvars[ddx] = dx_t
                dx = ddx
                x_t = dx_t
            end
            diff_to_var[dv] = nothing
        end
    end
end
""""""
function generate_derivative_variables!(
        ts::TearingState, neweqs, var_eq_matching, full_var_eq_matching,
        var_sccs, mm::Union{Nothing, SparseMatrixCLIL}, iv::Union{SymbolicT, Nothing})
    @unpack fullvars, sys, structure = ts
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    eq_var_matching = invview(var_eq_matching)
    diff_to_var = invview(var_to_diff)
    is_discrete = is_only_discrete(structure)
    linear_eqs = Dict{Int, Int}()
    if mm !== nothing
        for (i, e) in enumerate(mm.nzrows)
            linear_eqs[e] = i
        end
    end
    v_to_scc = NTuple{2, Int}[]
    resize!(v_to_scc, ndsts(graph))
    for (i, scc) in enumerate(var_sccs), (j, v) in enumerate(scc)
        v_to_scc[v] = (i, j)
    end
    v_t_dvs = NTuple{2, Int}[]
    for v in 1:length(var_to_diff)
        dv = var_to_diff[v]
        dv isa Int || continue
        solved = var_eq_matching[dv] isa Int
        solved && continue
        dd = find_duplicate_dd(dv, solvable_graph, diff_to_var, linear_eqs, mm)
        if dd === nothing
            dx = fullvars[dv]
            order, lv = var_order(dv, diff_to_var)
            x_t = is_discrete ? lower_shift_varname_with_unit(fullvars[dv], iv) :
                  lower_varname_with_unit(fullvars[lv], iv, order)
            v_t = add_dd_variable!(structure, fullvars, x_t, dv)
            dummy_eq = add_dd_equation!(structure, neweqs, 0 ~ dx - x_t, dv, v_t)
            for e in 𝑑neighbors(graph, dv)
                add_edge!(graph, e, v_t)
            end
            push!(var_eq_matching, unassigned)
            push!(full_var_eq_matching, unassigned)
            dd = dummy_eq, v_t
        end
        dummy_eq, v_t = dd
        var_to_diff[v_t] = var_to_diff[dv]
        old_matched_eq = full_var_eq_matching[dv]
        full_var_eq_matching[dv] = var_eq_matching[dv] = dummy_eq
        full_var_eq_matching[v_t] = old_matched_eq
        eq_var_matching[dummy_eq] = dv
        push!(v_t_dvs, (v_t, dv))
    end
    sccs_to_insert = similar(v_t_dvs, Tuple{Int, Vector{Int}})
    idxs_to_remove = Dict{Int, Vector{Int}}()
    for (k, (v_t, dv)) in enumerate(v_t_dvs)
        i, j = v_to_scc[dv]
        var_sccs[i][j] = v_t
        if v_t <= length(v_to_scc)
            i_, j_ = v_to_scc[v_t]
            scc_del_idxs = get!(() -> Int[], idxs_to_remove, i_)
            push!(scc_del_idxs, j_)
        end
        sccs_to_insert[k] = (i, [dv])
    end
    sort!(sccs_to_insert, by = first)
    for (i, idxs) in idxs_to_remove
        deleteat!(var_sccs[i], idxs)
    end
    new_sccs = insert_sccs(var_sccs, sccs_to_insert)
    if mm !== nothing
        @set! mm.ncols = ndsts(graph)
    end
    return new_sccs
end
""""""
function insert_sccs(
        var_sccs::Vector{Vector{Int}}, sccs_to_insert::Vector{Tuple{Int, Vector{Int}}})
    old_idx = 1
    insert_idx = 1
    new_sccs = similar(var_sccs, length(var_sccs) + length(sccs_to_insert))
    for i in eachindex(new_sccs)
        if insert_idx <= length(sccs_to_insert) && sccs_to_insert[insert_idx][1] == old_idx
            new_sccs[i] = sccs_to_insert[insert_idx][2]
            insert_idx += 1
        else
            new_sccs[i] = copy(var_sccs[old_idx])
            old_idx += 1
        end
    end
    filter!(!isempty, new_sccs)
    return new_sccs
end
""""""
function find_duplicate_dd(dv, solvable_graph, diff_to_var, linear_eqs, mm)
    for eq in 𝑑neighbors(solvable_graph, dv)
        mi = get(linear_eqs, eq, 0)
        iszero(mi) && continue
        row = @view mm[mi, :]
        nzs = nonzeros(row)
        rvs = SparseArrays.nonzeroinds(row)
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
function add_dd_variable!(s::SystemStructure, fullvars, x_t, dv)
    push!(fullvars, simplify_shifts(x_t))
    v_t = length(fullvars)
    v_t_idx = add_vertex!(s.var_to_diff)
    add_vertex!(s.graph, DST)
    add_vertex!(s.solvable_graph, DST)
    s.var_to_diff[v_t] = s.var_to_diff[dv]
    v_t
end
""""""
function add_dd_equation!(s::SystemStructure, neweqs, eq, dv, v_t)
    push!(neweqs, eq)
    add_vertex!(s.graph, SRC)
    dummy_eq = length(neweqs)
    add_edge!(s.graph, dummy_eq, dv)
    add_edge!(s.graph, dummy_eq, v_t)
    add_vertex!(s.solvable_graph, SRC)
    add_edge!(s.solvable_graph, dummy_eq, dv)
    dummy_eq
end
""""""
function generate_system_equations!(state::TearingState, neweqs::Vector{Equation},
        var_eq_matching::Matching, full_var_eq_matching::Matching,
        var_sccs::Vector{Vector{Int}}, extra_eqs_vars::NTuple{2, Vector{Int}},
        iv::Union{SymbolicT, Nothing}, D::Union{Differential, Shift, Nothing};
        simplify::Bool = false)
    @unpack fullvars, sys, structure = state
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    eq_var_matching = invview(var_eq_matching)
    diff_to_var = invview(var_to_diff)
    extra_eqs, extra_vars = extra_eqs_vars
    total_sub = Dict{SymbolicT, SymbolicT}()
    if is_only_discrete(structure)
        for (i, v) in enumerate(fullvars)
            Moshi.Match.@match v begin
                BSImpl.Term(; f) && if f isa Shift && f.steps < 0 end => begin
                    lowered = lower_shift_varname_with_unit(v, iv)
                    total_sub[v] = lowered
                    fullvars[i] = lowered
                end
                _ => nothing
            end
        end
    end
    eq_generator = EquationGenerator(state, total_sub, D, iv)
    for eq in extra_eqs
        var = eq_var_matching[eq]
        var isa Int || continue
        codegen_equation!(eq_generator, neweqs[eq], eq, var; simplify)
    end
    ispresent = let var_to_diff = var_to_diff, graph = graph
        i -> (!isempty(𝑑neighbors(graph, i)) ||
              (var_to_diff[i] !== nothing && !isempty(𝑑neighbors(graph, var_to_diff[i]))))
    end
    digraph = DiCMOBiGraph{false}(graph, var_eq_matching)
    for (i, scc) in enumerate(var_sccs)
        vscc, escc = get_sorted_scc(digraph, full_var_eq_matching, var_eq_matching, scc)
        var_sccs[i] = vscc
        if length(escc) != length(vscc)
            isempty(escc) && continue
            escc = setdiff(escc, extra_eqs)
            isempty(escc) && continue
            vscc = setdiff(vscc, extra_vars)
            isempty(vscc) && continue
        end
        for ieq in escc
            iv = eq_var_matching[ieq]
            neq = neweqs[ieq]
            codegen_equation!(eq_generator, neq, ieq, iv; simplify)
        end
    end
    for eq in extra_eqs
        var = eq_var_matching[eq]
        var isa Int && continue
        codegen_equation!(eq_generator, neweqs[eq], eq, var; simplify)
    end
    @unpack neweqs′, eq_ordering, var_ordering, solved_eqs, solved_vars = eq_generator
    is_diff_eq = .!iszero.(var_ordering)
    diff_vars = var_ordering[is_diff_eq]
    diff_vars_set = BitSet(diff_vars)
    if length(diff_vars_set) != length(diff_vars)
        error("Tearing internal error: lowering DAE into semi-implicit ODE failed!")
    end
    solved_vars_set = BitSet(solved_vars)
    offset = 1
    findnextfn = let diff_vars_set = diff_vars_set, solved_vars_set = solved_vars_set,
        diff_to_var = diff_to_var, ispresent = ispresent
        j -> !(j in diff_vars_set || j in solved_vars_set) && diff_to_var[j] === nothing &&
            ispresent(j)
    end
    for (i, v) in enumerate(var_ordering)
        v == 0 || continue
        index = findnext(findnextfn, 1:ndsts(graph), offset)
        index === nothing && break
        var_ordering[i] = index
        offset = index + 1
    end
    filter!(!iszero, var_ordering)
    var_ordering = [var_ordering; setdiff(1:ndsts(graph), var_ordering, solved_vars_set)]
    neweqs = neweqs′
    return neweqs, solved_eqs, eq_ordering, var_ordering, length(solved_vars),
    length(solved_vars_set)
end
""""""
function get_sorted_scc(
        digraph::DiCMOBiGraph, full_var_eq_matching::Matching, var_eq_matching::Matching, scc::Vector{Int})
    eq_var_matching = invview(var_eq_matching)
    scc_eqs = Int[]
    scc_solved_eqs = Int[]
    for v in scc
        e = full_var_eq_matching[v]
        if e isa Int
            push!(scc_eqs, e)
        end
        e = var_eq_matching[v]
        if e isa Int
            push!(scc_solved_eqs, e)
        end
    end
    subgraph, varmap = Graphs.induced_subgraph(digraph, scc_solved_eqs)
    scc_eqs = [varmap[reverse(topological_sort(subgraph))];
               setdiff(scc_eqs, scc_solved_eqs)]
    scc_vars = Int[]
    for e in scc_eqs
        v = eq_var_matching[e]
        if v isa Int
            push!(scc_vars, v)
        end
    end
    append!(scc_vars, setdiff(scc, scc_vars))
    return scc_vars, scc_eqs
end
""""""
struct EquationGenerator{S}
    """
    `TearingState` of the system.
    """
    state::S
    """
    Substitutions to perform in all subsequent equations. For each differential equation
    `D(x) ~ f(..)`, the substitution `D(x) => f(..)` is added to the rules.
    """
    total_sub::Dict{SymbolicT, SymbolicT}
    """
    The differential operator, or `nothing` if not applicable.
    """
    D::Union{Differential, Shift, Nothing}
    """
    The independent variable, or `nothing` if not applicable.
    """
    idep::Union{SymbolicT, Nothing}
    """
    The new generated equations of the system.
    """
    neweqs′::Vector{Equation}
    """
    `eq_ordering[i]` is the index `neweqs′[i]` was originally at in the untorn equations of
    the system. This is used to permute the state of the system into BLT sorted form.
    """
    eq_ordering::Vector{Int}
    """
    `var_ordering[i]` is the index in `state.fullvars` of the variable at the `i`th index in
    the BLT sorted form.
    """
    var_ordering::Vector{Int}
    """
    List of linearly solved (observed) equations.
    """
    solved_eqs::Vector{Equation}
    """
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
function is_dervar(eg::EquationGenerator, iv::Int)
    diff_to_var = invview(eg.state.structure.var_to_diff)
    diff_to_var[iv] !== nothing
end
""""""
function codegen_equation!(eg::EquationGenerator,
        eq::Equation, ieq::Int, iv::Union{Int, Unassigned}; simplify = false)
    @unpack state, total_sub, neweqs′, eq_ordering, var_ordering = eg
    @unpack solved_eqs, solved_vars, D, idep = eg
    @unpack fullvars, sys, structure = state
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    diff_to_var = invview(var_to_diff)
    issolvable = is_solvable(eg, ieq, iv)
    isdervar = issolvable && is_dervar(eg, iv)
    isdisc = is_only_discrete(structure)
    is_highest_diff = iv isa Int && isdervar && var_to_diff[iv] === nothing
    if issolvable && isdervar && (!isdisc || !is_highest_diff)
        var = fullvars[iv]
        isnothing(D) && throw(UnexpectedDifferentialError(equations(sys)[ieq]))
        order, lv = var_order(iv, diff_to_var)
        dx = D(simplify_shifts(fullvars[lv]))
        neweq = make_differential_equation(var, dx, eq, total_sub)
        for e in 𝑑neighbors(graph, iv)
            e == ieq && continue
            rem_edge!(graph, e, iv)
        end
        total_sub[simplify_shifts(neweq.lhs)] = neweq.rhs
        if is_only_discrete(structure)
            var_to_diff[iv] === nothing && (total_sub[var] = neweq.rhs)
        end
        push!(neweqs′, neweq)
        push!(eq_ordering, ieq)
        push!(var_ordering, diff_to_var[iv])
    elseif issolvable
        var = fullvars[iv]
        neweq = make_solved_equation(var, eq, total_sub; simplify)
        if neweq !== nothing
            if isdisc
                neweq = backshift_expr(neweq, idep)
            end
            push!(solved_eqs, neweq)
            push!(solved_vars, iv)
        end
    else
        neweq = make_algebraic_equation(eq, total_sub)
        if isdisc
            neweq = backshift_expr(neweq, idep)
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
end
""""""
function make_differential_equation(var, dx, eq, total_sub)
    v1 = Symbolics.symbolic_linear_solve(eq, var)::SymbolicT
    v2 = Symbolics.fixpoint_sub(v1, total_sub; operator = ModelingToolkit.Shift)
    v3 = simplify_shifts(v2)
    dx ~ v3
end
""""""
function make_algebraic_equation(eq, total_sub)
    rhs = eq.rhs - eq.lhs
    0 ~ simplify_shifts(Symbolics.fixpoint_sub(rhs, total_sub))
end
""""""
function make_solved_equation(var, eq, total_sub; simplify = false)
    residual = eq.lhs - eq.rhs
    a, b, islinear = linear_expansion(residual, var)
    @assert islinear
    if ModelingToolkit._iszero(a)
        @warn "Tearing: solving $eq for $var is singular!"
        return nothing
    else
        rhs = -b / a
        return var ~ simplify_shifts(Symbolics.fixpoint_sub(
            simplify ?
            Symbolics.simplify(rhs) : rhs,
            total_sub; operator = ModelingToolkit.Shift))
    end
end
""""""
function reorder_vars!(state::TearingState, var_eq_matching, var_sccs, eq_ordering,
        var_ordering, nsolved_eq, nsolved_var)
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = state.structure
    eqsperm = zeros(Int, nsrcs(graph))
    for (i, v) in enumerate(eq_ordering)
        eqsperm[v] = i
    end
    varsperm = zeros(Int, ndsts(graph))
    for (i, v) in enumerate(var_ordering)
        varsperm[v] = i
    end
    new_graph = contract_variables(graph, var_eq_matching, varsperm, eqsperm,
        nsolved_eq, nsolved_var)
    new_var_to_diff = complete(DiffGraph(length(var_ordering)))
    for (v, d) in enumerate(var_to_diff)
        v′ = varsperm[v]
        (v′ > 0 && d !== nothing) || continue
        d′ = varsperm[d]
        new_var_to_diff[v′] = d′ > 0 ? d′ : nothing
    end
    new_eq_to_diff = complete(DiffGraph(length(eq_ordering)))
    for (v, d) in enumerate(eq_to_diff)
        v′ = eqsperm[v]
        (v′ > 0 && d !== nothing) || continue
        d′ = eqsperm[d]
        new_eq_to_diff[v′] = d′ > 0 ? d′ : nothing
    end
    new_fullvars = state.fullvars[var_ordering]
    var_ordering_set = BitSet(var_ordering)
    for scc in var_sccs
        map!(Base.Fix1(getindex, varsperm), scc, scc)
        filter!(!iszero, scc)
    end
    filter!(!isempty, var_sccs)
    @set! state.structure.graph = complete(new_graph)
    @set! state.structure.var_to_diff = new_var_to_diff
    @set! state.structure.eq_to_diff = new_eq_to_diff
    @set! state.fullvars = new_fullvars
    state
end
""""""
function update_simplified_system!(
        state::TearingState, neweqs::Vector{Equation}, solved_eqs::Vector{Equation},
        dummy_sub::Dict{SymbolicT, SymbolicT}, var_sccs::Vector{Vector{Int}},
        extra_unknowns::Vector{SymbolicT}, iv::Union{SymbolicT, Nothing},
        D::Union{Differential, Shift, Nothing}; array_hack = true)
    @unpack fullvars, structure = state
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    diff_to_var = invview(var_to_diff)
    solved_vars = Set{SymbolicT}()
    if is_only_discrete(structure)
        iv = iv::SymbolicT
        D = D::Shift
        for eq in solved_eqs
            var = eq.lhs
            if isequal(eq.lhs, eq.rhs)
                var = lower_shift_varname_with_unit(D(eq.lhs), iv)
            end
            push!(solved_vars, var)
        end
        filter!(eq -> !isequal(eq.lhs, eq.rhs), solved_eqs)
    end
    ispresent = let var_to_diff = var_to_diff, graph = graph
        i -> (!isempty(𝑑neighbors(graph, i)) ||
              (var_to_diff[i] !== nothing && !isempty(𝑑neighbors(graph, var_to_diff[i]))))
    end
    sys = state.sys
    obs_sub = dummy_sub
    for eq in neweqs
        isdiffeq(eq) || continue
        obs_sub[eq.lhs] = eq.rhs
    end
    obs = [substitute(observed(sys), obs_sub); solved_eqs;
           substitute(state.additional_observed, obs_sub)]
    filterer = let diff_to_var = diff_to_var, ispresent = ispresent, fullvars = fullvars,
        solved_vars = solved_vars
        i -> diff_to_var[i] === nothing && ispresent(i) && !(fullvars[i] in solved_vars)
    end
    unknown_idxs = filter(filterer, eachindex(state.fullvars))
    unknowns = state.fullvars[unknown_idxs]
    unknowns = [unknowns; extra_unknowns]
    if is_only_discrete(structure)
        _unknowns = SymbolicT[]
        for var in unknowns
            Moshi.Match.@match var begin
                BSImpl.Term(; f, args, type, shape, metadata) && if f isa Shift && f.steps == 1 end => begin
                    push!(_unknowns, setio(args[1], false, false))
                end
                _ => push!(_unknowns, var)
            end
        end
        unknowns = _unknowns
    end
    @set! sys.unknowns = unknowns
    obs = tearing_hacks(sys, obs, unknowns, neweqs; array = array_hack)
    @set! sys.eqs = neweqs
    @set! sys.observed = obs
    if ModelingToolkit.has_schedule(sys)
        unknowns_set = BitSet(unknown_idxs)
        for scc in var_sccs
            intersect!(scc, unknowns_set)
        end
        filter!(!isempty, var_sccs)
        @set! sys.schedule = Schedule(var_sccs, dummy_sub)
    end
    if ModelingToolkit.has_isscheduled(sys)
        @set! sys.isscheduled = true
    end
    return sys
end
""""""
function var_order(dv, diff_to_var)
    order = 0
    while (dv′ = diff_to_var[dv]) !== nothing
        order += 1
        dv = dv′
    end
    order, dv
end
""""""
function tearing_reassemble(state::TearingState, var_eq_matching::Matching,
        full_var_eq_matching::Matching, var_sccs::Vector{Vector{Int}}; simplify = false, mm,
        array_hack = true, fully_determined = true)
    extra_eqs_vars = get_extra_eqs_vars(
        state, var_eq_matching, full_var_eq_matching, fully_determined)
    neweqs = collect(equations(state))
    dummy_sub = Dict{SymbolicT, SymbolicT}()
    if ModelingToolkit.has_iv(state.sys)
        iv = get_iv(state.sys)
        if !is_only_discrete(state.structure)
            D = Differential(iv)
        else
            D = Shift(iv, 1)
        end
    else
        iv = D = nothing
    end
    extra_unknowns = state.fullvars[extra_eqs_vars[2]]
    if is_only_discrete(state.structure)
        var_sccs = add_additional_history!(
            state, var_eq_matching, full_var_eq_matching, var_sccs, iv)
    end
    substitute_derivatives_algevars!(state, neweqs, var_eq_matching, dummy_sub, iv, D)
    var_sccs = generate_derivative_variables!(
        state, neweqs, var_eq_matching, full_var_eq_matching, var_sccs, mm, iv)
    neweqs, solved_eqs,
    eq_ordering,
    var_ordering,
    nelim_eq,
    nelim_var = generate_system_equations!(
        state, neweqs, var_eq_matching, full_var_eq_matching,
        var_sccs, extra_eqs_vars, iv, D; simplify)
    state = reorder_vars!(
        state, var_eq_matching, var_sccs, eq_ordering, var_ordering, nelim_eq, nelim_var)
    sys = update_simplified_system!(state, neweqs, solved_eqs, dummy_sub, var_sccs,
        extra_unknowns, iv, D; array_hack)
    @set! state.sys = sys
    @set! sys.tearing_state = state
    return invalidate_cache!(sys)
end
""""""
function add_additional_history!(
        state::TearingState, var_eq_matching::Matching,
        full_var_eq_matching::Matching, var_sccs::Vector{Vector{Int}}, iv::Union{SymbolicT, Nothing})
    iv === nothing && return var_sccs
    iv = iv::SymbolicT
    @unpack fullvars, sys, structure = state
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    diff_to_var = invview(var_to_diff)
    v_to_scc = NTuple{2, Int}[]
    resize!(v_to_scc, ndsts(graph))
    for (i, scc) in enumerate(var_sccs), (j, v) in enumerate(scc)
        v_to_scc[v] = (i, j)
    end
    vars_to_backshift = BitSet()
    for ivar in 1:length(fullvars)
        ieq = var_eq_matching[ivar]
        ieq isa SelectedState || continue
        diff_to_var[ivar] === nothing || continue
        push!(vars_to_backshift, ivar)
    end
    inserts = Tuple{Int, Vector{Int}}[]
    for var in vars_to_backshift
        add_backshifted_var!(state, var, iv)
        push!(var_eq_matching, SelectedState())
        push!(full_var_eq_matching, unassigned)
        push!(inserts, (v_to_scc[var][1], [length(fullvars)]))
    end
    sort!(inserts, by = first)
    new_sccs = insert_sccs(var_sccs, inserts)
    return new_sccs
end
""""""
function add_backshifted_var!(state::TearingState, ivar::Int, iv)
    @unpack fullvars, structure = state
    @unpack var_to_diff, graph, solvable_graph = structure
    var = fullvars[ivar]
    newvar = simplify_shifts(Shift(iv, -1)(var))
    push!(fullvars, newvar)
    inewvar = add_vertex!(var_to_diff)
    add_edge!(var_to_diff, inewvar, ivar)
    add_vertex!(graph, DST)
    add_vertex!(solvable_graph, DST)
    return inewvar
end
""""""
function backshift_expr(ex, iv)
    ex isa SymbolicT || return ex
    return descend_lower_shift_varname_with_unit(
        simplify_shifts(distribute_shift(Shift(iv, -1)(ex))), iv)
end
function backshift_expr(ex::Equation, iv)
    return backshift_expr(ex.lhs, iv) ~ backshift_expr(ex.rhs, iv)
end
""""""
function get_extra_eqs_vars(
        state::TearingState, var_eq_matching::Matching, full_var_eq_matching::Matching, fully_determined::Bool)
    fully_determined && return Int[], Int[]
    extra_eqs = Int[]
    extra_vars = Int[]
    full_eq_var_matching = invview(full_var_eq_matching)
    for v in 𝑑vertices(state.structure.graph)
        eq = full_var_eq_matching[v]
        eq isa Int && continue
        var_eq_matching[v] === unassigned || continue
        push!(extra_vars, v)
    end
    for eq in 𝑠vertices(state.structure.graph)
        v = full_eq_var_matching[eq]
        v isa Int && continue
        push!(extra_eqs, eq)
    end
    return extra_eqs, extra_vars
end
""""""
function tearing_hacks(sys, obs, unknowns, neweqs; array = true)
    arr_obs_occurrences = Dict{SymbolicT, Int}()
    for (i, eq) in enumerate(obs)
        lhs = eq.lhs
        rhs = eq.rhs
        array || continue
        iscall(lhs) || continue
        operation(lhs) === getindex || continue
        SU.shape(lhs) isa SU.Unknown && continue
        arg1 = arguments(lhs)[1]
        cnt = get(arr_obs_occurrences, arg1, 0)
        arr_obs_occurrences[arg1] = cnt + 1
        continue
    end
    for sym in unknowns
        iscall(sym) || continue
        operation(sym) === getindex || continue
        SU.shape(sym) isa SU.Unknown && continue
        arg1 = arguments(sym)[1]
        cnt = get(arr_obs_occurrences, arg1, 0)
        cnt == 0 && continue
        arr_obs_occurrences[arg1] = cnt + 1
    end
    obs_arr_eqs = Equation[]
    for (arrvar, cnt) in arr_obs_occurrences
        cnt == length(arrvar) || continue
        firstind = Tuple(first(eachindex(arrvar)))
        scal = [arrvar[i] for i in eachindex(arrvar)]
        rhs = change_origin(firstind, scal)
        push!(obs_arr_eqs, arrvar ~ rhs)
    end
    append!(obs, obs_arr_eqs)
    return obs
end
function change_origin(origin, arr)
    if all(isone, origin)
        return arr
    end
    return Origin(origin)(arr)
end
@register_array_symbolic change_origin(origin::Any, arr::AbstractArray) begin
    size = size(arr)
    eltype = eltype(arr)
    ndims = ndims(arr)
end
function tearing(state::TearingState; kwargs...)
    state.structure.solvable_graph === nothing && find_solvables!(state; kwargs...)
    complete!(state.structure)
    tearing_with_dummy_derivatives(state.structure, ())
end
""""""
function tearing(sys::AbstractSystem, state = TearingState(sys); mm = nothing,
        simplify = false, array_hack = true, fully_determined = true, kwargs...)
    var_eq_matching, full_var_eq_matching, var_sccs, can_eliminate = tearing(state)
    invalidate_cache!(tearing_reassemble(
        state, var_eq_matching, full_var_eq_matching, var_sccs; mm,
        simplify, array_hack, fully_determined))
end
""""""
function dummy_derivative(sys, state = TearingState(sys); simplify = false,
        mm = nothing, array_hack = true, fully_determined = true, kwargs...)
    jac = let state = state
        (eqs, vars) -> begin
            symeqs = EquationsView(state)[eqs]
            Symbolics.jacobian((x -> x.rhs).(symeqs), state.fullvars[vars])
        end
    end
    state_priority = let state = state
        var -> begin
            p = 0.0
            var_to_diff = state.structure.var_to_diff
            diff_to_var = invview(var_to_diff)
            while var_to_diff[var] !== nothing
                var = var_to_diff[var]
            end
            while true
                p = max(p, ModelingToolkit.state_priority(state.fullvars[var]))
                (var = diff_to_var[var]) === nothing && break
            end
            p
        end
    end
    var_eq_matching, full_var_eq_matching, var_sccs,
    can_eliminate, summary = dummy_derivative_graph!(
        state, jac; state_priority,
        kwargs...)
    tearing_reassemble(state, var_eq_matching, full_var_eq_matching, var_sccs;
        simplify, mm, array_hack, fully_determined)
end
