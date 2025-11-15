using OffsetArrays: Origin

# N.B. assumes `slist` and `dlist` are unique
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
    # the new equation is created by differentiating `eq`
    eq_diff = add_vertex!(s.eq_to_diff)
    add_edge!(s.eq_to_diff, eq, eq_diff)
    return eq_diff
end

function eq_derivative!(ts::TearingState, ieq::Int; kwargs...)
    s = ts.structure

    eq_diff = eq_derivative_graph!(s, ieq)

    sys = ts.sys
    eq = equations(ts)[ieq]
    eq = 0 ~ fast_substitute(
        ModelingToolkit.derivative(
            eq.rhs - eq.lhs, get_iv(sys); throw_no_derivative = true), ts.param_derivative_map)

    vs = ModelingToolkit.vars(eq.rhs)
    for v in vs
        # parameters with unknown derivatives have a value of `nothing` in the map,
        # so use `missing` as the default.
        get(ts.param_derivative_map, v, missing) === nothing || continue
        _original_eq = equations(ts)[ieq]
        error("""
        Encountered derivative of discrete variable `$(only(arguments(v)))` when \
        differentiating equation `$(_original_eq)`. This may indicate a model error or a \
        missing equation of the form `$v ~ ...` that defines this derivative.
        """)
    end

    push!(equations(ts), eq)
    # Analyze the new equation and update the graph/solvable_graph
    # First, copy the previous incidence and add the derivative terms.
    # That's a superset of all possible occurrences. find_solvables! will
    # remove those that doesn't actually occur.
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
    # @set! sys.substitutions = nothing
    @set! sys.schedule = nothing
end

function solve_equation(eq, var, simplify)
    rhs = value(symbolic_linear_solve(eq, var; simplify = simplify, check = false))
    occursin(var, rhs) && throw(EquationSolveErrors(eq, var, rhs))
    var ~ rhs
end

function substitute_vars!(structure, subs, cache = Int[], callback! = nothing;
        exclude = ())
    @unpack graph, solvable_graph = structure
    for su in subs
        su === nothing && continue
        v, v′ = su
        eqs = 𝑑neighbors(graph, v)
        # Note that the iterator is not robust under deletion and
        # insertion. Hence, we have a copy here.
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
    if rhs isa Symbolic
        # Check if the RHS is solvable in all unknown variable derivatives and if those
        # the linear terms for them are all zero. If so, move them to the
        # LHS.
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
        # 0 ~ a * D(x) + b
        # D(x) ~ -b/a
        a, b, islinear = linear_expansion(rhs, var)
        if !islinear
            return (0 ~ rhs), nothing
        end
        new_rhs = -b / a
        return (new_lhs ~ new_rhs), invview(var_to_diff)[dervar]
    else # a number
        if abs(rhs) > 100eps(float(rhs))
            @warn "The equation $eq is not consistent. It simplified to 0 == $rhs."
        end
        return nothing
    end
end

#=
function check_diff_graph(var_to_diff, fullvars)
    diff_to_var = invview(var_to_diff)
    for (iv, v) in enumerate(fullvars)
        ov, order = var_from_nested_derivative(v)
        graph_order = 0
        vv = iv
        while true
            vv = diff_to_var[vv]
            vv === nothing && break
            graph_order += 1
        end
        @assert graph_order==order "graph_order: $graph_order, order: $order for variable $v"
    end
end
=#

"""
Replace derivatives of non-selected unknown variables by dummy derivatives.

State selection may determine that some differential variables are
algebraic variables in disguise. The derivative of such variables are
called dummy derivatives.

`SelectedState` information is no longer needed after this function is called.
State selection is done. All non-differentiated variables are algebraic
variables, and all variables that appear differentiated are differential variables.
"""
function substitute_derivatives_algevars!(
        ts::TearingState, neweqs, var_eq_matching, dummy_sub; iv = nothing, D = nothing)
    @unpack fullvars, sys, structure = ts
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    diff_to_var = invview(var_to_diff)

    for var in 1:length(fullvars)
        dv = var_to_diff[var]
        dv === nothing && continue
        if var_eq_matching[var] !== SelectedState()
            dd = fullvars[dv]
            v_t = setio(diff2term_with_unit(unwrap(dd), unwrap(iv)), false, false)
            for eq in 𝑑neighbors(graph, dv)
                dummy_sub[dd] = v_t
                neweqs[eq] = fast_substitute(neweqs[eq], dd => v_t)
            end
            fullvars[dv] = v_t
            # If we have:
            # x -> D(x) -> D(D(x))
            # We need to to transform it to:
            # x   x_t -> D(x_t)
            # update the structural information
            dx = dv
            x_t = v_t
            while (ddx = var_to_diff[dx]) !== nothing
                dx_t = D(x_t)
                for eq in 𝑑neighbors(graph, ddx)
                    neweqs[eq] = fast_substitute(neweqs[eq], fullvars[ddx] => dx_t)
                end
                fullvars[ddx] = dx_t
                dx = ddx
                x_t = dx_t
            end
            diff_to_var[dv] = nothing
        end
    end
end

#=
There are three cases where we want to generate new variables to convert
the system into first order (semi-implicit) ODEs.

1. To first order:
Whenever higher order differentiated variable like `D(D(D(x)))` appears,
we introduce new variables `x_t`, `x_tt`, and `x_ttt` and new equations
```
D(x_tt) = x_ttt
D(x_t) = x_tt
D(x) = x_t
```
and replace `D(x)` to `x_t`, `D(D(x))` to `x_tt`, and `D(D(D(x)))` to
`x_ttt`.

2. To implicit to semi-implicit ODEs:
2.1: Unsolvable derivative:
If one derivative variable `D(x)` is unsolvable in all the equations it
appears in, then we introduce a new variable `x_t`, a new equation
```
D(x) ~ x_t
```
and replace all other `D(x)` to `x_t`.

2.2: Solvable derivative:
If one derivative variable `D(x)` is solvable in at least one of the
equations it appears in, then we introduce a new variable `x_t`. One of
the solvable equations must be in the form of `0 ~ L(D(x), u...)` and
there exists a function `l` such that `D(x) ~ l(u...)`. We should replace
it to
```
0 ~ x_t - l(u...)
D(x) ~ x_t
```
and replace all other `D(x)` to `x_t`.

Observe that we don't need to actually introduce a new variable `x_t`, as
the above equations can be lowered to
```
x_t := l(u...)
D(x) ~ x_t
```
where `:=` denotes assignment.

As a final note, in all the above cases where we need to introduce new
variables and equations, don't add them when they already exist.

###### DISCRETE SYSTEMS #######

Documenting the differences to structural simplification for discrete systems:

In discrete systems everything gets shifted forward a timestep by `shift_discrete_system`
in order to properly generate the difference equations.

In the system x(k) ~ x(k-1) + x(k-2), becomes Shift(t, 1)(x(t)) ~ x(t) + Shift(t, -1)(x(t))

The lowest-order term is Shift(t, k)(x(t)), instead of x(t). As such we actually want
dummy variables for the k-1 lowest order terms instead of the k-1 highest order terms.

Shift(t, -1)(x(t)) -> x\_{t-1}(t)

Since Shift(t, -1)(x) is not a derivative, it is directly substituted in `fullvars`.
No equation or variable is added for it.

For ODESystems D(D(D(x))) in equations is recursively substituted as D(x) ~ x_t, D(x_t) ~ x_tt, etc.
The analogue for discrete systems, Shift(t, 1)(Shift(t,1)(Shift(t,1)(Shift(t, -3)(x(t)))))
does not actually appear. So `total_sub` in generate_system_equations` is directly
initialized with all of the lowered variables `Shift(t, -3)(x) -> x_t-3(t)`, etc.
=#
"""
Generate new derivative variables for the system.

Effects on the system structure:
- fullvars: add the new derivative variables x_t
- neweqs: add the identity equations for the new variables, D(x) ~ x_t
- graph: update graph with the new equations and variables, and their connections
- solvable_graph: mark the new equation as solvable for `D(x)`
- var_eq_matching: match D(x) to the added identity equation `D(x) ~ x_t`
- full_var_eq_matching: match `x_t` to the equation that `D(x)` used to match to, and
  match `D(x)` to `D(x) ~ x_t`
- var_sccs: Replace `D(x)` in its SCC by `x_t`, and add `D(x)` in its own SCC. Return
  the new list of SCCs.
"""
function generate_derivative_variables!(
        ts::TearingState, neweqs, var_eq_matching, full_var_eq_matching,
        var_sccs; mm, iv = nothing, D = nothing)
    @unpack fullvars, sys, structure = ts
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    eq_var_matching = invview(var_eq_matching)
    diff_to_var = invview(var_to_diff)
    is_discrete = is_only_discrete(structure)
    linear_eqs = mm === nothing ? Dict{Int, Int}() :
                 Dict(reverse(en) for en in enumerate(mm.nzrows))

    # We need the inverse mapping of `var_sccs` to update it efficiently later.
    v_to_scc = Vector{NTuple{2, Int}}(undef, ndsts(graph))
    for (i, scc) in enumerate(var_sccs), (j, v) in enumerate(scc)

        v_to_scc[v] = (i, j)
    end
    # Pairs of `(x_t, dx)` added below
    v_t_dvs = NTuple{2, Int}[]

    # For variable x, make dummy derivative x_t if the
    # derivative is in the system
    for v in 1:length(var_to_diff)
        dv = var_to_diff[v]
        # if the variable is not differentiated, there is nothing to do
        dv isa Int || continue
        # if we will solve for the differentiated variable, there is nothing to do
        solved = var_eq_matching[dv] isa Int
        solved && continue

        # If there's `D(x) = x_t` already, update mappings and continue without
        # adding new equations/variables
        dd = find_duplicate_dd(dv, solvable_graph, diff_to_var, linear_eqs, mm)
        if dd === nothing
            # there is no such pre-existing equation
            # generate the dummy derivative variable
            dx = fullvars[dv]
            order, lv = var_order(dv, diff_to_var)
            x_t = is_discrete ? lower_shift_varname_with_unit(fullvars[dv], iv) :
                  lower_varname_with_unit(fullvars[lv], iv, order)

            # Add `x_t` to the graph
            v_t = add_dd_variable!(structure, fullvars, x_t, dv)
            # Add `D(x) - x_t ~ 0` to the graph
            dummy_eq = add_dd_equation!(structure, neweqs, 0 ~ dx - x_t, dv, v_t)
            # Update graph to say, all the equations featuring D(x) also feature x_t
            for e in 𝑑neighbors(graph, dv)
                add_edge!(graph, e, v_t)
            end
            # Update matching
            push!(var_eq_matching, unassigned)
            push!(full_var_eq_matching, unassigned)

            # We also need to substitute all occurrences of `D(x)` with `x_t` in all equations
            # except `dummy_eq`, but that is handled in `generate_system_equations!` since
            # we will solve for `D(x) ~ x_t` and add it to the substitution map.
            dd = dummy_eq, v_t
        end
        # there is a duplicate `D(x) ~ x_t` equation
        # `dummy_eq` is the index of the equation
        # `v_t` is the dummy derivative variable
        dummy_eq, v_t = dd
        var_to_diff[v_t] = var_to_diff[dv]
        old_matched_eq = full_var_eq_matching[dv]
        full_var_eq_matching[dv] = var_eq_matching[dv] = dummy_eq
        full_var_eq_matching[v_t] = old_matched_eq
        eq_var_matching[dummy_eq] = dv
        push!(v_t_dvs, (v_t, dv))
    end

    # tuples of (index, scc) indicating that `scc` has to be inserted at
    # index `index` in `var_sccs`. Same length as `v_t_dvs` because we will
    # have one new SCC per new variable.
    sccs_to_insert = similar(v_t_dvs, Tuple{Int, Vector{Int}})
    # mapping of SCC index to indexes in the SCC to delete
    idxs_to_remove = Dict{Int, Vector{Int}}()
    for (k, (v_t, dv)) in enumerate(v_t_dvs)
        # replace `dv` with `v_t`
        i, j = v_to_scc[dv]
        var_sccs[i][j] = v_t
        if v_t <= length(v_to_scc)
            # v_t wasn't added by this process, it was already present. Which
            # means we need to remove it from whatever SCC it is in, since it is
            # now in this one
            i_, j_ = v_to_scc[v_t]
            scc_del_idxs = get!(() -> Int[], idxs_to_remove, i_)
            push!(scc_del_idxs, j_)
        end
        # `dv` still needs to be present in some SCC. Since we solve for `dv` from
        # `0 ~ D(x) - x_t`, it is in its own SCC. This new singleton SCC is solved
        # immediately before the one that `dv` used to be in (`i`)
        sccs_to_insert[k] = (i, [dv])
    end
    sort!(sccs_to_insert, by = first)
    # remove the idxs we need to remove
    for (i, idxs) in idxs_to_remove
        deleteat!(var_sccs[i], idxs)
    end
    new_sccs = insert_sccs(var_sccs, sccs_to_insert)

    if mm !== nothing
        @set! mm.ncols = ndsts(graph)
    end

    return new_sccs
end

"""
    $(TYPEDSIGNATURES)

Given a list of SCCs and a list of SCCs to insert at specific indices, insert them and
return the new SCC vector.
"""
function insert_sccs(
        var_sccs::Vector{Vector{Int}}, sccs_to_insert::Vector{Tuple{Int, Vector{Int}}})
    # insert the new SCCs, accounting for the fact that we might have multiple entries
    # in `sccs_to_insert` to be inserted at the same index.
    old_idx = 1
    insert_idx = 1
    new_sccs = similar(var_sccs, length(var_sccs) + length(sccs_to_insert))
    for i in eachindex(new_sccs)
        # if we have SCCs to insert, and the index we have to insert them at is the current
        # one in the old list of SCCs
        if insert_idx <= length(sccs_to_insert) && sccs_to_insert[insert_idx][1] == old_idx
            # insert it
            new_sccs[i] = sccs_to_insert[insert_idx][2]
            insert_idx += 1
        else
            # otherwise, insert the old SCC
            new_sccs[i] = copy(var_sccs[old_idx])
            old_idx += 1
        end
    end

    filter!(!isempty, new_sccs)
    return new_sccs
end

"""
Check if there's `D(x) ~ x_t` already.
"""
function find_duplicate_dd(dv, solvable_graph, diff_to_var, linear_eqs, mm)
    for eq in 𝑑neighbors(solvable_graph, dv)
        mi = get(linear_eqs, eq, 0)
        iszero(mi) && continue
        row = @view mm[mi, :]
        nzs = nonzeros(row)
        rvs = SparseArrays.nonzeroinds(row)
        # note that `v_t` must not be differentiated
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

"""
Add a dummy derivative variable x_t corresponding to symbolic variable D(x)
which has index dv in `fullvars`. Return the new index of x_t.
"""
function add_dd_variable!(s::SystemStructure, fullvars, x_t, dv)
    push!(fullvars, simplify_shifts(x_t))
    v_t = length(fullvars)
    v_t_idx = add_vertex!(s.var_to_diff)
    add_vertex!(s.graph, DST)
    # TODO: do we care about solvable_graph? We don't use them after
    # `dummy_derivative_graph`.
    add_vertex!(s.solvable_graph, DST)
    s.var_to_diff[v_t] = s.var_to_diff[dv]
    v_t
end

"""
Add the equation D(x) - x_t ~ 0 to `neweqs`. `dv` and `v_t` are the indices
of the higher-order derivative variable and the newly-introduced dummy
derivative variable. Return the index of the new equation in `neweqs`.
"""
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

"""
Solve the equations in `neweqs` to obtain the final equations of the
system.

For each equation of `neweqs`, do one of the following:
   1. If the equation is solvable for a differentiated variable D(x),
      then solve for D(x), and add D(x) ~ sol as a differential equation
      of the system.
   2. If the equation is solvable for an un-differentiated variable x,
      solve for x and then add x ~ sol as a solved equation. These will
      become observables.
   3. If the equation is not solvable, add it as an algebraic equation.

Solved equations are added to `total_sub`. Occurrences of differential
or solved variables on the RHS of the final equations will get substituted.
The topological sort of the equations ensures that variables are solved for
before they appear in equations.

Reorder the equations and unknowns to be in the BLT sorted form.

Return the new equations, the solved equations,
the new orderings, and the number of solved variables and equations.
"""
function generate_system_equations!(state::TearingState, neweqs, var_eq_matching,
        full_var_eq_matching, var_sccs, extra_eqs_vars;
        simplify = false, iv = nothing, D = nothing)
    @unpack fullvars, sys, structure = state
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    eq_var_matching = invview(var_eq_matching)
    full_eq_var_matching = invview(full_var_eq_matching)
    diff_to_var = invview(var_to_diff)
    extra_eqs, extra_vars = extra_eqs_vars

    total_sub = Dict()
    if is_only_discrete(structure)
        for (i, v) in enumerate(fullvars)
            op = operation(v)
            op isa Shift && (op.steps < 0) &&
                begin
                    lowered = lower_shift_varname_with_unit(v, iv)
                    total_sub[v] = lowered
                    fullvars[i] = lowered
                end
        end
    end

    eq_generator = EquationGenerator(state, total_sub, D, iv)

    # We need to solve extra equations before everything to repsect
    # topological order.
    for eq in extra_eqs
        var = eq_var_matching[eq]
        var isa Int || continue
        codegen_equation!(eq_generator, neweqs[eq], eq, var; simplify)
    end

    # if the variable is present in the equations either as-is or differentiated
    ispresent = let var_to_diff = var_to_diff, graph = graph
        i -> (!isempty(𝑑neighbors(graph, i)) ||
              (var_to_diff[i] !== nothing && !isempty(𝑑neighbors(graph, var_to_diff[i]))))
    end

    digraph = DiCMOBiGraph{false}(graph, var_eq_matching)
    idep = iv
    for (i, scc) in enumerate(var_sccs)
        # note that the `vscc <-> escc` relation is a set-to-set mapping, and not
        # point-to-point.
        vscc, escc = get_sorted_scc(digraph, full_var_eq_matching, var_eq_matching, scc)
        var_sccs[i] = vscc

        if length(escc) != length(vscc)
            isempty(escc) && continue
            escc = setdiff(escc, extra_eqs)
            isempty(escc) && continue
            vscc = setdiff(vscc, extra_vars)
            isempty(vscc) && continue
        end

        offset = 1
        for ieq in escc
            iv = eq_var_matching[ieq]
            eq = neweqs[ieq]
            codegen_equation!(eq_generator, neweqs[ieq], ieq, iv; simplify)
        end
    end

    for eq in extra_eqs
        var = eq_var_matching[eq]
        var isa Int && continue
        codegen_equation!(eq_generator, neweqs[eq], eq, var; simplify)
    end

    @unpack neweqs′, eq_ordering, var_ordering, solved_eqs, solved_vars = eq_generator

    is_diff_eq = .!iszero.(var_ordering)
    # Generate new equations and orderings
    diff_vars = var_ordering[is_diff_eq]
    diff_vars_set = BitSet(diff_vars)
    if length(diff_vars_set) != length(diff_vars)
        error("Tearing internal error: lowering DAE into semi-implicit ODE failed!")
    end
    solved_vars_set = BitSet(solved_vars)
    # We filled zeros for algebraic variables, so fill them properly here
    offset = 1
    for (i, v) in enumerate(var_ordering)
        v == 0 || continue
        # find the next variable which is not differential or solved, is not the
        # derivative of another variable and is present in the equations
        index = findnext(1:ndsts(graph), offset) do j
            !(j in diff_vars_set || j in solved_vars_set) && diff_to_var[j] === nothing &&
                ispresent(j)
        end
        # in case of overdetermined systems, this may not be present
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

"""
    $(TYPEDSIGNATURES)

Sort the provided SCC `scc`, given the `digraph` of the system constructed using
`var_eq_matching` along with both the matchings of the system.
"""
function get_sorted_scc(
        digraph::DiCMOBiGraph, full_var_eq_matching::Matching, var_eq_matching::Matching, scc::Vector{Int})
    eq_var_matching = invview(var_eq_matching)
    full_eq_var_matching = invview(full_var_eq_matching)
    # obtain the matched equations in the SCC
    scc_eqs = Int[full_var_eq_matching[v] for v in scc if full_var_eq_matching[v] isa Int]
    # obtain the equations in the SCC that are linearly solvable
    scc_solved_eqs = Int[var_eq_matching[v] for v in scc if var_eq_matching[v] isa Int]
    # obtain the subgraph of the contracted graph involving the solved equations
    subgraph, varmap = Graphs.induced_subgraph(digraph, scc_solved_eqs)
    # topologically sort the solved equations and append the remainder
    scc_eqs = [varmap[reverse(topological_sort(subgraph))];
               setdiff(scc_eqs, scc_solved_eqs)]
    # the variables of the SCC are obtained by inverse mapping the sorted equations
    # and appending the rest
    scc_vars = [eq_var_matching[e] for e in scc_eqs if eq_var_matching[e] isa Int]
    append!(scc_vars, setdiff(scc, scc_vars))
    return scc_vars, scc_eqs
end

"""
    $(TYPEDSIGNATURES)

Struct containing the information required to generate equations of a system, as well as
the generated equations and associated metadata.
"""
struct EquationGenerator{S, D, I}
    """
    `TearingState` of the system.
    """
    state::S
    """
    Substitutions to perform in all subsequent equations. For each differential equation
    `D(x) ~ f(..)`, the substitution `D(x) => f(..)` is added to the rules.
    """
    total_sub::Dict{Any, Any}
    """
    The differential operator, or `nothing` if not applicable.
    """
    D::D
    """
    The independent variable, or `nothing` if not applicable.
    """
    idep::I
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

"""
    $(TYPEDSIGNATURES)

Check if equation at index `ieq` is linearly solvable for variable at index `iv`.
"""
function is_solvable(eg::EquationGenerator, ieq, iv)
    solvable_graph = eg.state.structure.solvable_graph
    return ieq isa Int && iv isa Int && BipartiteEdge(ieq, iv) in solvable_graph
end

"""
    $(TYPEDSIGNATURES)

    If `iv` is like D(x) or Shift(t, 1)(x)
"""
function is_dervar(eg::EquationGenerator, iv::Int)
    diff_to_var = invview(eg.state.structure.var_to_diff)
    diff_to_var[iv] !== nothing
end

"""
    $(TYPEDSIGNATURES)

Appropriately codegen the given equation `eq`, which occurs at index `ieq` in the untorn
list of equations and is matched to variable at index `iv`.
"""
function codegen_equation!(eg::EquationGenerator,
        eq::Equation, ieq::Int, iv::Union{Int, Unassigned}; simplify = false)
    # We generate equations ordered by the matched variables
    #   Solvable equations of differential variables D(x) become differential equations
    #   Solvable equations of non-differential variables become observable equations
    #   Non-solvable equations become algebraic equations.
    @unpack state, total_sub, neweqs′, eq_ordering, var_ordering = eg
    @unpack solved_eqs, solved_vars, D, idep = eg
    @unpack fullvars, sys, structure = state
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    diff_to_var = invview(var_to_diff)

    issolvable = is_solvable(eg, ieq, iv)
    isdervar = issolvable && is_dervar(eg, iv)
    isdisc = is_only_discrete(structure)
    # The variable is derivative variable and the "most differentiated"
    # This is only used for discrete systems, and basically refers to
    # `Shift(t, 1)(x(k))` in `Shift(t, 1)(x(k)) ~ x(k) + x(k-1)`. As illustrated in
    # the docstring for `add_additional_history!`, this is an exception and needs to be
    # treated like a solved equation rather than a differential equation.
    is_highest_diff = iv isa Int && isdervar && var_to_diff[iv] === nothing
    if issolvable && isdervar && (!isdisc || !is_highest_diff)
        var = fullvars[iv]
        isnothing(D) && throw(UnexpectedDifferentialError(equations(sys)[ieq]))
        order, lv = var_order(iv, diff_to_var)
        dx = D(simplify_shifts(fullvars[lv]))

        neweq = make_differential_equation(var, dx, eq, total_sub)
        # We will add `neweq.lhs` to `total_sub`, so any equation involving it won't be
        # incident on it. Remove the edges incident on `iv` from the graph, and add
        # the replacement vertices from `ieq` so that the incidence is still correct.
        for e in 𝑑neighbors(graph, iv)
            e == ieq && continue
            for v in 𝑠neighbors(graph, ieq)
                add_edge!(graph, e, v)
            end
            rem_edge!(graph, e, iv)
        end

        total_sub[simplify_shifts(neweq.lhs)] = neweq.rhs
        # Substitute unshifted variables x(k), y(k) on RHS of implicit equations
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
            # backshift solved equations to calculate the value of the variable at the
            # current time. This works because we added one additional history element
            # in `add_additional_history!`.
            if isdisc
                neweq = backshift_expr(neweq, idep)
            end
            push!(solved_eqs, neweq)
            push!(solved_vars, iv)
        end
    else
        neweq = make_algebraic_equation(eq, total_sub)
        # For the same reason as solved equations (they are effectively the same)
        if isdisc
            neweq = backshift_expr(neweq, idep)
        end
        push!(neweqs′, neweq)
        push!(eq_ordering, ieq)
        # we push a dummy to `var_ordering` here because `iv` is `unassigned`
        push!(var_ordering, 0)
    end
end

"""
Occurs when a variable D(x) occurs in a non-differential system.
"""
struct UnexpectedDifferentialError
    eq::Equation
end

function Base.showerror(io::IO, err::UnexpectedDifferentialError)
    error("Differential found in a non-differential system. Likely this is a bug in the construction of an initialization system. Please report this issue with a reproducible example. Offending equation: $(err.eq)")
end

"""
Generate a first-order differential equation whose LHS is `dx`.

`var` and `dx` represent the same variable, but `var` may be a higher-order differential and `dx` is always first-order. For example, if `var` is D(D(x)), then `dx` would be `D(x_t)`. Solve `eq` for `var`, substitute previously solved variables, and return the differential equation.
"""
function make_differential_equation(var, dx, eq, total_sub)
    dx ~ simplify_shifts(Symbolics.fixpoint_sub(
        Symbolics.symbolic_linear_solve(eq, var),
        total_sub; operator = ModelingToolkit.Shift))
end

"""
Generate an algebraic equation. Substitute solved variables into `eq` and return the equation.
"""
function make_algebraic_equation(eq, total_sub)
    rhs = eq.rhs
    if !(eq.lhs isa Number && eq.lhs == 0)
        rhs = eq.rhs - eq.lhs
    end
    0 ~ simplify_shifts(Symbolics.fixpoint_sub(rhs, total_sub))
end

"""
Solve equation `eq` for `var`, substitute previously solved variables, and return the solved equation.
"""
function make_solved_equation(var, eq, total_sub; simplify = false)
    residual = eq.lhs - eq.rhs
    a, b, islinear = linear_expansion(residual, var)
    @assert islinear
    # 0 ~ a * var + b
    # var ~ -b/a
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

"""
Given the ordering returned by `generate_system_equations!`, update the
tearing state to account for the new order. Permute the variables and equations.
Eliminate the solved variables and equations from the graph and permute the
graph's vertices to account for the new variable/equation ordering.
"""
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

    # Contract the vertices in the structure graph to make the structure match
    # the new reality of the system we've just created.
    new_graph = contract_variables(graph, var_eq_matching, varsperm, eqsperm,
        nsolved_eq, nsolved_var)
    new_solvable_graph = contract_variables(solvable_graph, var_eq_matching, varsperm, eqsperm,
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

    # Update the SCCs
    var_ordering_set = BitSet(var_ordering)
    for scc in var_sccs
        # Map variables to their new indices
        map!(v -> varsperm[v], scc, scc)
        # Remove variables not in the reduced set
        filter!(!iszero, scc)
    end
    # Remove empty SCCs
    filter!(!isempty, var_sccs)

    # Update system structure
    @set! state.structure.graph = complete(new_graph)
    @set! state.structure.solvable_graph = complete(new_solvable_graph)
    @set! state.structure.var_to_diff = new_var_to_diff
    @set! state.structure.eq_to_diff = new_eq_to_diff
    @set! state.fullvars = new_fullvars
    state
end

"""
Update the system equations, unknowns, and observables after simplification.
"""
function update_simplified_system!(
        state::TearingState, neweqs, solved_eqs, dummy_sub, var_sccs, extra_unknowns;
        array_hack = true, D = nothing, iv = nothing)
    @unpack fullvars, structure = state
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    diff_to_var = invview(var_to_diff)
    # Since we solved the highest order derivative variable in discrete systems,
    # we make a list of the solved variables and avoid including them in the
    # unknowns.
    solved_vars = Set()
    if is_only_discrete(structure)
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
    # TODO: compute the dependency correctly so that we don't have to do this
    obs = [fast_substitute(observed(sys), obs_sub); solved_eqs;
           fast_substitute(state.additional_observed, obs_sub)]

    unknown_idxs = filter(
        i -> diff_to_var[i] === nothing && ispresent(i) && !(fullvars[i] in solved_vars), eachindex(state.fullvars))
    unknowns = state.fullvars[unknown_idxs]
    unknowns = [unknowns; extra_unknowns]
    if is_only_discrete(structure)
        # Algebraic variables are shifted forward by one, so we backshift them.
        unknowns = map(enumerate(unknowns)) do (i, var)
            if iscall(var) && operation(var) isa Shift && operation(var).steps == 1
                # We might have shifted a variable with io metadata. That is irrelevant now
                # because we handled io variables earlier in `_mtkcompile!` so just ignore
                # it here.
                setio(backshift_expr(var, iv), false, false)
            else
                var
            end
        end
    end
    @set! sys.unknowns = unknowns

    obs = tearing_hacks(sys, obs, unknowns, neweqs; array = array_hack)

    @set! sys.eqs = neweqs
    @set! sys.observed = obs

    # Only makes sense for time-dependent
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

"""
Give the order of the variable indexed by dv.
"""
function var_order(dv, diff_to_var)
    order = 0
    while (dv′ = diff_to_var[dv]) !== nothing
        order += 1
        dv = dv′
    end
    order, dv
end

"""
Main internal function for structural simplification for DAE systems and discrete systems.
Generate dummy derivative variables, new equations in terms of variables, return updated
system and tearing state.

Terminology and Definition:

A general DAE is in the form of `F(u'(t), u(t), p, t) == 0`. We can
characterize variables in `u(t)` into two classes: differential variables
(denoted `v(t)`) and algebraic variables (denoted `z(t)`). Differential
variables are marked as `SelectedState` and they are differentiated in the
DAE system, i.e. `v'(t)` are all the variables in `u'(t)` that actually
appear in the system. Algebraic variables are variables that are not
differential variables.

# Arguments

- `state`: The `TearingState` of the system.
- `var_eq_matching`: The maximal matching after state selection.
- `full_var_eq_matching`: The maximal matching prior to state selection.
- `var_sccs`: The topologically sorted strongly connected components of the system
  according to `full_var_eq_matching`.
"""
@kwdef struct DefaultReassembleAlgorithm <: ReassembleAlgorithm
    simplify::Bool = false
    array_hack::Bool = true
end

function (alg::DefaultReassembleAlgorithm)(state::TearingState, tearing_result::TearingResult, mm::Union{SparseMatrixCLIL,  Nothing}; fully_determined::Bool = true, kw...)
    @unpack simplify, array_hack = alg
    @unpack var_eq_matching, full_var_eq_matching, var_sccs = tearing_result

    extra_eqs_vars = get_extra_eqs_vars(
        state, var_eq_matching, full_var_eq_matching, fully_determined)
    neweqs = collect(equations(state))
    dummy_sub = Dict()

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
            state, neweqs, var_eq_matching, full_var_eq_matching, var_sccs; iv, D)
    end

    # Structural simplification
    substitute_derivatives_algevars!(state, neweqs, var_eq_matching, dummy_sub; iv, D)

    var_sccs = generate_derivative_variables!(
        state, neweqs, var_eq_matching, full_var_eq_matching, var_sccs; mm, iv, D)

    neweqs, solved_eqs,
    eq_ordering,
    var_ordering,
    nelim_eq,
    nelim_var = generate_system_equations!(
        state, neweqs, var_eq_matching, full_var_eq_matching,
        var_sccs, extra_eqs_vars; simplify, iv, D)

    state = reorder_vars!(
        state, var_eq_matching, var_sccs, eq_ordering, var_ordering, nelim_eq, nelim_var)
    # var_eq_matching and full_var_eq_matching are now invalidated

    sys = update_simplified_system!(state, neweqs, solved_eqs, dummy_sub, var_sccs,
        extra_unknowns; array_hack, iv, D)

    @set! state.sys = sys
    @set! sys.tearing_state = state
    return invalidate_cache!(sys)
end

"""
    $(TYPEDSIGNATURES)

Add one more history equation for discrete systems. For example, if we have

```julia
Shift(t, 1)(x(k-1)) ~ x(k)
Shift(t, 1)(x(k)) ~ x(k) + x(k-1)
```

This turns it into

```julia
Shift(t, 1)(x(k-2)) ~ x(k-1)
Shift(t, 1)(x(k-1)) ~ x(k)
Shift(t, 1)(x(k)) ~ x(k) + x(k-1)
```

Thus adding an additional unknown as well. Later, the highest derivative equation will
be backshifted by one and turned into an observed equation, resulting in:

```julia
Shift(t, 1)(x(k-2)) ~ x(k-1)
Shift(t, 1)(x(k-1)) ~ x(k)

x(k) ~ x(k-1) + x(k-2)
```

Where the last equation is the observed equation.
"""
function add_additional_history!(
        state::TearingState, neweqs::Vector, var_eq_matching::Matching,
        full_var_eq_matching::Matching, var_sccs::Vector{Vector{Int}}; iv, D)
    @unpack fullvars, sys, structure = state
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    eq_var_matching = invview(var_eq_matching)
    diff_to_var = invview(var_to_diff)
    is_discrete = is_only_discrete(structure)
    digraph = DiCMOBiGraph{false}(graph, var_eq_matching)

    # We need the inverse mapping of `var_sccs` to update it efficiently later.
    v_to_scc = Vector{NTuple{2, Int}}(undef, ndsts(graph))
    for (i, scc) in enumerate(var_sccs), (j, v) in enumerate(scc)

        v_to_scc[v] = (i, j)
    end

    vars_to_backshift = BitSet()
    eqs_to_backshift = BitSet()
    # add history for differential variables
    for ivar in 1:length(fullvars)
        ieq = var_eq_matching[ivar]
        # the variable to backshift is a state variable which is not the
        # derivative of any other one.
        ieq isa SelectedState || continue
        diff_to_var[ivar] === nothing || continue
        push!(vars_to_backshift, ivar)
    end

    inserts = Tuple{Int, Vector{Int}}[]

    for var in vars_to_backshift
        add_backshifted_var!(state, var, iv)
        # all backshifted vars are differential vars, hence SelectedState
        push!(var_eq_matching, SelectedState())
        push!(full_var_eq_matching, unassigned)
        # add to the SCCs right before the variable that was backshifted
        push!(inserts, (v_to_scc[var][1], [length(fullvars)]))
    end

    sort!(inserts, by = first)
    new_sccs = insert_sccs(var_sccs, inserts)
    return new_sccs
end

"""
    $(TYPEDSIGNATURES)

Add the backshifted version of variable `ivar` to the system.
"""
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

"""
    $(TYPEDSIGNATURES)

Backshift the given expression `ex`.
"""
function backshift_expr(ex, iv)
    ex isa Symbolic || return ex
    return descend_lower_shift_varname_with_unit(
        simplify_shifts(distribute_shift(Shift(iv, -1)(ex))), iv)
end

function backshift_expr(ex::Equation, iv)
    return backshift_expr(ex.lhs, iv) ~ backshift_expr(ex.rhs, iv)
end

"""
    $(TYPEDSIGNATURES)

Return a 2-tuple of integer vectors containing indices of extra equations and variables
respectively. For fully-determined systems, both of these are empty. Overdetermined systems
have extra equations, and underdetermined systems have extra variables.
"""
function get_extra_eqs_vars(
        state::TearingState, var_eq_matching::Matching, full_var_eq_matching::Matching, fully_determined::Bool)
    fully_determined && return Int[], Int[]

    extra_eqs = Int[]
    extra_vars = Int[]
    full_eq_var_matching = invview(full_var_eq_matching)

    for v in 𝑑vertices(state.structure.graph)
        eq = full_var_eq_matching[v]
        eq isa Int && continue
        # Only if the variable is also unmatched in `var_eq_matching`.
        # Otherwise, `SelectedState` differential variables from order lowering
        # are also considered "extra"
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

"""
# HACK

Add equations for array observed variables. If `p[i] ~ (...)` are equations, add an
equation `p ~ [p[1], p[2], ...]` allow topsort to reorder them only add the new equation
if all `p[i]` are present and the unscalarized form is used in any equation (observed or
not) we first count the number of times the scalarized form of each observed variable
occurs in observed equations (and unknowns if it's split).
"""
function tearing_hacks(sys, obs, unknowns, neweqs; array = true)
    # map of array observed variable (unscalarized) to number of its
    # scalarized terms that appear in observed equations
    arr_obs_occurrences = Dict()
    for (i, eq) in enumerate(obs)
        lhs = eq.lhs
        rhs = eq.rhs

        array || continue
        iscall(lhs) || continue
        operation(lhs) === getindex || continue
        Symbolics.shape(lhs) != Symbolics.Unknown() || continue
        arg1 = arguments(lhs)[1]
        cnt = get(arr_obs_occurrences, arg1, 0)
        arr_obs_occurrences[arg1] = cnt + 1
        continue
    end

    # count variables in unknowns if they are scalarized forms of variables
    # also present as observed. e.g. if `x[1]` is an unknown and `x[2] ~ (..)`
    # is an observed equation.
    for sym in unknowns
        iscall(sym) || continue
        operation(sym) === getindex || continue
        Symbolics.shape(sym) != Symbolics.Unknown() || continue
        arg1 = arguments(sym)[1]
        cnt = get(arr_obs_occurrences, arg1, 0)
        cnt == 0 && continue
        arr_obs_occurrences[arg1] = cnt + 1
    end

    obs_arr_eqs = Equation[]
    for (arrvar, cnt) in arr_obs_occurrences
        cnt == length(arrvar) || continue
        # firstindex returns 1 for multidimensional array symbolics
        firstind = Tuple(first(eachindex(arrvar)))
        scal = [arrvar[i] for i in eachindex(arrvar)]
        # respect non-1-indexed arrays
        # TODO: get rid of this hack together with the above hack, then remove OffsetArrays dependency
        # `change_origin` is required because `Origin(firstind)(scal)` makes codegen
        # try to `create_array(OffsetArray{...}, ...)` which errors.
        # `term(Origin(firstind), scal)` doesn't retain the `symtype` and `size`
        # of `scal`.
        rhs = change_origin(firstind, scal)
        push!(obs_arr_eqs, arrvar ~ rhs)
    end
    append!(obs, obs_arr_eqs)

    return obs
end

# PART OF HACK
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

function tearing(state::TearingState; tearing_alg::TearingAlgorithm = DummyDerivativeTearing(),
                 kwargs...)
    state.structure.solvable_graph === nothing && find_solvables!(state; kwargs...)
    complete!(state.structure)
    tearing_alg(state.structure)
end

"""
    tearing(sys)

Tear the nonlinear equations in system. When `simplify=true`, we simplify the
new residual equations after tearing. End users are encouraged to call [`mtkcompile`](@ref)
instead, which calls this function internally.
"""
function tearing(sys::AbstractSystem, state = TearingState(sys); mm = nothing,
        reassemble_alg::ReassembleAlgorithm = DefaultReassembleAlgorithm(),
        fully_determined = true, kwargs...)
    tearing_result, extras = tearing(state; kwargs...)
    invalidate_cache!(reassemble_alg(state, tearing_result, mm; fully_determined))
end

"""
    dummy_derivative(sys)

Perform index reduction and use the dummy derivative technique to ensure that
the system is balanced.
"""
function dummy_derivative(sys, state = TearingState(sys);
        reassemble_alg::ReassembleAlgorithm = DefaultReassembleAlgorithm(),
        mm = nothing, fully_determined = true, kwargs...)
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
    tearing_result, extras = dummy_derivative_graph!(
        state, jac; state_priority, kwargs...)
    reassemble_alg(state, tearing_result, mm; fully_determined)
end
