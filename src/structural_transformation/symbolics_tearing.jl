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

function var_derivative!(ts::TearingState{ODESystem}, v::Int)
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

function eq_derivative!(ts::TearingState{ODESystem}, ieq::Int; kwargs...)
    s = ts.structure

    eq_diff = eq_derivative_graph!(s, ieq)

    sys = ts.sys
    eq = equations(ts)[ieq]
    eq = 0 ~ ModelingToolkit.derivative(eq.rhs - eq.lhs, get_iv(sys))
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

function tearing_sub(expr, dict, s)
    expr = Symbolics.fixpoint_sub(expr, dict)
    s ? simplify(expr) : expr
end

"""
$(TYPEDSIGNATURES)

Like `equations(sys)`, but includes substitutions done by the tearing process.
These equations matches generated numerical code.

See also [`equations`](@ref) and [`ModelingToolkit.get_eqs`](@ref).
"""
function full_equations(sys::AbstractSystem; simplify = false)
    empty_substitutions(sys) && return equations(sys)
    substitutions = get_substitutions(sys)
    substitutions.subed_eqs === nothing || return substitutions.subed_eqs
    @unpack subs = substitutions
    solved = Dict(eq.lhs => eq.rhs for eq in subs)
    neweqs = map(equations(sys)) do eq
        if iscall(eq.lhs) && operation(eq.lhs) isa Union{Shift, Differential}
            return tearing_sub(eq.lhs, solved, simplify) ~ tearing_sub(eq.rhs, solved,
                simplify)
        else
            if !(eq.lhs isa Number && eq.lhs == 0)
                eq = 0 ~ eq.rhs - eq.lhs
            end
            rhs = tearing_sub(eq.rhs, solved, simplify)
            if rhs isa Symbolic
                return 0 ~ rhs
            else # a number
                error("tearing failed because the system is singular")
            end
        end
        eq
    end
    substitutions.subed_eqs = neweqs
    return neweqs
end

function tearing_substitution(sys::AbstractSystem; kwargs...)
    neweqs = full_equations(sys::AbstractSystem; kwargs...)
    @set! sys.eqs = neweqs
    @set! sys.substitutions = nothing
end

function tearing_assignments(sys::AbstractSystem)
    if empty_substitutions(sys)
        assignments = []
        deps = Int[]
        sol_states = Code.LazyState()
    else
        @unpack subs, deps = get_substitutions(sys)
        assignments = [Assignment(eq.lhs, eq.rhs) for eq in subs]
        sol_states = Code.NameState(Dict(eq.lhs => Symbol(eq.lhs) for eq in subs))
    end
    return assignments, deps, sol_states
end

function solve_equation(eq, var, simplify)
    rhs = value(solve_for(eq, var; simplify = simplify, check = false))
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

function tearing_reassemble(state::TearingState, var_eq_matching,
        full_var_eq_matching = nothing; simplify = false, mm = nothing)
    @unpack fullvars, sys, structure = state
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    extra_vars = Int[]
    if full_var_eq_matching !== nothing
        for v in 𝑑vertices(state.structure.graph)
            eq = full_var_eq_matching[v]
            eq isa Int && continue
            push!(extra_vars, v)
        end
    end

    neweqs = collect(equations(state))
    # Terminology and Definition:
    #
    # A general DAE is in the form of `F(u'(t), u(t), p, t) == 0`. We can
    # characterize variables in `u(t)` into two classes: differential variables
    # (denoted `v(t)`) and algebraic variables (denoted `z(t)`). Differential
    # variables are marked as `SelectedState` and they are differentiated in the
    # DAE system, i.e. `v'(t)` are all the variables in `u'(t)` that actually
    # appear in the system. Algebraic variables are variables that are not
    # differential variables.
    #
    # Dummy derivatives may determine that some differential variables are
    # algebraic variables in disguise. The derivative of such variables are
    # called dummy derivatives.

    # Step 1:
    # Replace derivatives of non-selected unknown variables by dummy derivatives

    if ModelingToolkit.has_iv(state.sys)
        iv = get_iv(state.sys)
        if is_only_discrete(state.structure)
            D = Shift(iv, 1)
        else
            D = Differential(iv)
        end
    else
        iv = D = nothing
    end
    diff_to_var = invview(var_to_diff)
    dummy_sub = Dict()
    for var in 1:length(fullvars)
        dv = var_to_diff[var]
        dv === nothing && continue
        if var_eq_matching[var] !== SelectedState()
            dd = fullvars[dv]
            v_t = setio(diff2term(unwrap(dd)), false, false)
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

    # `SelectedState` information is no longer needed past here. State selection
    # is done. All non-differentiated variables are algebraic variables, and all
    # variables that appear differentiated are differential variables.

    ### extract partition information
    is_solvable = let solvable_graph = solvable_graph
        (eq, iv) -> eq isa Int && iv isa Int && BipartiteEdge(eq, iv) in solvable_graph
    end

    # if var is like D(x)
    isdervar = let diff_to_var = diff_to_var
        var -> diff_to_var[var] !== nothing
    end
    var_order = let diff_to_var = diff_to_var
        dv -> begin
            order = 0
            while (dv′ = diff_to_var[dv]) !== nothing
                order += 1
                dv = dv′
            end
            order, dv
        end
    end

    #retear = BitSet()
    # There are three cases where we want to generate new variables to convert
    # the system into first order (semi-implicit) ODEs.
    #
    # 1. To first order:
    # Whenever higher order differentiated variable like `D(D(D(x)))` appears,
    # we introduce new variables `x_t`, `x_tt`, and `x_ttt` and new equations
    # ```
    # D(x_tt) = x_ttt
    # D(x_t) = x_tt
    # D(x) = x_t
    # ```
    # and replace `D(x)` to `x_t`, `D(D(x))` to `x_tt`, and `D(D(D(x)))` to
    # `x_ttt`.
    #
    # 2. To implicit to semi-implicit ODEs:
    # 2.1: Unsolvable derivative:
    # If one derivative variable `D(x)` is unsolvable in all the equations it
    # appears in, then we introduce a new variable `x_t`, a new equation
    # ```
    # D(x) ~ x_t
    # ```
    # and replace all other `D(x)` to `x_t`.
    #
    # 2.2: Solvable derivative:
    # If one derivative variable `D(x)` is solvable in at least one of the
    # equations it appears in, then we introduce a new variable `x_t`. One of
    # the solvable equations must be in the form of `0 ~ L(D(x), u...)` and
    # there exists a function `l` such that `D(x) ~ l(u...)`. We should replace
    # it to
    # ```
    # 0 ~ x_t - l(u...)
    # D(x) ~ x_t
    # ```
    # and replace all other `D(x)` to `x_t`.
    #
    # Observe that we don't need to actually introduce a new variable `x_t`, as
    # the above equations can be lowered to
    # ```
    # x_t := l(u...)
    # D(x) ~ x_t
    # ```
    # where `:=` denotes assignment.
    #
    # As a final note, in all the above cases where we need to introduce new
    # variables and equations, don't add them when they already exist.
    eq_var_matching = invview(var_eq_matching)
    linear_eqs = mm === nothing ? Dict{Int, Int}() :
                 Dict(reverse(en) for en in enumerate(mm.nzrows))
    for v in 1:length(var_to_diff)
        dv = var_to_diff[v]
        dv isa Int || continue
        solved = var_eq_matching[dv] isa Int
        solved && continue
        # check if there's `D(x) = x_t` already
        local v_t, dummy_eq
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
                dummy_eq = eq
                @goto FOUND_DUMMY_EQ
            end
        end
        dx = fullvars[dv]
        # add `x_t`
        order, lv = var_order(dv)
        x_t = lower_varname_withshift(fullvars[lv], iv, order)
        push!(fullvars, simplify_shifts(x_t))
        v_t = length(fullvars)
        v_t_idx = add_vertex!(var_to_diff)
        add_vertex!(graph, DST)
        # TODO: do we care about solvable_graph? We don't use them after
        # `dummy_derivative_graph`.
        add_vertex!(solvable_graph, DST)
        # var_eq_matching is a bit odd.
        # length(var_eq_matching) == length(invview(var_eq_matching))
        push!(var_eq_matching, unassigned)
        @assert v_t_idx == ndsts(graph) == ndsts(solvable_graph) == length(fullvars) ==
                length(var_eq_matching)
        # add `D(x) - x_t ~ 0`
        push!(neweqs, 0 ~ x_t - dx)
        add_vertex!(graph, SRC)
        dummy_eq = length(neweqs)
        add_edge!(graph, dummy_eq, dv)
        add_edge!(graph, dummy_eq, v_t)
        add_vertex!(solvable_graph, SRC)
        add_edge!(solvable_graph, dummy_eq, dv)
        @assert nsrcs(graph) == nsrcs(solvable_graph) == dummy_eq
        @label FOUND_DUMMY_EQ
        var_to_diff[v_t] = var_to_diff[dv]
        var_eq_matching[dv] = unassigned
        eq_var_matching[dummy_eq] = dv
    end

    # Will reorder equations and unknowns to be:
    # [diffeqs; ...]
    # [diffvars; ...]
    # such that the mass matrix is:
    # [I  0
    #  0  0].
    diffeq_idxs = Int[]
    algeeq_idxs = Int[]
    diff_eqs = Equation[]
    alge_eqs = Equation[]
    diff_vars = Int[]
    subeqs = Equation[]
    solved_equations = Int[]
    solved_variables = Int[]
    # Solve solvable equations
    toporder = topological_sort(DiCMOBiGraph{false}(graph, var_eq_matching))
    eqs = Iterators.reverse(toporder)
    total_sub = Dict()
    idep = iv
    for ieq in eqs
        iv = eq_var_matching[ieq]
        if is_solvable(ieq, iv)
            # We don't solve differential equations, but we will need to try to
            # convert it into the mass matrix form.
            # We cannot solve the differential variable like D(x)
            if isdervar(iv)
                order, lv = var_order(iv)
                dx = D(simplify_shifts(lower_varname_withshift(
                    fullvars[lv], idep, order - 1)))
                eq = dx ~ simplify_shifts(Symbolics.fixpoint_sub(
                    Symbolics.solve_for(neweqs[ieq],
                        fullvars[iv]),
                    total_sub; operator = ModelingToolkit.Shift))
                for e in 𝑑neighbors(graph, iv)
                    e == ieq && continue
                    for v in 𝑠neighbors(graph, e)
                        add_edge!(graph, e, v)
                    end
                    rem_edge!(graph, e, iv)
                end
                push!(diff_eqs, eq)
                total_sub[simplify_shifts(eq.lhs)] = eq.rhs
                push!(diffeq_idxs, ieq)
                push!(diff_vars, diff_to_var[iv])
                continue
            end
            eq = neweqs[ieq]
            var = fullvars[iv]
            residual = eq.lhs - eq.rhs
            a, b, islinear = linear_expansion(residual, var)
            @assert islinear
            # 0 ~ a * var + b
            # var ~ -b/a
            if ModelingToolkit._iszero(a)
                @warn "Tearing: solving $eq for $var is singular!"
            else
                rhs = -b / a
                neweq = var ~ Symbolics.fixpoint_sub(
                    simplify ?
                    Symbolics.simplify(rhs) : rhs,
                    total_sub; operator = ModelingToolkit.Shift)
                push!(subeqs, neweq)
                push!(solved_equations, ieq)
                push!(solved_variables, iv)
            end
        else
            eq = neweqs[ieq]
            rhs = eq.rhs
            if !(eq.lhs isa Number && eq.lhs == 0)
                rhs = eq.rhs - eq.lhs
            end
            push!(alge_eqs, 0 ~ Symbolics.fixpoint_sub(rhs, total_sub))
            push!(algeeq_idxs, ieq)
        end
    end
    # TODO: BLT sorting
    neweqs = [diff_eqs; alge_eqs]
    inveqsperm = [diffeq_idxs; algeeq_idxs]
    eqsperm = zeros(Int, nsrcs(graph))
    for (i, v) in enumerate(inveqsperm)
        eqsperm[v] = i
    end
    diff_vars_set = BitSet(diff_vars)
    if length(diff_vars_set) != length(diff_vars)
        error("Tearing internal error: lowering DAE into semi-implicit ODE failed!")
    end
    solved_variables_set = BitSet(solved_variables)
    invvarsperm = [diff_vars;
                   setdiff!(setdiff(1:ndsts(graph), diff_vars_set),
                       solved_variables_set)]
    varsperm = zeros(Int, ndsts(graph))
    for (i, v) in enumerate(invvarsperm)
        varsperm[v] = i
    end

    deps = Vector{Int}[i == 1 ? Int[] : collect(1:(i - 1))
                       for i in 1:length(solved_equations)]
    # Contract the vertices in the structure graph to make the structure match
    # the new reality of the system we've just created.
    graph = contract_variables(graph, var_eq_matching, varsperm, eqsperm,
        length(solved_variables), length(solved_variables_set))

    # Update system
    new_var_to_diff = complete(DiffGraph(length(invvarsperm)))
    for (v, d) in enumerate(var_to_diff)
        v′ = varsperm[v]
        (v′ > 0 && d !== nothing) || continue
        d′ = varsperm[d]
        new_var_to_diff[v′] = d′ > 0 ? d′ : nothing
    end
    new_eq_to_diff = complete(DiffGraph(length(inveqsperm)))
    for (v, d) in enumerate(eq_to_diff)
        v′ = eqsperm[v]
        (v′ > 0 && d !== nothing) || continue
        d′ = eqsperm[d]
        new_eq_to_diff[v′] = d′ > 0 ? d′ : nothing
    end

    var_to_diff = new_var_to_diff
    eq_to_diff = new_eq_to_diff
    diff_to_var = invview(var_to_diff)

    old_fullvars = fullvars
    @set! state.structure.graph = complete(graph)
    @set! state.structure.var_to_diff = var_to_diff
    @set! state.structure.eq_to_diff = eq_to_diff
    @set! state.fullvars = fullvars = fullvars[invvarsperm]
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
    obs = [fast_substitute(observed(sys), obs_sub); subeqs]

    # HACK: Substitute non-scalarized symbolic arrays of observed variables
    # E.g. if `p[1] ~ (...)` and `p[2] ~ (...)` then substitute `p => [p[1], p[2]]` in all equations
    # ideally, we want to support equations such as `p ~ [p[1], p[2]]` which will then be handled
    # by the topological sorting and dependency identification pieces
    obs_arr_subs = Dict()

    for eq in obs
        lhs = eq.lhs
        iscall(lhs) || continue
        operation(lhs) === getindex || continue
        Symbolics.shape(lhs) !== Symbolics.Unknown() || continue
        arg1 = arguments(lhs)[1]
        haskey(obs_arr_subs, arg1) && continue
        obs_arr_subs[arg1] = [arg1[i] for i in eachindex(arg1)]
    end
    for i in eachindex(neweqs)
        neweqs[i] = fast_substitute(neweqs[i], obs_arr_subs; operator = Symbolics.Operator)
    end
    for i in eachindex(obs)
        obs[i] = fast_substitute(obs[i], obs_arr_subs; operator = Symbolics.Operator)
    end
    for i in eachindex(subeqs)
        subeqs[i] = fast_substitute(subeqs[i], obs_arr_subs; operator = Symbolics.Operator)
    end

    @set! sys.eqs = neweqs
    @set! sys.observed = obs

    unknowns = Any[v
                   for (i, v) in enumerate(fullvars)
                   if diff_to_var[i] === nothing && ispresent(i)]
    if !isempty(extra_vars)
        for v in extra_vars
            push!(unknowns, old_fullvars[v])
        end
    end
    @set! sys.unknowns = unknowns
    @set! sys.substitutions = Substitutions(subeqs, deps)

    # Only makes sense for time-dependent
    # TODO: generalize to SDE
    if sys isa ODESystem
        @set! sys.schedule = Schedule(var_eq_matching, dummy_sub)
    end
    @set! state.sys = sys
    @set! sys.tearing_state = state
    return invalidate_cache!(sys)
end

function tearing(state::TearingState; kwargs...)
    state.structure.solvable_graph === nothing && find_solvables!(state; kwargs...)
    complete!(state.structure)
    tearing_with_dummy_derivatives(state.structure, ())
end

"""
    tearing(sys; simplify=false)

Tear the nonlinear equations in system. When `simplify=true`, we simplify the
new residual equations after tearing. End users are encouraged to call [`structural_simplify`](@ref)
instead, which calls this function internally.
"""
function tearing(sys::AbstractSystem, state = TearingState(sys); mm = nothing,
        simplify = false, kwargs...)
    var_eq_matching, full_var_eq_matching = tearing(state)
    invalidate_cache!(tearing_reassemble(
        state, var_eq_matching, full_var_eq_matching; mm, simplify))
end

"""
    partial_state_selection(sys; simplify=false)

Perform partial state selection and tearing.
"""
function partial_state_selection(sys; simplify = false)
    state = TearingState(sys)
    var_eq_matching = partial_state_selection_graph!(state)

    tearing_reassemble(state, var_eq_matching; simplify = simplify)
end

"""
    dummy_derivative(sys)

Perform index reduction and use the dummy derivative technique to ensure that
the system is balanced.
"""
function dummy_derivative(sys, state = TearingState(sys); simplify = false,
        mm = nothing, kwargs...)
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
    var_eq_matching = dummy_derivative_graph!(state, jac; state_priority,
        kwargs...)
    tearing_reassemble(state, var_eq_matching; simplify, mm)
end
