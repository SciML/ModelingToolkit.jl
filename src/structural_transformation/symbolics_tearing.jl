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

function tearing_substitute_expr(sys::AbstractSystem, expr; simplify = false)
    empty_substitutions(sys) && return expr
    substitutions = get_substitutions(sys)
    @unpack subs = substitutions
    solved = Dict(eq.lhs => eq.rhs for eq in subs)
    return tearing_sub(expr, solved, simplify)
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
    @set! sys.schedule = nothing
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
function substitute_derivatives_algevars!(ts::TearingState, neweqs, var_eq_matching, dummy_sub; iv = nothing, D = nothing)
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

The lowest-order term is Shift(t, k)(x(t)), instead of x(t).
As such we actually want dummy variables for the k-1 lowest order terms instead of the k-1 highest order terms.

Shift(t, -1)(x(t)) -> x\_{t-1}(t)

Since Shift(t, -1)(x) is not a derivative, it is directly substituted in `fullvars`. No equation or variable is added for it. 

For ODESystems D(D(D(x))) in equations is recursively substituted as D(x) ~ x_t, D(x_t) ~ x_tt, etc. The analogue for discrete systems, Shift(t, 1)(Shift(t,1)(Shift(t,1)(Shift(t, -3)(x(t))))) does not actually appear. So `total_sub` in generate_system_equations` is directly initialized with all of the lowered variables `Shift(t, -3)(x) -> x_t-3(t)`, etc. 
=#
"""
Generate new derivative variables for the system.

Effects on the system structure: 
- fullvars: add the new derivative variables x_t
- neweqs: add the identity equations for the new variables, D(x) ~ x_t
- graph: update graph with the new equations and variables, and their connections
- solvable_graph:
- var_eq_matching: match D(x) to the added identity equation
"""
function generate_derivative_variables!(ts::TearingState, neweqs, var_eq_matching; mm = nothing, iv = nothing, D = nothing)
    @unpack fullvars, sys, structure = ts
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    eq_var_matching = invview(var_eq_matching)
    diff_to_var = invview(var_to_diff)
    is_discrete = is_only_discrete(structure)
    lower_varname = is_discrete ? lower_shift_varname : lower_varname_with_unit
    linear_eqs = mm === nothing ? Dict{Int, Int}() :
                 Dict(reverse(en) for en in enumerate(mm.nzrows))

    # For variable x, make dummy derivative x_t if the
    # derivative is in the system
    for v in 1:length(var_to_diff)
        dv = var_to_diff[v]
        # For discrete systems, directly substitute lowest-order shift 
        if is_discrete && diff_to_var[v] == nothing
            operation(fullvars[v]) isa Shift && (fullvars[v] = lower_varname(fullvars[v], iv))
        end
        dv isa Int || continue
        solved = var_eq_matching[dv] isa Int
        solved && continue

        # If there's `D(x) = x_t` already, update mappings and continue without
        # adding new equations/variables
        dd = find_duplicate_dd(dv, solvable_graph, diff_to_var, linear_eqs, mm)
        if !isnothing(dd)
            dummy_eq, v_t = dd
            var_to_diff[v_t] = var_to_diff[dv]
            var_eq_matching[dv] = unassigned
            eq_var_matching[dummy_eq] = dv
            continue
        end

        dx = fullvars[dv]
        order, lv = var_order(dv, diff_to_var)
        x_t = is_discrete ? lower_varname(fullvars[dv], iv) : lower_varname(fullvars[lv], iv, order)

        # Add `x_t` to the graph
        v_t = add_dd_variable!(structure, fullvars, x_t, dv)
        # Add `D(x) - x_t ~ 0` to the graph
        dummy_eq = add_dd_equation!(structure, neweqs, 0 ~ dx - x_t, dv, v_t)

        # Update matching
        push!(var_eq_matching, unassigned)
        var_eq_matching[dv] = unassigned
        eq_var_matching[dummy_eq] = dv 
    end
end

"""
Check if there's `D(x) = x_t` already.
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

# dv = index of D(x), v_t = index of x_t
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
Solve the solvable equations of the system and generate differential (or discrete)
equations in terms of the selected states.
"""
function generate_system_equations!(state::TearingState, neweqs, var_eq_matching; simplify = false, iv = nothing, D = nothing)
    @unpack fullvars, sys, structure = state 
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = structure
    eq_var_matching = invview(var_eq_matching)
    diff_to_var = invview(var_to_diff)

    total_sub = Dict()
    if is_only_discrete(structure)
        for v in fullvars
            op = operation(v)
            op isa Shift && (op.steps < 0) && (total_sub[v] = lower_shift_varname(v, iv))
        end
   end

    # if var is like D(x) or Shift(t, 1)(x)
    isdervar = let diff_to_var = diff_to_var
        var -> diff_to_var[var] !== nothing
    end

    # Extract partition information
    is_solvable = let solvable_graph = solvable_graph
        (eq, iv) -> eq isa Int && iv isa Int && BipartiteEdge(eq, iv) in solvable_graph
    end

    diff_eqs = Equation[]
    diffeq_idxs = Int[]
    diff_vars = Int[]
    alge_eqs = Equation[]
    algeeq_idxs = Int[]
    solved_eqs = Equation[]
    solvedeq_idxs = Int[]
    solved_vars = Int[]

    toporder = topological_sort(DiCMOBiGraph{false}(graph, var_eq_matching))
    eqs = Iterators.reverse(toporder)
    idep = iv

    # Generate equations.
    #   Solvable equations of differential variables D(x) become differential equations
    #   Solvable equations of non-differential variables become observable equations
    #   Non-solvable equations become algebraic equations.
    for ieq in eqs
        iv = eq_var_matching[ieq]
        eq = neweqs[ieq]

        if is_solvable(ieq, iv) && isdervar(iv)
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
            push!(diff_eqs, neweq)
            push!(diffeq_idxs, ieq)
            push!(diff_vars, diff_to_var[iv])
        elseif is_solvable(ieq, iv)
            var = fullvars[iv]
            neweq = make_solved_equation(var, eq, total_sub; simplify)
            !isnothing(neweq) && begin
                push!(solved_eqs, neweq)
                push!(solvedeq_idxs, ieq)
                push!(solved_vars, iv)
            end
        else
            neweq = make_algebraic_equation(eq, total_sub)
            push!(alge_eqs, neweq)
            push!(algeeq_idxs, ieq)
        end
    end

    # Generate new equations and orderings 
    neweqs = [diff_eqs; alge_eqs]
    eq_ordering = [diffeq_idxs; algeeq_idxs]
    diff_vars_set = BitSet(diff_vars)
    if length(diff_vars_set) != length(diff_vars)
        error("Tearing internal error: lowering DAE into semi-implicit ODE failed!")
    end
    solved_vars_set = BitSet(solved_vars)
    var_ordering = [diff_vars;
                   setdiff!(setdiff(1:ndsts(graph), diff_vars_set),
                       solved_vars_set)]

    return neweqs, solved_eqs, eq_ordering, var_ordering, length(solved_vars), length(solved_vars_set)
end

struct UnexpectedDifferentialError
    eq::Equation
end

function Base.showerror(io::IO, err::UnexpectedDifferentialError)
    error("Differential found in a non-differential system. Likely this is a bug in the construction of an initialization system. Please report this issue with a reproducible example. Offending equation: $(err.eq)")
end

function make_differential_equation(var, dx, eq, total_sub)
    dx ~ simplify_shifts(Symbolics.fixpoint_sub(
        Symbolics.symbolic_linear_solve(eq, var),
        total_sub; operator = ModelingToolkit.Shift))
end

function make_algebraic_equation(eq, total_sub)
    rhs = eq.rhs
    if !(eq.lhs isa Number && eq.lhs == 0)
        rhs = eq.rhs - eq.lhs
    end
    0 ~ simplify_shifts(Symbolics.fixpoint_sub(rhs, total_sub))
end

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
Reorder the equations and unknowns to be:
   [diffeqs; ...]
   [diffvars; ...]
such that the mass matrix is:
   [I  0
    0  0].

Update the state to account for the new ordering and equations.
"""
# TODO: BLT sorting
function reorder_vars!(state::TearingState, var_eq_matching, eq_ordering, var_ordering, nelim_eq, nelim_var)
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
        nelim_eq, nelim_var)

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

    # Update system structure 
    @set! state.structure.graph = complete(new_graph)
    @set! state.structure.var_to_diff = new_var_to_diff
    @set! state.structure.eq_to_diff = new_eq_to_diff
    @set! state.fullvars = new_fullvars
    state
end

"""
Set the system equations, unknowns, observables post-tearing.
"""
function update_simplified_system!(state::TearingState, neweqs, solved_eqs, dummy_sub, var_eq_matching, extra_unknowns; 
        cse_hack = true, array_hack = true)
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = state.structure
    diff_to_var = invview(var_to_diff)

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
    obs = [fast_substitute(observed(sys), obs_sub); solved_eqs]

    unknowns = Any[v
                   for (i, v) in enumerate(state.fullvars)
                   if diff_to_var[i] === nothing && ispresent(i)]
    unknowns = [unknowns; extra_unknowns]
    @set! sys.unknowns = unknowns

    obs, subeqs, deps = cse_and_array_hacks(
        sys, obs, solved_eqs, unknowns, neweqs; cse = cse_hack, array = array_hack)

    @set! sys.eqs = neweqs
    @set! sys.observed = obs
    @set! sys.substitutions = Substitutions(subeqs, deps)

    # Only makes sense for time-dependent
    # TODO: generalize to SDE
    if sys isa ODESystem
        @set! sys.schedule = Schedule(var_eq_matching, dummy_sub)
    end
    sys = schedule(sys)
end

# Terminology and Definition:
# A general DAE is in the form of `F(u'(t), u(t), p, t) == 0`. We can
# characterize variables in `u(t)` into two classes: differential variables
# (denoted `v(t)`) and algebraic variables (denoted `z(t)`). Differential
# variables are marked as `SelectedState` and they are differentiated in the
# DAE system, i.e. `v'(t)` are all the variables in `u'(t)` that actually
# appear in the system. Algebraic variables are variables that are not
# differential variables.
    
# Give the order of the variable indexed by dv
function var_order(dv, diff_to_var)
    order = 0
    while (dv′ = diff_to_var[dv]) !== nothing
        order += 1
        dv = dv′
    end
    order, dv
end

function tearing_reassemble(state::TearingState, var_eq_matching,
        full_var_eq_matching = nothing; simplify = false, mm = nothing, cse_hack = true, array_hack = true)
    extra_vars = Int[]
    if full_var_eq_matching !== nothing
        for v in 𝑑vertices(state.structure.graph)
            eq = full_var_eq_matching[v]
            eq isa Int && continue
            push!(extra_vars, v)
        end
    end
    extra_unknowns = state.fullvars[extra_vars]
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

    # Structural simplification 
    substitute_derivatives_algevars!(state, neweqs, var_eq_matching, dummy_sub; iv, D)

    generate_derivative_variables!(state, neweqs, var_eq_matching; mm, iv, D)

    neweqs, solved_eqs, eq_ordering, var_ordering, nelim_eq, nelim_var = 
        generate_system_equations!(state, neweqs, var_eq_matching; simplify, iv, D)

    state = reorder_vars!(state, var_eq_matching, eq_ordering, var_ordering, nelim_eq, nelim_var)

    sys = update_simplified_system!(state, neweqs, solved_eqs, dummy_sub, var_eq_matching, extra_unknowns; cse_hack, array_hack)

    @set! state.sys = sys
    @set! sys.tearing_state = state
    return invalidate_cache!(sys)
end

"""
# HACK 1

Since we don't support array equations, any equation of the sort `x[1:n] ~ f(...)[1:n]`
gets turned into `x[1] ~ f(...)[1], x[2] ~ f(...)[2]`. Repeatedly calling `f` gets
_very_ expensive. this hack performs a limited form of CSE specifically for this case to
avoid the unnecessary cost. This and the below hack are implemented simultaneously

# HACK 2

Add equations for array observed variables. If `p[i] ~ (...)` are equations, add an
equation `p ~ [p[1], p[2], ...]` allow topsort to reorder them only add the new equation
if all `p[i]` are present and the unscalarized form is used in any equation (observed or
not) we first count the number of times the scalarized form of each observed variable
occurs in observed equations (and unknowns if it's split).
"""
function cse_and_array_hacks(sys, obs, subeqs, unknowns, neweqs; cse = true, array = true)
    # HACK 1
    # mapping of rhs to temporary CSE variable
    # `f(...) => tmpvar` in above example
    rhs_to_tempvar = Dict()

    # HACK 2
    # map of array observed variable (unscalarized) to number of its
    # scalarized terms that appear in observed equations
    arr_obs_occurrences = Dict()
    # to check if array variables occur in unscalarized form anywhere
    all_vars = Set()
    for (i, eq) in enumerate(obs)
        lhs = eq.lhs
        rhs = eq.rhs
        vars!(all_vars, rhs)

        # HACK 1
        if cse && is_getindexed_array(rhs)
            rhs_arr = arguments(rhs)[1]
            if !haskey(rhs_to_tempvar, rhs_arr)
                tempvar = gensym(Symbol(lhs))
                N = length(rhs_arr)
                tempvar = unwrap(Symbolics.variable(
                    tempvar; T = Symbolics.symtype(rhs_arr)))
                tempvar = setmetadata(
                    tempvar, Symbolics.ArrayShapeCtx, Symbolics.shape(rhs_arr))
                tempeq = tempvar ~ rhs_arr
                rhs_to_tempvar[rhs_arr] = tempvar
                push!(obs, tempeq)
                push!(subeqs, tempeq)
            end

            # getindex_wrapper is used because `observed2graph` treats `x` and `x[i]` as different,
            # so it doesn't find a dependency between this equation and `tempvar ~ rhs_arr`
            # which fails the topological sort
            neweq = lhs ~ getindex_wrapper(
                rhs_to_tempvar[rhs_arr], Tuple(arguments(rhs)[2:end]))
            obs[i] = neweq
            subeqi = findfirst(isequal(eq), subeqs)
            if subeqi !== nothing
                subeqs[subeqi] = neweq
            end
        end
        # end HACK 1

        array || continue
        iscall(lhs) || continue
        operation(lhs) === getindex || continue
        Symbolics.shape(lhs) != Symbolics.Unknown() || continue
        arg1 = arguments(lhs)[1]
        cnt = get(arr_obs_occurrences, arg1, 0)
        arr_obs_occurrences[arg1] = cnt + 1
        continue
    end

    # Also do CSE for `equations(sys)`
    if cse
        for (i, eq) in enumerate(neweqs)
            (; lhs, rhs) = eq
            is_getindexed_array(rhs) || continue
            rhs_arr = arguments(rhs)[1]
            if !haskey(rhs_to_tempvar, rhs_arr)
                tempvar = gensym(Symbol(lhs))
                N = length(rhs_arr)
                tempvar = unwrap(Symbolics.variable(
                    tempvar; T = Symbolics.symtype(rhs_arr)))
                tempvar = setmetadata(
                    tempvar, Symbolics.ArrayShapeCtx, Symbolics.shape(rhs_arr))
                vars!(all_vars, rhs_arr)
                tempeq = tempvar ~ rhs_arr
                rhs_to_tempvar[rhs_arr] = tempvar
                push!(obs, tempeq)
                push!(subeqs, tempeq)
            end
            # don't need getindex_wrapper, but do it anyway to know that this
            # hack took place
            neweq = lhs ~ getindex_wrapper(
                rhs_to_tempvar[rhs_arr], Tuple(arguments(rhs)[2:end]))
            neweqs[i] = neweq
        end
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
        arr_obs_occurrences[arg1] = cnt + 1
    end
    for eq in neweqs
        vars!(all_vars, eq.rhs)
    end

    # also count unscalarized variables used in callbacks
    for ev in Iterators.flatten((continuous_events(sys), discrete_events(sys)))
        vars!(all_vars, ev)
    end
    obs_arr_eqs = Equation[]
    for (arrvar, cnt) in arr_obs_occurrences
        cnt == length(arrvar) || continue
        arrvar in all_vars || continue
        # firstindex returns 1 for multidimensional array symbolics
        firstind = first(eachindex(arrvar))
        scal = [arrvar[i] for i in eachindex(arrvar)]
        # respect non-1-indexed arrays
        # TODO: get rid of this hack together with the above hack, then remove OffsetArrays dependency
        # `change_origin` is required because `Origin(firstind)(scal)` makes codegen
        # try to `create_array(OffsetArray{...}, ...)` which errors.
        # `term(Origin(firstind), scal)` doesn't retain the `symtype` and `size`
        # of `scal`.
        rhs = scal
        rhs = change_origin(firstind, rhs)
        push!(obs_arr_eqs, arrvar ~ rhs)
    end
    append!(obs, obs_arr_eqs)
    append!(subeqs, obs_arr_eqs)

    # need to re-sort subeqs
    subeqs = ModelingToolkit.topsort_equations(subeqs, [eq.lhs for eq in subeqs])

    deps = Vector{Int}[i == 1 ? Int[] : collect(1:(i - 1))
                       for i in 1:length(subeqs)]

    return obs, subeqs, deps
end

function is_getindexed_array(rhs)
    (!ModelingToolkit.isvariable(rhs) || ModelingToolkit.iscalledparameter(rhs)) &&
        iscall(rhs) && operation(rhs) === getindex &&
        Symbolics.shape(rhs) != Symbolics.Unknown()
end

# PART OF HACK 1
getindex_wrapper(x, i) = x[i...]

@register_symbolic getindex_wrapper(x::AbstractArray, i::Tuple{Vararg{Int}})

# PART OF HACK 2
function change_origin(origin, arr)
    if all(isone, Tuple(origin))
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

"""
    tearing(sys; simplify=false)

Tear the nonlinear equations in system. When `simplify=true`, we simplify the
new residual equations after tearing. End users are encouraged to call [`structural_simplify`](@ref)
instead, which calls this function internally.
"""
function tearing(sys::AbstractSystem, state = TearingState(sys); mm = nothing,
        simplify = false, cse_hack = true, array_hack = true, kwargs...)
    var_eq_matching, full_var_eq_matching = tearing(state)
    invalidate_cache!(tearing_reassemble(
        state, var_eq_matching, full_var_eq_matching; mm, simplify, cse_hack, array_hack))
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
        mm = nothing, cse_hack = true, array_hack = true, kwargs...)
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
    tearing_reassemble(state, var_eq_matching; simplify, mm, cse_hack, array_hack)
end
