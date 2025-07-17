###
### Bipartite graph utilities
###

"""
    maximal_matching(s::SystemStructure, eqfilter=eq->true, varfilter=v->true) -> Matching

Find equation-variable maximal bipartite matching. `s.graph` is a bipartite graph.
"""
function BipartiteGraphs.maximal_matching(s::SystemStructure, eqfilter = eq -> true,
        varfilter = v -> true)
    maximal_matching(s.graph, eqfilter, varfilter)
end

n_concrete_eqs(state::TransformationState) = n_concrete_eqs(state.structure)
n_concrete_eqs(structure::SystemStructure) = n_concrete_eqs(structure.graph)
function n_concrete_eqs(graph::BipartiteGraph)
    neqs = count(e -> !isempty(𝑠neighbors(graph, e)), 𝑠vertices(graph))
end

function error_reporting(state, bad_idxs, n_highest_vars, iseqs, orig_inputs)
    io = IOBuffer()
    neqs = n_concrete_eqs(state)
    if iseqs
        error_title = "More equations than variables, here are the potential extra equation(s):\n"
        out_arr = has_equations(state) ? equations(state)[bad_idxs] : bad_idxs
    else
        error_title = "More variables than equations, here are the potential extra variable(s):\n"
        out_arr = get_fullvars(state)[bad_idxs]
        unset_inputs = intersect(out_arr, orig_inputs)
        n_missing_eqs = n_highest_vars - neqs
        n_unset_inputs = length(unset_inputs)
        if n_unset_inputs > 0
            println(io, "In particular, the unset input(s) are:")
            Base.print_array(io, unset_inputs)
            println(io)
            println(io, "The rest of potentially unset variable(s) are:")
        end
    end

    Base.print_array(io, out_arr)
    msg = String(take!(io))
    if iseqs
        throw(ExtraEquationsSystemException("The system is unbalanced. There are " *
                                            "$n_highest_vars highest order derivative variables "
                                            * "and $neqs equations.\n"
                                            * error_title
                                            * msg))
    else
        throw(ExtraVariablesSystemException("The system is unbalanced. There are " *
                                            "$n_highest_vars highest order derivative variables "
                                            * "and $neqs equations.\n"
                                            * error_title
                                            * msg))
    end
end

###
### Structural check
###

"""
    $(TYPEDSIGNATURES)

Check if the `state` represents a singular system, and return the unmatched variables.
"""
function singular_check(state::TransformationState)
    @unpack graph, var_to_diff = state.structure
    fullvars = get_fullvars(state)
    # This is defined to check if Pantelides algorithm terminates. For more
    # details, check the equation (15) of the original paper.
    extended_graph = (@set graph.fadjlist = Vector{Int}[graph.fadjlist;
                                                        map(collect, edges(var_to_diff))])
    extended_var_eq_matching = maximal_matching(extended_graph)

    nvars = ndsts(graph)
    unassigned_var = []
    for (vj, eq) in enumerate(extended_var_eq_matching)
        vj > nvars && break
        if eq === unassigned && !isempty(𝑑neighbors(graph, vj))
            push!(unassigned_var, fullvars[vj])
        end
    end
    return unassigned_var
end

"""
    $(TYPEDSIGNATURES)

Check the consistency of `state`, given the inputs `orig_inputs`. If `nothrow == false`,
throws an error if the system is under-/over-determined or singular. In this case, if the
function returns it will return `true`. If `nothrow == true`, it will return `false`
instead of throwing an error. The singular case will print a warning.
"""
function check_consistency(state::TransformationState, orig_inputs; nothrow = false)
    fullvars = get_fullvars(state)
    neqs = n_concrete_eqs(state)
    @unpack graph, var_to_diff = state.structure
    highest_vars = computed_highest_diff_variables(complete!(state.structure))
    n_highest_vars = 0
    for (v, h) in enumerate(highest_vars)
        h || continue
        isempty(𝑑neighbors(graph, v)) && continue
        n_highest_vars += 1
    end
    is_balanced = n_highest_vars == neqs

    if neqs > 0 && !is_balanced
        nothrow && return false
        varwhitelist = var_to_diff .== nothing
        var_eq_matching = maximal_matching(graph, eq -> true, v -> varwhitelist[v]) # not assigned
        # Just use `error_reporting` to do conditional
        iseqs = n_highest_vars < neqs
        if iseqs
            eq_var_matching = invview(complete(var_eq_matching, nsrcs(graph))) # extra equations
            bad_idxs = findall(isequal(unassigned), @view eq_var_matching[1:nsrcs(graph)])
        else
            bad_idxs = findall(isequal(unassigned), var_eq_matching)
        end
        error_reporting(state, bad_idxs, n_highest_vars, iseqs, orig_inputs)
    end

    unassigned_var = singular_check(state)

    if !isempty(unassigned_var) || !is_balanced
        if nothrow
            return false
        end
        io = IOBuffer()
        Base.print_array(io, unassigned_var)
        unassigned_var_str = String(take!(io))
        errmsg = "The system is structurally singular! " *
                 "Here are the problematic variables: \n" *
                 unassigned_var_str
        throw(InvalidSystemException(errmsg))
    end

    return true
end

###
### BLT ordering
###

"""
    find_var_sccs(g::BipartiteGraph, assign=nothing)

Find strongly connected components of the variables defined by `g`. `assign`
gives the undirected bipartite graph a direction. When `assign === nothing`, we
assume that the ``i``-th variable is assigned to the ``i``-th equation.
"""
function find_var_sccs(g::BipartiteGraph, assign = nothing)
    cmog = DiCMOBiGraph{true}(g,
        Matching(assign === nothing ? Base.OneTo(nsrcs(g)) : assign))
    sccs = Graphs.strongly_connected_components(cmog)
    cgraph = MatchedCondensationGraph(cmog, sccs)
    toporder = topological_sort(cgraph)
    permute!(sccs, toporder)
    foreach(sort!, sccs)
    return sccs
end

function sorted_incidence_matrix(ts::TransformationState, val = true; only_algeqs = false,
        only_algvars = false)
    var_eq_matching, var_scc = algebraic_variables_scc(ts)
    s = ts.structure
    graph = ts.structure.graph
    varsmap = zeros(Int, ndsts(graph))
    eqsmap = zeros(Int, nsrcs(graph))
    varidx = 0
    eqidx = 0
    for vs in var_scc, v in vs

        eq = var_eq_matching[v]
        if eq !== unassigned
            eqsmap[eq] = (eqidx += 1)
            varsmap[v] = (varidx += 1)
        end
    end
    for i in diffvars_range(s)
        varsmap[i] = (varidx += 1)
    end
    for i in dervars_range(s)
        varsmap[i] = (varidx += 1)
    end
    for i in 1:nsrcs(graph)
        if eqsmap[i] == 0
            eqsmap[i] = (eqidx += 1)
        end
    end

    I = Int[]
    J = Int[]
    algeqs_set = algeqs(s)
    for eq in 𝑠vertices(graph)
        only_algeqs && (eq in algeqs_set || continue)
        for var in 𝑠neighbors(graph, eq)
            only_algvars && (isalgvar(s, var) || continue)
            i = eqsmap[eq]
            j = varsmap[var]
            (iszero(i) || iszero(j)) && continue
            push!(I, i)
            push!(J, j)
        end
    end
    sparse(I, J, val, nsrcs(graph), ndsts(graph))
end

"""
    $(TYPEDSIGNATURES)

Obtain the incidence matrix of the system sorted by the SCCs. Requires that the system is
simplified and has a `schedule`.
"""
function sorted_incidence_matrix(sys::AbstractSystem)
    if !iscomplete(sys) || get_tearing_state(sys) === nothing ||
       get_schedule(sys) === nothing
        error("A simplified `System` is required. Call `mtkcompile` on the system before creating an `SCCNonlinearProblem`.")
    end
    sched = get_schedule(sys)
    var_sccs = sched.var_sccs

    ts = get_tearing_state(sys)
    imat = Graphs.incidence_matrix(ts.structure.graph)
    buffer = similar(imat)
    permute!(buffer, imat, 1:size(imat, 2), reduce(vcat, var_sccs))
    buffer
end

###
### Structural and symbolic utilities
###

function find_eq_solvables!(state::TearingState, ieq, to_rm = Int[], coeffs = nothing;
        may_be_zero = false,
        allow_symbolic = false, allow_parameter = true,
        conservative = false,
        kwargs...)
    fullvars = state.fullvars
    @unpack graph, solvable_graph = state.structure
    eq = equations(state)[ieq]
    term = value(eq.rhs - eq.lhs)
    all_int_vars = true
    coeffs === nothing || empty!(coeffs)
    empty!(to_rm)
    for j in 𝑠neighbors(graph, ieq)
        var = fullvars[j]
        isirreducible(var) && (all_int_vars = false; continue)
        a, b, islinear = linear_expansion(term, var)
        a, b = unwrap(a), unwrap(b)
        islinear || (all_int_vars = false; continue)
        if a isa Symbolic
            all_int_vars = false
            if !allow_symbolic
                if allow_parameter
                    # if any of the variables in `a` are present in fullvars (taking into account arrays)
                    if any(
                        v -> any(isequal(v), fullvars) ||
                             symbolic_type(v) == ArraySymbolic() &&
                             Symbolics.shape(v) != Symbolics.Unknown() &&
                             any(x -> any(isequal(x), fullvars), collect(v)),
                        vars(
                            a; op = Union{Differential, Shift, Pre, Sample, Hold, Initial}))
                        continue
                    end
                else
                    continue
                end
            end
            add_edge!(solvable_graph, ieq, j)
            continue
        end
        if !(a isa Number)
            all_int_vars = false
            continue
        end
        # When the expression is linear with numeric `a`, then we can safely
        # only consider `b` for the following iterations.
        term = b
        if isone(abs(a))
            coeffs === nothing || push!(coeffs, convert(Int, a))
        else
            all_int_vars = false
            conservative && continue
        end
        if a != 0
            add_edge!(solvable_graph, ieq, j)
        else
            if may_be_zero
                push!(to_rm, j)
            else
                @warn "Internal error: Variable $var was marked as being in $eq, but was actually zero"
            end
        end
    end
    for j in to_rm
        rem_edge!(graph, ieq, j)
    end
    all_int_vars, term
end

function find_solvables!(state::TearingState; kwargs...)
    @assert state.structure.solvable_graph === nothing
    eqs = equations(state)
    graph = state.structure.graph
    state.structure.solvable_graph = BipartiteGraph(nsrcs(graph), ndsts(graph))
    to_rm = Int[]
    for ieq in 1:length(eqs)
        find_eq_solvables!(state, ieq, to_rm; kwargs...)
    end
    return nothing
end

function linear_subsys_adjmat!(state::TransformationState; kwargs...)
    graph = state.structure.graph
    if state.structure.solvable_graph === nothing
        state.structure.solvable_graph = BipartiteGraph(nsrcs(graph), ndsts(graph))
    end
    linear_equations = Int[]
    eqs = equations(state.sys)
    eadj = Vector{Int}[]
    cadj = Vector{Int}[]
    coeffs = Int[]
    to_rm = Int[]
    for i in eachindex(eqs)
        all_int_vars, rhs = find_eq_solvables!(state, i, to_rm, coeffs; kwargs...)

        # Check if all unknowns in the equation is both linear and homogeneous,
        # i.e. it is in the form of
        #
        #       ``∑ c_i * v_i = 0``,
        #
        # where ``c_i`` ∈ ℤ and ``v_i`` denotes unknowns.
        if all_int_vars && Symbolics._iszero(rhs)
            push!(linear_equations, i)
            push!(eadj, copy(𝑠neighbors(graph, i)))
            push!(cadj, copy(coeffs))
        end
    end

    mm = SparseMatrixCLIL(nsrcs(graph),
        ndsts(graph),
        linear_equations, eadj, cadj)
    return mm
end

highest_order_variable_mask(ts) =
    let v2d = ts.structure.var_to_diff
        v -> isempty(outneighbors(v2d, v))
    end

lowest_order_variable_mask(ts) =
    let v2d = ts.structure.var_to_diff
        v -> isempty(inneighbors(v2d, v))
    end

function but_ordered_incidence(ts::TearingState, varmask = highest_order_variable_mask(ts))
    graph = complete(ts.structure.graph)
    var_eq_matching = complete(maximal_matching(graph, _ -> true, varmask))
    scc = find_var_sccs(graph, var_eq_matching)
    vordering = Vector{Int}(undef, 0)
    bb = Int[1]
    sizehint!(vordering, ndsts(graph))
    sizehint!(bb, ndsts(graph))
    l = 1
    for c in scc
        isemptyc = true
        for v in c
            if varmask(v)
                push!(vordering, v)
                l += 1
                isemptyc = false
            end
        end
        isemptyc || push!(bb, l)
    end
    mm = incidence_matrix(graph)
    reverse!(vordering)
    mm[[var_eq_matching[v] for v in vordering if var_eq_matching[v] isa Int], vordering], bb
end

"""
    $(TYPEDSIGNATURES)

Given a system `sys` and torn variable-equation matching `torn_matching`, return the sparse
incidence matrix of the system with SCCs grouped together, and each SCC sorted to contain
the analytically solved equations/variables before the unsolved ones.
"""
function reordered_matrix(sys::System, torn_matching)
    s = TearingState(sys)
    complete!(s.structure)
    @unpack graph = s.structure
    eqs = equations(sys)
    nvars = ndsts(graph)
    max_matching = complete(maximal_matching(graph))
    torn_matching = complete(torn_matching)
    sccs = find_var_sccs(graph, max_matching)
    I, J = Int[], Int[]
    ii = 0
    M = Int[]
    solved = BitSet(findall(torn_matching .!== unassigned))
    for vars in sccs
        append!(M, filter(in(solved), vars))
        append!(M, filter(!in(solved), vars))
    end
    M = invperm(vcat(M, setdiff(1:nvars, M)))
    for vars in sccs
        e_solved = [torn_matching[v] for v in vars if torn_matching[v] !== unassigned]
        for es in e_solved
            isdiffeq(eqs[es]) && continue
            ii += 1
            js = [M[x] for x in 𝑠neighbors(graph, es) if isalgvar(s.structure, x)]
            append!(I, fill(ii, length(js)))
            append!(J, js)
        end

        e_residual = setdiff(
            [max_matching[v]
             for v in vars if max_matching[v] !== unassigned], e_solved)
        for er in e_residual
            isdiffeq(eqs[er]) && continue
            ii += 1
            js = [M[x] for x in 𝑠neighbors(graph, er) if isalgvar(s.structure, x)]
            append!(I, fill(ii, length(js)))
            append!(J, js)
        end
    end
    # only plot algebraic variables and equations
    sparse(I, J, true)
end

"""
    uneven_invmap(n::Int, list)

returns an uneven inv map with length `n`.
"""
function uneven_invmap(n::Int, list)
    rename = zeros(Int, n)
    for (i, v) in enumerate(list)
        rename[v] = i
    end
    return rename
end

function torn_system_jacobian_sparsity(sys)
    state = get_tearing_state(sys)
    state isa TearingState || return nothing
    @unpack structure = state
    @unpack graph, var_to_diff = structure

    neqs = nsrcs(graph)
    nsts = ndsts(graph)
    states_idxs = findall(!Base.Fix1(isdervar, structure), 1:nsts)
    var2idx = uneven_invmap(nsts, states_idxs)
    I = Int[]
    J = Int[]
    for ieq in 1:neqs
        for ivar in 𝑠neighbors(graph, ieq)
            nivar = get(var2idx, ivar, 0)
            nivar == 0 && continue
            push!(I, ieq)
            push!(J, nivar)
        end
    end
    return sparse(I, J, true, neqs, neqs)
end

###
### Misc
###

"""
Handle renaming variable names for discrete structural simplification. Three cases: 
- positive shift: do nothing
- zero shift: x(t) => Shift(t, 0)(x(t))
- negative shift: rename the variable
"""
function lower_shift_varname(var, iv)
    op = operation(var)
    if op isa Shift
        return shift2term(var)
    else
        return Shift(iv, 0)(var, true)
    end
end

function descend_lower_shift_varname_with_unit(var, iv)
    symbolic_type(var) == NotSymbolic() && return var
    ModelingToolkit._with_unit(descend_lower_shift_varname, var, iv, iv)
end
function descend_lower_shift_varname(var, iv)
    iscall(var) || return var
    op = operation(var)
    if op isa Shift
        return shift2term(var)
    else
        args = arguments(var)
        args = map(Base.Fix2(descend_lower_shift_varname, iv), args)
        return maketerm(typeof(var), op, args, Symbolics.metadata(var))
    end
end

"""
Rename a Shift variable with negative shift, Shift(t, k)(x(t)) to xₜ₋ₖ(t).
"""
function shift2term(var)
    iscall(var) || return var
    op = operation(var)
    op isa Shift || return var
    iv = op.t
    arg = only(arguments(var))
    if operation(arg) === getindex
        idxs = arguments(arg)[2:end]
        newvar = shift2term(op(first(arguments(arg))))[idxs...]
        unshifted = ModelingToolkit.getunshifted(newvar)[idxs...]
        newvar = setmetadata(newvar, ModelingToolkit.VariableUnshifted, unshifted)
        return newvar
    end
    is_lowered = !isnothing(ModelingToolkit.getunshifted(arg))

    backshift = is_lowered ? op.steps + ModelingToolkit.getshift(arg) : op.steps

    # Char(0x208b) = ₋ (subscripted minus)
    # Char(0x208a) = ₊ (subscripted plus)
    pm = backshift > 0 ? Char(0x208a) : Char(0x208b)
    # subscripted number, e.g. ₁
    num = join(Char(0x2080 + d) for d in reverse!(digits(abs(backshift))))
    # Char(0x209c) = ₜ
    # ds = ₜ₋₁
    ds = join([Char(0x209c), pm, num])

    O = is_lowered ? ModelingToolkit.getunshifted(arg) : arg
    oldop = operation(O)
    newname = backshift != 0 ? Symbol(string(nameof(oldop)), ds) :
              Symbol(string(nameof(oldop)))

    newvar = maketerm(typeof(O), Symbolics.rename(oldop, newname),
        Symbolics.children(O), Symbolics.metadata(O))
    newvar = setmetadata(newvar, Symbolics.VariableSource, (:variables, newname))
    newvar = setmetadata(newvar, ModelingToolkit.VariableUnshifted, O)
    newvar = setmetadata(newvar, ModelingToolkit.VariableShift, backshift)
    return newvar
end

function isdoubleshift(var)
    return ModelingToolkit.isoperator(var, ModelingToolkit.Shift) &&
           ModelingToolkit.isoperator(arguments(var)[1], ModelingToolkit.Shift)
end

"""
Simplify multiple shifts: Shift(t, k1)(Shift(t, k2)(x)) becomes Shift(t, k1+k2)(x).
"""
function simplify_shifts(var)
    ModelingToolkit.hasshift(var) || return var
    var isa Equation && return simplify_shifts(var.lhs) ~ simplify_shifts(var.rhs)
    (op = operation(var)) isa Shift && op.steps == 0 && return first(arguments(var))
    if isdoubleshift(var)
        op1 = operation(var)
        vv1 = arguments(var)[1]
        op2 = operation(vv1)
        vv2 = arguments(vv1)[1]
        s1 = op1.steps
        s2 = op2.steps
        t1 = op1.t
        t2 = op2.t
        return simplify_shifts(ModelingToolkit.Shift(t1 === nothing ? t2 : t1, s1 + s2)(vv2))
    else
        return maketerm(typeof(var), operation(var), simplify_shifts.(arguments(var)),
            unwrap(var).metadata)
    end
end

"""
Distribute a shift applied to a whole expression or equation. 
Shift(t, 1)(x + y) will become Shift(t, 1)(x) + Shift(t, 1)(y).
Only shifts variables whose independent variable is the same t that appears in the Shift (i.e. constants, time-independent parameters, etc. do not get shifted).
"""
function distribute_shift(var)
    var = unwrap(var)
    var isa Equation && return distribute_shift(var.lhs) ~ distribute_shift(var.rhs)

    ModelingToolkit.hasshift(var) || return var
    shift = operation(var)
    shift isa Shift || return var

    shift = operation(var)
    expr = only(arguments(var))
    if expr isa Equation
        return distribute_shift(shift(expr.lhs)) ~ distribute_shift(shift(expr.rhs))
    end
    shiftexpr = _distribute_shift(expr, shift)
    return simplify_shifts(shiftexpr)
end

function _distribute_shift(expr, shift)
    if iscall(expr)
        op = operation(expr)
        (op isa Union{Pre, Initial, Sample, Hold}) && return expr
        args = arguments(expr)

        if ModelingToolkit.isvariable(expr) && operation(expr) !== getindex &&
           !ModelingToolkit.iscalledparameter(expr)
            (length(args) == 1 && isequal(shift.t, only(args))) ? (return shift(expr)) :
            (return expr)
        elseif op isa Shift
            return shift(expr)
        else
            return maketerm(
                typeof(expr), operation(expr), Base.Fix2(_distribute_shift, shift).(args),
                unwrap(expr).metadata)
        end
    else
        return expr
    end
end
