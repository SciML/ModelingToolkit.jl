function BipartiteGraphs.maximal_matching(s::SystemStructure, eqfilter = eq -> true,
        varfilter = v -> true)
end
function n_concrete_eqs(graph::BipartiteGraph)
end
function error_reporting(state, bad_idxs, n_highest_vars, iseqs, orig_inputs)
    if iseqs
        if n_unset_inputs > 0
        end
    end
    if iseqs
        throw(ExtraEquationsSystemException("The system is unbalanced. There are " *
                                            "$n_highest_vars highest order derivative variables "
                                            * msg))
        throw(ExtraVariablesSystemException("The system is unbalanced. There are " *
                                            "$n_highest_vars highest order derivative variables "
                                            * msg))
    end
    extended_graph = (@set graph.fadjlist = Vector{Int}[graph.fadjlist;
                                                        map(collect, edges(var_to_diff))])
    for (vj, eq) in enumerate(extended_var_eq_matching)
    end
end
function check_consistency(state::TransformationState, orig_inputs; nothrow = false)
    highest_vars = computed_highest_diff_variables(complete!(state.structure))
    for (v, h) in enumerate(highest_vars)
        if iseqs
        end
    end
end
function find_var_sccs(g::BipartiteGraph, assign = nothing)
    cmog = DiCMOBiGraph{true}(g,
        Matching(assign === nothing ? Base.OneTo(nsrcs(g)) : assign))
    for vs in var_scc, v in vs
        if eq !== unassigned
        end
    end
    for i in diffvars_range(s)
    end
end
function find_eq_solvables!(state::TearingState, ieq, to_rm = Int[], coeffs = nothing;
        kwargs...)
    fullvars = state.fullvars
    @unpack graph, solvable_graph = state.structure
    all_int_vars = true
    __indexed_fullvar_is_var = let fullvars = fullvars
        function indexed_fullvar_is_var(x::SymbolicT)
            for v in fullvars
                Moshi.Match.@match v begin
                end
            end
        end
    end
    __allow_sym_par_cond = let fullvars = fullvars, is_atomic = ModelingToolkit.OperatorIsAtomic{Union{Differential, Shift, Pre, Sample, Hold, Initial}}(), __indexed_fullvar_is_var = __indexed_fullvar_is_var
        function allow_sym_par_cond(v)
                symbolic_type(v) == ArraySymbolic() && (SU.shape(v) isa SU.Unknown ||
                __indexed_fullvar_is_var(v))
         end
    end
    for j in 𝑠neighbors(graph, ieq)
        if !SU.isconst(a)
            if !allow_symbolic
                if allow_parameter
                end
            end
        end
        if SU._isone(abs(unwrap_const(a)))
        end
        if !SU._iszero(a)
        end
    end
    all_int_vars, term
end
function find_solvables!(state::TearingState; kwargs...)
    for ieq in 1:length(eqs)
    end
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
        if all_int_vars && Symbolics._iszero(rhs)
        end
    end
    mm = SparseMatrixCLIL(nsrcs(graph),
        ndsts(graph),
        linear_equations, eadj, cadj)
    end
lowest_order_variable_mask(ts) =
    let v2d = ts.structure.var_to_diff
    for vars in sccs
        for es in e_solved
        end
        e_residual = setdiff(
            [max_matching[v]
             for v in vars if max_matching[v] !== unassigned], e_solved)
    end
end
function torn_system_jacobian_sparsity(sys)
    for ieq in 1:neqs
        for ivar in 𝑠neighbors(graph, ieq)
        end
    end
end
function lower_shift_varname(var, iv)
    if op isa Shift
    end
    Moshi.Match.@match var begin
        BSImpl.Term(f, args) && if f isa Shift end => begin
            Moshi.Match.@match arg begin
                BSImpl.Term(; f, args, type, shape, metadata) && if f === getindex end => begin
                    if metadata === nothing
                    end
                end
            end
            Moshi.Match.@match vv1 begin
                BSImpl.Term(; f = op2, args = a2) && if op2 isa Shift end => begin
                end
            end
        end
    end
end
function simplify_shifts(var::SymbolicT)
end
function distribute_shift(var)
end
function _distribute_shift(expr, shift)
    if iscall(expr)
        if ModelingToolkit.isvariable(expr) && operation(expr) !== getindex &&
            return maketerm(
                unwrap(expr).metadata)
        end
    end
end
