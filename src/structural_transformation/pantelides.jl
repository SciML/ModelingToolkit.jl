
const NOTHING_EQ = nothing ~ nothing
function pantelides_reassemble(state::TearingState, var_eq_matching)
    fullvars = state.fullvars
    @unpack var_to_diff, eq_to_diff = state.structure
    sys = state.sys
    in_eqs = equations(sys)
    out_eqs = Vector{Equation}(undef, nv(eq_to_diff))
    fill!(out_eqs, NOTHING_EQ)
    out_eqs[1:length(in_eqs)] .= in_eqs
    out_vars = Vector{SymbolicT}(undef, nv(var_to_diff))
    fill!(out_vars, ModelingToolkit.COMMON_NOTHING)
    out_vars[1:length(fullvars)] .= fullvars
    iv = get_iv(sys)
    D = Differential(iv)
    for (varidx, diff) in edges(var_to_diff)
        vi = out_vars[varidx]
        @assert vi!==ModelingToolkit.COMMON_NOTHING "Something went wrong on reconstructing unknowns from variable association list"
        if isdifferential(vi)
            vi = out_vars[varidx] = diff2term_with_unit(vi, iv)
        end
        out_vars[diff] = D(vi)
    end
    d_dict = Dict{SymbolicT, Int}(zip(fullvars, 1:length(fullvars)))
    for (eqidx, diff) in edges(eq_to_diff)
        eq = out_eqs[eqidx]
        lhs = if SU.isconst(eq.lhs)
            Symbolics.COMMON_ZERO
        elseif isdiffeq(eq)
            lhsarg1 = arguments(eq.lhs)[1]
            @assert !(lhsarg1 isa Differential) "The equation $eq is not first order"
            i = get(d_dict, lhsarg1, nothing)
            if i === nothing
                D(eq.lhs)
            else
                lhs = ModelingToolkit.COMMON_NOTHING
            end
        else
            D(eq.lhs)
        end
        rhs = ModelingToolkit.expand_derivatives(D(eq.rhs))
        rhs = substitute(rhs, state.param_derivative_map)
        substitution_dict = Dict(x.lhs => x.rhs
        for x in out_eqs if x !== NOTHING_EQ && !SU.isconst(eq.lhs))
        sub_rhs = substitute(rhs, substitution_dict)
        out_eqs[diff] = lhs ~ sub_rhs
    end
    final_vars = unique(filter(x -> !(operation(x) isa Differential), fullvars))
    final_eqs = map(identity,
        filter(x -> x.lhs !== ModelingToolkit.COMMON_NOTHING,
            out_eqs[sort(filter(x -> x !== unassigned, var_eq_matching))]))
    @set! sys.eqs = final_eqs
    @set! sys.unknowns = final_vars
    return sys
end
""""""
function computed_highest_diff_variables(structure)
    @unpack graph, var_to_diff = structure
    nvars = length(var_to_diff)
    varwhitelist = falses(nvars)
    for var in 1:nvars
        if var_to_diff[var] === nothing && !varwhitelist[var]
            while isempty(𝑑neighbors(graph, var))
                var′ = invview(var_to_diff)[var]
                var′ === nothing && break
                var = var′
            end
            varwhitelist[var] = true
        end
    end
    for var in 1:nvars
        varwhitelist[var] || continue
        var′ = var
        while (var′ = var_to_diff[var′]) !== nothing
            if varwhitelist[var′]
                varwhitelist[var] = false
                break
            end
        end
    end
    return varwhitelist
end
""""""
function pantelides!(
        state::TransformationState; finalize = true, maxiters = 8000, kwargs...)
    @unpack graph, solvable_graph, var_to_diff, eq_to_diff = state.structure
    neqs = nsrcs(graph)
    nvars = nv(var_to_diff)
    vcolor = falses(nvars)
    ecolor = falses(neqs)
    var_eq_matching = Matching(nvars)
    neqs′ = neqs
    nnonemptyeqs = count(
        eq -> !isempty(𝑠neighbors(graph, eq)) && eq_to_diff[eq] === nothing,
        1:neqs′)
    varwhitelist = computed_highest_diff_variables(state.structure)
    if nnonemptyeqs > count(varwhitelist)
        throw(InvalidSystemException("System is structurally singular"))
    end
    for k in 1:neqs′
        eq′ = k
        eq_to_diff[eq′] === nothing || continue
        isempty(𝑠neighbors(graph, eq′)) && continue
        pathfound = false
        for iii in 1:maxiters
            resize!(vcolor, nvars)
            fill!(vcolor, false)
            resize!(ecolor, neqs)
            fill!(ecolor, false)
            pathfound = construct_augmenting_path!(var_eq_matching, graph, eq′,
                v -> varwhitelist[v], vcolor, ecolor)
            pathfound && break
            if is_only_discrete(state.structure)
                error("The discrete system has high structural index. This is not supported.")
            end
            for var in eachindex(vcolor)
                vcolor[var] || continue
                if var_to_diff[var] === nothing
                    nvars += 1
                    var_diff = var_derivative!(state, var)
                    push!(var_eq_matching, unassigned)
                    push!(varwhitelist, false)
                    @assert length(var_eq_matching) == var_diff
                end
                varwhitelist[var] = false
                varwhitelist[var_to_diff[var]] = true
            end
            for eq in eachindex(ecolor)
                ecolor[eq] || continue
                neqs += 1
                eq_derivative!(state, eq; kwargs...)
            end
            for var in eachindex(vcolor)
                vcolor[var] || continue
                var_eq_matching[var_to_diff[var]] = eq_to_diff[var_eq_matching[var]]
            end
            eq′ = eq_to_diff[eq′]
        end
        pathfound ||
            error("maxiters=$maxiters reached! File a bug report if your system has a reasonable index (<100), and you are using the default `maxiters`. Try to increase the maxiters by `pantelides(sys::System; maxiters=1_000_000)` if your system has an incredibly high index and it is truly extremely large.")
    end
    finalize && for var in 1:ndsts(graph)
        varwhitelist[var] && continue
        var_eq_matching[var] = unassigned
    end
    return var_eq_matching
end
""""""
function dae_index_lowering(sys::System; kwargs...)
    state = TearingState(sys)
    var_eq_matching = pantelides!(state; finalize = false, kwargs...)
    return invalidate_cache!(pantelides_reassemble(state, var_eq_matching))
end
