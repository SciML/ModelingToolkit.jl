###
### Reassemble: structural information -> system
###

function pantelides_reassemble(state::TearingState, var_eq_matching)
    fullvars = state.fullvars
    @unpack var_to_diff, eq_to_diff = state.structure
    sys = state.sys
    # Step 1: write derivative equations
    in_eqs = equations(sys)
    out_eqs = Vector{Any}(undef, nv(eq_to_diff))
    fill!(out_eqs, nothing)
    out_eqs[1:length(in_eqs)] .= in_eqs

    out_vars = Vector{Any}(undef, nv(var_to_diff))
    fill!(out_vars, nothing)
    out_vars[1:length(fullvars)] .= fullvars

    iv = get_iv(sys)
    D = Differential(iv)

    for (varidx, diff) in edges(var_to_diff)
        # fullvars[diff] = D(fullvars[var])
        vi = out_vars[varidx]
        @assert vi!==nothing "Something went wrong on reconstructing unknowns from variable association list"
        # `fullvars[i]` needs to be not a `D(...)`, because we want the DAE to be
        # first-order.
        if isdifferential(vi)
            vi = out_vars[varidx] = diff2term_with_unit(vi, iv)
        end
        out_vars[diff] = D(vi)
    end

    d_dict = Dict(zip(fullvars, 1:length(fullvars)))
    lhss = Set{Any}([x.lhs for x in in_eqs if isdiffeq(x)])
    for (eqidx, diff) in edges(eq_to_diff)
        # LHS variable is looked up from var_to_diff
        # the var_to_diff[i]-th variable is the differentiated version of var at i
        eq = out_eqs[eqidx]
        lhs = if !(eq.lhs isa Symbolic)
            0
        elseif isdiffeq(eq)
            # look up the variable that represents D(lhs)
            lhsarg1 = arguments(eq.lhs)[1]
            @assert !(lhsarg1 isa Differential) "The equation $eq is not first order"
            i = get(d_dict, lhsarg1, nothing)
            if i === nothing
                D(eq.lhs)
            else
                # remove clashing equations
                lhs = Num(nothing)
            end
        else
            D(eq.lhs)
        end
        rhs = ModelingToolkit.expand_derivatives(D(eq.rhs))
        substitution_dict = Dict(x.lhs => x.rhs
        for x in out_eqs if x !== nothing && x.lhs isa Symbolic)
        sub_rhs = substitute(rhs, substitution_dict)
        out_eqs[diff] = lhs ~ sub_rhs
    end

    final_vars = unique(filter(x -> !(operation(x) isa Differential), fullvars))
    final_eqs = map(identity,
        filter(x -> value(x.lhs) !== nothing,
            out_eqs[sort(filter(x -> x !== unassigned, var_eq_matching))]))

    @set! sys.eqs = final_eqs
    @set! sys.unknowns = final_vars
    return sys
end

"""
    computed_highest_diff_variables(structure)

Computes which variables are the "highest-differentiated" for purposes of
pantelides. Ordinarily this is relatively straightforward. However, in our
case, there is one complicating condition:

    We allow variables in the structure graph that don't appear in the
    system at all. What we are interested in is the highest-differentiated
    variable that actually appears in the system.

This function takes care of these complications are returns a boolean array
for every variable, indicating whether it is considered "highest-differentiated".
"""
function computed_highest_diff_variables(structure)
    @unpack graph, var_to_diff = structure

    nvars = length(var_to_diff)
    varwhitelist = falses(nvars)
    for var in 1:nvars
        if var_to_diff[var] === nothing && !varwhitelist[var]
            # This variable is structurally highest-differentiated, but may not actually appear in the
            # system (complication 1 above). Ascend the differentiation graph to find the highest
            # differentiated variable that does appear in the system or the alias graph).
            while isempty(𝑑neighbors(graph, var))
                var′ = invview(var_to_diff)[var]
                var′ === nothing && break
                var = var′
            end
            varwhitelist[var] = true
        end
    end

    # Remove any variables from the varwhitelist for whom a higher-differentiated
    # var is already on the whitelist.
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

"""
    pantelides!(state::TransformationState; kwargs...)

Perform Pantelides algorithm.
"""
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
        # In practice, `maxiters=8000` should never be reached, otherwise, the
        # index would be on the order of thousands.
        for iii in 1:maxiters
            # run matching on (dx, y) variables
            #
            # the derivatives and algebraic variables are zeros in the variable
            # association list
            resize!(vcolor, nvars)
            fill!(vcolor, false)
            resize!(ecolor, neqs)
            fill!(ecolor, false)
            pathfound = construct_augmenting_path!(var_eq_matching, graph, eq′,
                v -> varwhitelist[v], vcolor, ecolor)
            pathfound && break # terminating condition
            if is_only_discrete(state.structure)
                error("The discrete system has high structural index. This is not supported.")
            end
            for var in eachindex(vcolor)
                vcolor[var] || continue
                if var_to_diff[var] === nothing
                    # introduce a new variable
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
                # introduce a new equation
                neqs += 1
                eq_derivative!(state, eq; kwargs...)
            end

            for var in eachindex(vcolor)
                vcolor[var] || continue
                # the newly introduced `var`s and `eq`s have the inherits
                # assignment
                var_eq_matching[var_to_diff[var]] = eq_to_diff[var_eq_matching[var]]
            end
            eq′ = eq_to_diff[eq′]
        end # for _ in 1:maxiters
        pathfound ||
            error("maxiters=$maxiters reached! File a bug report if your system has a reasonable index (<100), and you are using the default `maxiters`. Try to increase the maxiters by `pantelides(sys::ODESystem; maxiters=1_000_000)` if your system has an incredibly high index and it is truly extremely large.")
    end # for k in 1:neqs′

    finalize && for var in 1:ndsts(graph)
        varwhitelist[var] && continue
        var_eq_matching[var] = unassigned
    end
    return var_eq_matching
end

"""
    dae_index_lowering(sys::ODESystem; kwargs...) -> ODESystem

Perform the Pantelides algorithm to transform a higher index DAE to an index 1
DAE. `kwargs` are forwarded to [`pantelides!`](@ref). End users are encouraged to call [`structural_simplify`](@ref)
instead, which calls this function internally.
"""
function dae_index_lowering(sys::ODESystem; kwargs...)
    state = TearingState(sys)
    var_eq_matching = pantelides!(state; finalize = false, kwargs...)
    return invalidate_cache!(pantelides_reassemble(state, var_eq_matching))
end
