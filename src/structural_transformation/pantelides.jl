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
