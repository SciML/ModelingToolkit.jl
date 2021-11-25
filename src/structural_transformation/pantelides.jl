###
### Reassemble: structural information -> system
###

function pantelides_reassemble(sys::ODESystem, eq_to_diff, assign)
    s = structure(sys)
    @unpack fullvars, var_to_diff = s
    # Step 1: write derivative equations
    in_eqs = equations(sys)
    out_eqs = Vector{Any}(undef, nv(eq_to_diff))
    fill!(out_eqs, nothing)
    out_eqs[1:length(in_eqs)] .= in_eqs

    out_vars = Vector{Any}(undef, nv(var_to_diff))
    fill!(out_vars, nothing)
    out_vars[1:length(fullvars)] .= fullvars

    D = Differential(get_iv(sys))

    for (varidx, diff) in edges(var_to_diff)
        # fullvars[diff] = D(fullvars[var])
        vi = out_vars[varidx]
        @assert vi !== nothing "Something went wrong on reconstructing states from variable association list"
        # `fullvars[i]` needs to be not a `D(...)`, because we want the DAE to be
        # first-order.
        if isdifferential(vi)
            vi = out_vars[varidx] = diff2term(vi)
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
        substitution_dict = Dict(x.lhs => x.rhs for x in out_eqs if x !== nothing && x.lhs isa Symbolic)
        sub_rhs = substitute(rhs, substitution_dict)
        out_eqs[diff] = lhs ~ sub_rhs
    end

    final_vars = unique(filter(x->!(operation(x) isa Differential), fullvars))
    final_eqs = map(identity, filter(x->value(x.lhs) !== nothing, out_eqs[sort(filter(x->x !== unassigned, assign))]))

    @set! sys.eqs = final_eqs
    @set! sys.states = final_vars
    @set! sys.structure = nothing
    return sys
end

"""
    pantelides!(sys::ODESystem; kwargs...)

Perform Pantelides algorithm.
"""
function pantelides!(sys::ODESystem; maxiters = 8000)
    s = structure(sys)
    # D(j) = assoc[j]
    @unpack graph, var_to_diff = s
    return (sys, pantelides!(graph, var_to_diff)...)
end

function pantelides!(graph, var_to_diff; maxiters = 8000)
    neqs = nsrcs(graph)
    nvars = nv(var_to_diff)
    vcolor = falses(nvars)
    ecolor = falses(neqs)
    var_eq_matching = Matching(nvars)
    eq_to_diff = DiffGraph(neqs)
    neqs′ = neqs
    for k in 1:neqs′
        eq′ = k
        pathfound = false
        # In practice, `maxiters=8000` should never be reached, otherwise, the
        # index would be on the order of thousands.
        for iii in 1:maxiters
            # run matching on (dx, y) variables
            #
            # the derivatives and algebraic variables are zeros in the variable
            # association list
            varwhitelist = var_to_diff .== nothing
            resize!(vcolor, nvars)
            fill!(vcolor, false)
            resize!(ecolor, neqs)
            fill!(ecolor, false)
            pathfound = construct_augmenting_path!(var_eq_matching, graph, eq′, v->varwhitelist[v], vcolor, ecolor)
            pathfound && break # terminating condition
            for var in eachindex(vcolor); vcolor[var] || continue
                # introduce a new variable
                nvars += 1
                add_vertex!(graph, DST)
                # the new variable is the derivative of `var`

                add_edge!(var_to_diff, var, add_vertex!(var_to_diff))
                push!(var_eq_matching, unassigned)
            end

            for eq in eachindex(ecolor); ecolor[eq] || continue
                # introduce a new equation
                neqs += 1
                add_vertex!(graph, SRC)
                # the new equation is created by differentiating `eq`
                eq_diff = add_vertex!(eq_to_diff)
                add_edge!(eq_to_diff, eq, eq_diff)
                for var in 𝑠neighbors(graph, eq)
                    add_edge!(graph, eq_diff, var)
                    add_edge!(graph, eq_diff, var_to_diff[var])
                end
            end

            for var in eachindex(vcolor); vcolor[var] || continue
                # the newly introduced `var`s and `eq`s have the inherits
                # assignment
                var_eq_matching[var_to_diff[var]] = eq_to_diff[var_eq_matching[var]]
            end
            eq′ = eq_to_diff[eq′]
        end # for _ in 1:maxiters
        pathfound || error("maxiters=$maxiters reached! File a bug report if your system has a reasonable index (<100), and you are using the default `maxiters`. Try to increase the maxiters by `pantelides(sys::ODESystem; maxiters=1_000_000)` if your system has an incredibly high index and it is truly extremely large.")
    end # for k in 1:neqs′
    return var_eq_matching, eq_to_diff
end

"""
    dae_index_lowering(sys::ODESystem; kwargs...) -> ODESystem

Perform the Pantelides algorithm to transform a higher index DAE to an index 1
DAE. `kwargs` are forwarded to [`pantelides!`](@ref). End users are encouraged to call [`structural_simplify`](@ref)
instead, which calls this function internally.
"""
function dae_index_lowering(sys::ODESystem; kwargs...)
    s = get_structure(sys)
    (s isa SystemStructure) || (sys = initialize_system_structure(sys))
    sys, var_eq_matching, eq_to_diff = pantelides!(sys; kwargs...)
    return pantelides_reassemble(sys, eq_to_diff, var_eq_matching)
end
