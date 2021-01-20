export alias_elimination, flatten

function flatten(sys::ODESystem)
    if isempty(sys.systems)
        return sys
    else
        return ODESystem(equations(sys),
                         independent_variable(sys),
                         states(sys),
                         parameters(sys),
                         observed=observed(sys))
    end
end


using SymbolicUtils: Rewriters

function fixpoint_sub(x, dict)
    y = substitute(x, dict)
    while !isequal(x, y)
        y = x
        x = substitute(y, dict)
    end

    return x
end

function substitute_aliases(diffeqs, dict)
    lhss(diffeqs) .~ fixpoint_sub.(rhss(diffeqs), (dict,))
end

# Note that we reduce parameters, too
# i.e. `2param = 3` will be reduced away
isvar(s::Sym) = true
isvar(s::Term) = isvar(operation(s))
isvar(s::Any) = false

function get_α_x(αx)
    if isvar(αx)
        return 1, αx
    elseif istree(αx) && operation(αx) === (*)
        args = arguments(αx)
        nums = []
        syms = []
        for arg in args
            isvar(arg) ? push!(syms, arg) : push!(nums, arg)
        end

        if length(syms) == 1
            return prod(nums), syms[1]
        end
    else
        return nothing
    end
end

function alias_elimination(sys::ODESystem)
    eqs = vcat(equations(sys), observed(sys))
    subs = Pair[]
    diff_vars = filter(!isnothing, map(eqs) do eq
            if isdiffeq(eq)
                arguments(eq.lhs)[1]
            else
                nothing
            end
        end) |> Set

    # only substitute when the variable is algebraic
    del = Int[]
    for (i, eq) in enumerate(eqs)
        isdiffeq(eq) && continue
        res_left = get_α_x(eq.lhs)
        if !isnothing(res_left) && !(res_left[2] in diff_vars)
            # `α x = rhs` => `x = rhs / α`
            α, x = res_left
            push!(subs, x => _isone(α) ? eq.rhs : eq.rhs / α)
            push!(del, i)
        else
            res_right = get_α_x(eq.rhs)
            if !isnothing(res_right) && !(res_right[2] in diff_vars)
                # `lhs = β y` => `y = lhs / β`
                β, y = res_right
                push!(subs, y => _isone(β) ? eq.lhs : β * eq.lhs)
                push!(del, i)
            end
        end
    end
    deleteat!(eqs, del)

    eqs′ = substitute_aliases(eqs, Dict(subs))
    alias_vars = first.(subs)

    newstates = setdiff(states(sys), alias_vars)
    ODESystem(eqs′, sys.iv, newstates, parameters(sys), observed=alias_vars .~ last.(subs))
end

"""
$(SIGNATURES)

Use Kahn's algorithm to topologically sort observed equations.

Example:
```julia
julia> @variables t x(t) y(t) z(t) k(t)
(t, x(t), y(t), z(t), k(t))

julia> eqs = [
           x ~ y + z
           z ~ 2
           y ~ 2z + k
       ];

julia> ModelingToolkit.topsort_observed(eqs, [x, y, z, k])
3-element Vector{Equation}:
 Equation(z(t), 2)
 Equation(y(t), k(t) + 2z(t))
 Equation(x(t), y(t) + z(t))
```
"""
function topsort_observed(eqs, states)
    graph, assigns = observed2graph(eqs, states)
    neqs = length(eqs)
    degrees = zeros(Int, neqs)

    for 𝑠eq in 1:length(eqs); var = assigns[𝑠eq]
        for 𝑑eq in 𝑑neighbors(graph, var)
            # 𝑠eq => 𝑑eq
            degrees[𝑑eq] += 1
        end
    end

    q = Queue{Int}(neqs)
    for (i, d) in enumerate(degrees)
        d == 0 && enqueue!(q, i)
    end

    idx = 0
    ordered_eqs = similar(eqs)
    while !isempty(q)
        𝑠eq = dequeue!(q)
        ordered_eqs[idx+=1] = eqs[𝑠eq]
        var = assigns[𝑠eq]
        for 𝑑eq in 𝑑neighbors(graph, var)
            degree = degrees[𝑑eq] = degrees[𝑑eq] - 1
            degree == 0 && enqueue!(q, 𝑑eq)
        end
    end

    idx == neqs || throw(ArgumentError("The obversed equations have at least one cycle."))

    return ordered_eqs
end

function observed2graph(eqs, states)
    graph = BipartiteGraph(length(eqs), length(states))
    v2j = Dict(states .=> 1:length(states))

    # `assigns: eq -> var`, `eq` defines `var`
    assigns = similar(eqs, Int)

    for (i, eq) in enumerate(eqs)
        lhs_j = get(v2j, eq.lhs, nothing)
        lhs_j === nothing && throw(ArgumentError("The lhs $lhs of $eq, doesn't appear in states."))
        assigns[i] = lhs_j
        vs = vars(eq.rhs)
        for v in vs
            j = get(v2j, v, nothing)
            j !== nothing && add_edge!(graph, i, j)
        end
    end

    return graph, assigns
end
