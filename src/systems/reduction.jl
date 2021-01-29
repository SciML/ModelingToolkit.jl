export alias_elimination, flatten

function flatten(sys::ODESystem)
    if isempty(sys.systems)
        return sys
    else
        return ODESystem(equations(sys),
                         independent_variable(sys),
                         states(sys),
                         parameters(sys),
                         pins=pins(sys),
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
isvar(s) = s isa Sym ? true :
           istree(s) ? isvar(operation(s)) :
                       false

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

function is_univariate_expr(ex, iv)
    count = 0
    for var in vars(ex)
        if !isequal(iv, var) && !isparameter(var)
            count += 1
            count > 1 && return false
        end
    end
    return count <= 1
end

function is_sub_candidate(ex, iv, conservative)
    conservative || return true
    isvar(ex) || ex isa Number || is_univariate_expr(ex, iv)
end

function maybe_alias(lhs, rhs, diff_vars, iv, conservative)
    is_sub_candidate(rhs, iv, conservative) || return false, nothing

    res_left = get_α_x(lhs)
    if res_left !== nothing && !(res_left[2] in diff_vars)
        α, x = res_left
        sub = x => _isone(α) ? rhs : rhs / α
        return true, sub
    else
        return false, nothing
    end
end

function alias_elimination(sys::ODESystem; conservative=true)
    iv = independent_variable(sys)
    eqs = vcat(equations(sys), observed(sys))
    diff_vars = filter(!isnothing, map(eqs) do eq
            if isdiffeq(eq)
                arguments(eq.lhs)[1]
            else
                nothing
            end
        end) |> Set

    deps = Set()
    subs = Pair[]
    neweqs = Equation[]; sizehint!(neweqs, length(eqs))

    for (i, eq) in enumerate(eqs)
        # only substitute when the variable is algebraic
        if isdiffeq(eq)
            push!(neweqs, eq)
            continue
        end

        # `α x = rhs` => `x = rhs / α`
        ma, sub = maybe_alias(eq.lhs, eq.rhs, diff_vars, iv, conservative)
        if !ma
            # `lhs = β y` => `y = lhs / β`
            ma, sub = maybe_alias(eq.rhs, eq.lhs, diff_vars, iv, conservative)
        end

        isalias = false
        if ma
            l, r = sub
            # alias equations shouldn't introduce cycles
            if !(l in deps) && isempty(intersect(deps, vars(r)))
                push!(deps, l)
                push!(subs, sub)
                isalias = true
            end
        end

        if !isalias
            neweq = _iszero(eq.lhs) ? eq : 0 ~ eq.rhs - eq.lhs
            push!(neweqs, neweq)
        end
    end

    alias_vars = first.(subs)
    sys_states = states(sys)
    alias_eqs = topsort_equations(alias_vars .~ last.(subs), sys_states)

    newstates = setdiff(sys_states, alias_vars)
    ODESystem(neweqs, sys.iv, newstates, parameters(sys), observed=alias_eqs)
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

julia> ModelingToolkit.topsort_equations(eqs, [x, y, z, k])
3-element Vector{Equation}:
 Equation(z(t), 2)
 Equation(y(t), k(t) + 2z(t))
 Equation(x(t), y(t) + z(t))
```
"""
function topsort_equations(eqs, states; check=true)
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
    ordered_eqs = similar(eqs, 0); sizehint!(ordered_eqs, neqs)
    while !isempty(q)
        𝑠eq = dequeue!(q)
        idx+=1
        push!(ordered_eqs, eqs[𝑠eq])
        var = assigns[𝑠eq]
        for 𝑑eq in 𝑑neighbors(graph, var)
            degree = degrees[𝑑eq] = degrees[𝑑eq] - 1
            degree == 0 && enqueue!(q, 𝑑eq)
        end
    end

    (check && idx != neqs) && throw(ArgumentError("The equations have at least one cycle."))

    return ordered_eqs
end

function observed2graph(eqs, states)
    graph = BipartiteGraph(length(eqs), length(states))
    v2j = Dict(states .=> 1:length(states))

    # `assigns: eq -> var`, `eq` defines `var`
    assigns = similar(eqs, Int)

    for (i, eq) in enumerate(eqs)
        lhs_j = get(v2j, eq.lhs, nothing)
        lhs_j === nothing && throw(ArgumentError("The lhs $(eq.lhs) of $eq, doesn't appear in states."))
        assigns[i] = lhs_j
        vs = vars(eq.rhs)
        for v in vs
            j = get(v2j, v, nothing)
            j !== nothing && add_edge!(graph, i, j)
        end
    end

    return graph, assigns
end
