using SymbolicUtils: Rewriters

function fixpoint_sub(x, dict)
    y = substitute(x, dict)
    while !isequal(x, y)
        y = x
        x = substitute(y, dict)
    end

    return x
end

function substitute_aliases(eqs, dict)
    sub = Base.Fix2(fixpoint_sub, dict)
    map(eq->eq.lhs ~ sub(eq.rhs), eqs)
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

function alias_elimination(sys)
    sys = flatten(sys)
    s = get_structure(sys)
    if !(s isa SystemStructure)
        sys = initialize_system_structure(sys)
        s = structure(sys)
    end
    iv = independent_variable(sys)
    eqs = equations(sys)
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
    sts = states(sys)
    fullsts = vcat(map(eq->eq.lhs, observed(sys)), sts, parameters(sys))
    alias_eqs = topsort_equations(alias_vars .~ last.(subs), fullsts)
    newstates = setdiff(sts, alias_vars)

    @set! sys.eqs = substitute_aliases(neweqs, Dict(subs))
    @set! sys.states = newstates
    @set! sys.observed = [observed(sys); alias_eqs]
    return 
end


function alias_elimination_2(sys)
    sys = flatten(sys)
    s = get_structure(sys)
    if !(s isa SystemStructure)
        sys = initialize_system_structure(sys)
        s = structure(sys)
    end
    find_solvables!(sys)
    @unpack graph, solvable_graph, is_linear_equations, varassoc = s

    is_not_potential_state = iszero.(varassoc)
    is_linear_variables = copy(is_not_potential_state)
    for i in 𝑠edges(graph); is_linear_equations[i] || continue
        for j in 𝑠vertices(graph, i)
            is_linear_variables[j] = false
        end
    end
    solvable_variables = findall(is_linear_variables)

    linear_equations = findall(is_linear_equations)

    offset = 1
    coeffs = solvable_graph.metadata
    old_coeffs = map(copy, coeffs)
    fadj = solvable_graph.fadjlist

    rank1 = bareiss!(
        (fadj, coeffs),
        old_coeffs, linear_equations, is_linear_variables, offset
       )

    v_solved = [fadj[i][1] for i in 1:rank1]
    v_null = setdiff(solvable_variables, v_solved)
    n_null_vars = length(v_null)

    v_types = fill(KEEP, ndsts(graph))
    for v in v_null
        v_types[v] = 0
    end

    rank2 = bareiss!(
        (fadj, coeffs),
        old_coeffs, linear_equations, is_not_potential_state, offset
       )
end

iszeroterm(v_types, v) = v_types[v] == 0
isirreducible(v_types, v) = v_types[v] == KEEP
isalias(v_types, v) = v_types[v] > 0 && !isirreducible(v_types, v)
alias(v_types, v) = v_types[v]
negalias(v_types, v) = -v_types[v]

function locally_structure_simplify!(
        (vars, coeffs),
        invvarassoc, v_null, v_types
       )
    while length(vars) > 1 && any(!isequal(KEEP), (v_types[v] in @view vars[2:end]))
        for vj in 2:length(vars)
            v = vars[vj]
            if isirreducible(v_types, v)
                continue
            elseif iszeroterm(v_types, v)
                deleteat!(vars, vj)
                deleteat!(coeffs, vj)
                break
            else
                coeff = coeffs[vj]
                if isalias(v_types, v)
                    v = alias(v_types, v)
                else
                    v = negalias(v_types, v)
                    coeff = -coeff
                end

                has_v = false
                for vi in 2:length(vars)
                    (vi !== vj && vars[vi] == v) || continue
                    has_v = true
                    c = (coeffs[vi] += coeff)
                    if c == 0
                        if vi < vj
                            deleteat!(vars, [vi, vj])
                            deleteat!(coeffs, [vi, vj])
                        else
                            deleteat!(vars, [vj, vi])
                            deleteat!(coeffs, [vj, vi])
                        end
                    end
                    break
                end # for vi

                if has_v
                    break
                else
                    vars[vj] = v
                    coeffs[vj] = coeff
                end # if
            end # else
        end # for
    end # while

    v = first(vars)
    if invvarassoc[v] == 0
        if length(nvars) == 1
            push!(v_null, v)
            v_types[v] = 0
            empty!(vars); empty!(coeffs)
            return true
        elseif length(vars) == 2 && abs(coeffs[1]) == abs(coeffs[2])
            if (coeffs[1] > 0 && coeffs[2] < 0) || (coeffs[1] < 0 && coeffs[2] > 0)
                # positive alias
                push!(v_null, v)
                v_types[v] = vars[2]
            else
                # negative alias
                push!(v_null, v)
                v_types[v] = -vars[2]
            end
            empty!(vars); empty!(coeffs)
            return true
        end
    end
    return false
end

"""
$(SIGNATURES)

Use Bareiss algorithm to compute the nullspace of an integer matrix exactly.
"""
function bareiss!(
        (fadj, coeffs),
        old_coeffs, linear_equations, is_linear_variables, offset
       )
    m = nsrcs(solvable_graph)
    # v = fadj[ei][vj]
    v = ei = vj = 0
    pivot = last_pivot = 1
    tmp_incidence = Int[]
    tmp_coeffs = Int[]
    vars = Set{Int}()

    # j -> vj
    # e -> ei
    # vj -> v
    # GcInt2 modified

    for k in offset:m
        ###
        ### Pivoting:
        ###
        ei, vj = find_first_linear_variable(solvable_graph, k:m, is_linear_variables, isequal(1))
        if vj == 0
            ei, vj = find_first_linear_variable(solvable_graph, k:m, is_linear_variables, isequal(2))
        else
            ei, vj = find_first_linear_variable(solvable_graph, k:m, is_linear_variables, _->true)
        end

        if vj > 0 # has a pivot
            pivot = coeffs[ei][vj]
            deleteat!(coeffs[ei] , vj)
            v = fadj[ei][vj]
            deleteat!(fadj[ei], vj)
            if ei != k
                swap!(coeffs, ei, k)
                swap!(old_coeffs, ei, k)
                swap!(fadj, ei, k)
                swap!(linear_equations, ei, k)
            end
        else # rank deficient
            return k-1
        end

        for ei in k+1
            # elimate `v`
            coeff = 0
            vars = fadj[ei]
            vj = findfirst(isequal(v), vars)
            if vj === nothing # `v` is not in in `e`
                continue
            else # remove `v`
                coeff = coeffs[ei][vj]
                deleteat!(coeffs[ei], vj)
                deleteat!(fadj[ei], vj)
            end

            # the pivot row
            kvars = fadj[k]
            kcoeffs = coeffs[k]
            # the elimination target
            ivars = fadj[ei]
            icoeffs = coeffs[ei]

            empty!(tmp_incidence)
            empty!(tmp_coeffs)
            empty!(vars)
            union!(vars, kvars, ivars)

            for v in vars
                ck = getcoeff(kvars, kcoeffs, v)
                ci = getcoeff(ivars, icoeffs, v)
                ci = (pivot*ci - coeff*ck) ÷ last_pivot
                if ci !== 0
                    push!(tmp_incidence, v)
                    push!(tmp_coeffs, ci)
                end
            end

            fadj[ei], tmp_incidence = tmp_incidence, fadj[ei]
            coeffs[ei], tmp_coeffs = tmp_coeffs, coeffs[ei]
        end
        last_pivot = pivot
        # add `v` in the front of the `k`-th equation
        pushfirst!(fadj[k], v)
        pushfirst!(coeffs[k], pivot)
    end

    return m # fully ranked
end

swap!(v, i, j) = ((v[i], v[j] = v[j], v[i]); nothing)

function getcoeff(vars, coeffs, var)
    for (vj, v) in enumerate(vars)
        v == var && return coeffs[vj]
    end
    return 0
end

"""
$(SIGNATURES)

Find the first linear variable such that `𝑠vertices(adj, i)[j]` is true given
the `constraint`.
"""
@inline function find_first_linear_variable(
        solvable_graph,
        range,
        is_linear_variables,
        constraint,
    )
    for i in range
        vertices = 𝑠vertices(solvable_graph, i)
        if constraint(length(vertices))
            for (j, v) in enumerate(vertices)
                is_linear_variables[v] && return i, j
            end
        end
    end
    return 0, 0
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
