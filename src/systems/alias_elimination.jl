using SymbolicUtils: Rewriters

const KEEP = typemin(Int)

function alias_eliminate_graph(sys)
    sys = flatten(sys)
    s = get_structure(sys)
    if !(s isa SystemStructure)
        sys = initialize_system_structure(sys)
        s = structure(sys)
    end

    @unpack graph, varassoc = s

    is_linear_equations, eadj, cadj = find_linear_equations(sys)
    old_cadj = map(copy, cadj)

    is_not_potential_state = iszero.(varassoc)
    is_linear_variables = copy(is_not_potential_state)
    for i in 𝑠edges(graph); is_linear_equations[i] || continue
        for j in 𝑠vertices(graph, i)
            is_linear_variables[j] = false
        end
    end
    solvable_variables = findall(is_linear_variables)

    linear_equations = findall(is_linear_equations)


    rank1 = bareiss!(
        (eadg, cadj),
        old_cadj, linear_equations, is_linear_variables, 1
       )

    v_solved = [eadg[i][1] for i in 1:rank1]
    v_eliminated = setdiff(solvable_variables, v_solved)
    n_null_vars = length(v_eliminated)

    v_types = fill(KEEP, ndsts(graph))
    for v in v_eliminated
        v_types[v] = 0
    end

    rank2 = bareiss!(
        (eadg, cadj),
        old_cadj, linear_equations, is_not_potential_state, rank1+1
       )

    rank3 = bareiss!(
        (eadg, cadj),
        old_cadj, linear_equations, nothing, rank2+1
       )

    # kind of like the backward substitution
    for ei in reverse(1:rank2)
        locally_structure_simplify!(
                                    (eadg[ei], cadj[ei]),
                                    invvarassoc, v_eliminated, v_types
                                   )
    end

    reduced = false
    for ei in 1:rank2
        if length(cadj[ei]) > length(old_cadj[ei])
            cadj[ei] = old_cadj[ei]
        else
            cadj[ei] = eadg[linear_equations[ei]]
            reduced |= locally_structure_simplify!(
                                                   (eadg[ei], cadj[ei]),
                                                   invvarassoc, v_eliminated, v_types
                                                  )
        end
    end

    while reduced
        for ei in 1:rank2
            if !isempty(eadg[ei])
                reduced |= locally_structure_simplify!(
                                                       (eadg[ei], cadj[ei]),
                                                       invvarassoc, v_eliminated, v_types
                                                      )
                reduced && break # go back to the begining of equations
            end
        end
    end

    for ei in rank2+1:length(linear_equations)
        eadg[ei] = old_cadj[ei]
    end

    for (ei, e) in enumerate(linear_equations)
        graph.eadglist[e] = eadg[ei]
    end

    degenerate_equations = rank3 < length(linear_equations) ? linear_equations[rank3+1:end] : Int[]
    return v_eliminated, v_types, n_null_vars, degenerate_equations
end

iszeroterm(v_types, v) = v_types[v] == 0
isirreducible(v_types, v) = v_types[v] == KEEP
isalias(v_types, v) = v_types[v] > 0 && !isirreducible(v_types, v)
alias(v_types, v) = v_types[v]
negalias(v_types, v) = -v_types[v]

function locally_structure_simplify!(
        (vars, coeffs),
        invvarassoc, v_eliminated, v_types
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
            push!(v_eliminated, v)
            v_types[v] = 0
            empty!(vars); empty!(coeffs)
            return true
        elseif length(vars) == 2 && abs(coeffs[1]) == abs(coeffs[2])
            if (coeffs[1] > 0 && coeffs[2] < 0) || (coeffs[1] < 0 && coeffs[2] > 0)
                # positive alias
                push!(v_eliminated, v)
                v_types[v] = vars[2]
            else
                # negative alias
                push!(v_eliminated, v)
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
        (eadg, cadj),
        old_cadj, linear_equations, is_linear_variables, offset
       )
    m = nsrcs(solvable_graph)
    # v = eadg[ei][vj]
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
            pivot = cadj[ei][vj]
            deleteat!(cadj[ei] , vj)
            v = eadg[ei][vj]
            deleteat!(eadg[ei], vj)
            if ei != k
                swap!(cadj, ei, k)
                swap!(old_cadj, ei, k)
                swap!(eadg, ei, k)
                swap!(linear_equations, ei, k)
            end
        else # rank deficient
            return k-1
        end

        for ei in k+1
            # elimate `v`
            coeff = 0
            vars = eadg[ei]
            vj = findfirst(isequal(v), vars)
            if vj === nothing # `v` is not in in `e`
                continue
            else # remove `v`
                coeff = cadj[ei][vj]
                deleteat!(cadj[ei], vj)
                deleteat!(eadg[ei], vj)
            end

            # the pivot row
            kvars = eadg[k]
            kcoeffs = cadj[k]
            # the elimination target
            ivars = eadg[ei]
            icoeffs = cadj[ei]

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

            eadg[ei], tmp_incidence = tmp_incidence, eadg[ei]
            cadj[ei], tmp_coeffs = tmp_coeffs, cadj[ei]
        end
        last_pivot = pivot
        # add `v` in the front of the `k`-th equation
        pushfirst!(eadg[k], v)
        pushfirst!(cadj[k], pivot)
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
        mask,
        constraint,
    )
    for i in range
        vertices = 𝑠vertices(solvable_graph, i)
        if constraint(length(vertices))
            for (j, v) in enumerate(vertices)
                (mask === nothing || mask[v]) && return i, j
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
