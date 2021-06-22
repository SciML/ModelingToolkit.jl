Base.@kwdef mutable struct TearingState
    td::TearingSetup
    e_solved::Vector{Int} = Int[]
    v_solved::Vector{Int} = Int[]
    e_residual::Vector{Int} = Int[]
    v_residual::Vector{Int} = Int[]
    der_v_solved::Vector{Int} = Int[]
    der_e_solved::Vector{Int} = Int[]
    der_v_tear::Vector{Int} = Int[]
    vc::Vector{Int} = Int[]
    ec::Vector{Int} = Int[]
    vc_candidates::Vector{Vector{Int}} = [Int[] for _ in 1:6]
    inspected::BitVector = falses(0)
end

function state_selection!(sys; kwargs...)
    s = get_structure(sys)
    if !(s isa SystemStructure)
        sys = initialize_system_structure(sys)
        s = structure(sys)
    end
    check_consistency(s)
    sys, assign, eqassoc = pantelides!(sys; kwargs...)
    s = get_structure(sys)
    eqs = equations(sys)
    eqs = [eqs; Vector{Equation}(undef, length(eqassoc) - length(eqs))]
    iv = independent_variable(sys)
    for (i, a) in enumerate(eqassoc); a > 0 || continue
        eqs[a] = derivative(eqs[i].lhs, iv) ~ derivative(eqs[i].rhs, iv)
    end
    @set! sys.eqs = eqs

    fullvars = [s.fullvars; Vector{Any}(undef, length(s.varassoc) - length(s.fullvars))]
    for (i, a) in enumerate(s.varassoc); a > 0 || continue
        fullvars[a] = value(derivative(fullvars[i], iv))
    end
    @set! s.fullvars = fullvars

    # TODO: use eqmask
    # When writing updating a graph we often have two choices:
    #   1. Construct a new graph with updated vertices and edges;
    #   2. Make all the functions take masks for equations and variables.
    #
    # We want to use the second choice as much as possible as it minimizes
    # memory footprint while potentially strengthen temporal locality.
    #
    # N.B.: assign changes meaning here. In Pantelides, `assign` is w.r.t. to
    # all the equations and variables with the highest derivative. Here,
    # `assign` is w.r.t. equations and variables with the highest derivative.
    highest_order_graph = BipartiteGraph(0, length(assign), Val(false))
    eq_reidx = Int[]
    for (i, es) in enumerate(eqassoc); es == 0 || continue
        vars = s.graph.fadjlist[i]
        # N.B.: this is an alias to the original graph
        push!(highest_order_graph.fadjlist, vars)
        highest_order_graph.ne += length(vars)
        push!(eq_reidx, i)
    end
    # matching is read-only, so aliasing is fine
    assign = matching(highest_order_graph, s.varmask)

    # Compute SCC and map the indices back to the original system
    scc = Vector{Int}[]
    for component in find_scc(highest_order_graph, assign)
        push!(scc, eq_reidx[component])
    end
    for (j, a) in enumerate(assign); a == UNASSIGNED && continue
        assign[j] = eq_reidx[a]
    end

    @set! s.inv_varassoc = inverse_mapping(s.varassoc)
    @set! s.inv_eqassoc = inverse_mapping(s.eqassoc)
    @set! s.inv_assign = inverse_mapping(assign)
    @set! s.assign = assign
    @set! s.scc = scc
    @set! sys.structure = s

    # Otter and Elmqvist (2017) 4.3 [1, 2]:
    for (i, (a, ia)) in enumerate(zip(s.varassoc, s.inv_varassoc))
        # lower differentiated and algebraic variables are potentially states
        s.varmask[i] = a > 0 || a == 0 == ia
    end

    eq_constraint_set = Vector{Vector{Int}}[]
    var_constraint_set = Vector{Vector{Int}}[]

    @unpack graph = s
    unknown_subgraph = BipartiteGraph(nsrcs(graph), ndsts(graph), Val(false))
    original_vars = trues(ndsts(graph))

    # Optimization
    tmp_var_constaint_set = falses(ndsts(graph))
    for c in scc
        check_component(c, eqassoc) || break # skip lower order equations
        eq_constraints, var_constraints = sorted_constraints(c, s)
        push!(eq_constraint_set, eq_constraints)
        push!(var_constraint_set, var_constraints)

        for (i, eqs) in enumerate(eq_constraints)
            for v in var_constraints[i]
                tmp_var_constaint_set[v] = true
            end
            for eq in eqs
                for var in 𝑠neighbors(graph, eq); tmp_var_constaint_set[var] || continue
                    add_edge!(unknown_subgraph, eq, var)
                end

                if s.inv_eqassoc[eq] != 0 # differentiated equations
                    for var in 𝑠neighbors(graph, eq)
                        original_vars[var] = false
                    end
                end
            end
            for v in var_constraints[i]
                tmp_var_constaint_set[v] = false
            end
        end
    end
    ts = TearingState(td = TearingSetup(unknown_subgraph.fadjlist, length(assign)))

    higher_der(vs, assoc) = map(vs) do v
        dv = assoc[v]
        dv < 1 && scc_throw()
        dv
    end
    issolvable = let eqs=equations(sys), fullvars=s.fullvars
        (e, v) -> linear_expansion(eqs[e], fullvars[v])[3] # islinear
    end
    islineareq = let eqs=equations(sys), fullvars=s.fullvars
        (e, v) -> begin
            a, _, islinear = linear_expansion(eqs[e], fullvars[v])
            return islinear, !(a isa Symbolic)
        end
    end

    obs = similar(eqs)
    empty!(s.partitions)

    for (bk, constraints) in enumerate(zip(eq_constraint_set, var_constraint_set))
        highest_order = length(first(constraints))
        println("=============== BLT block $bk ===============")
        for (order, (eq_constraint, var_constraint)) in enumerate(Iterators.reverse(zip(constraints...)))
            println("--------------- BLT: $bk.$order ---------------")
            lower_order = order < highest_order
            @info "" equations=eqs[eq_constraint] unknows=s.fullvars[var_constraint]
            ###
            ### Try to solve the scalar equation via a linear solver
            ###
            if length(eq_constraint) == 1 == length(var_constraint) && solve_equations!(obs, sys, eq_constraint, var_constraint)
                @info "Scalar equation solved"
                v = var_constraint[1]
                # We can reduce the state if it doesn't have the highest order, or it's algebraic
                if lower_order || isalgebraic(s, v)
                    @assert s.varmask[v]
                    s.varmask[v] = false
                end
                continue
            end

            ###
            ### Try tearing
            ###
            length(var_constraint) >= length(eq_constraint) || scc_throw()
            if order == 1
                empty!(ts.der_e_solved)
                empty!(ts.der_v_solved)
                empty!(ts.der_v_tear)

                ts.ec = eq_constraint
                ts.vc = var_constraint
            else
                ts.der_e_solved = higher_der(ts.e_solved, eqassoc)
                ts.der_v_solved = higher_der(ts.v_solved, s.varassoc)
                ts.der_v_tear = higher_der(ts.v_residual, s.varassoc)

                ts.ec = setdiff(eq_constraint, ts.der_e_solved)
                ts.vc = setdiff(var_constraint, ts.der_v_solved)
            end

            #@show tearing_with_candidates!(ts, issolvable, s.fullvars)

            for v in ts.v_solved; (lower_order || isalgebraic(s, v)) || continue
                @assert s.varmask[v]
                s.varmask[v] = false
            end

            @unpack e_solved, v_solved, e_residual, v_residual = ts
            if isempty(e_residual)
                @assert solve_equations!(obs, sys, e_solved, v_solved)
                push!(s.partitions, SystemPartition(;e_solved, v_solved, e_residual, v_residual))
            elseif length(eq_constraint) == length(var_constraint)
                # Use a linear system of equations solver to further reduce the
                # the number of states.
                # TODO: solving linear system of equations
                islinear, isconstant = islinearsystem(s, ts, islineareq, eq_constraint, var_constraint)
                @show order == highest_order
                @show islinear, isconstant
                if islinear && isconstant #TODO: solve the equations at runtime
                    @views obs[eq_constraint] .= s.fullvars[var_constraint] .~ solve_for(eqs[eq_constraint], s.fullvars[var_constraint], check=false)
                    for v in v_residual; (lower_order || isalgebraic(s, v)) || continue
                        @show s.fullvars[v]
                        @assert s.varmask[v]
                        s.varmask[v] = false
                    end
                    push!(s.partitions, SystemPartition(;e_solved=eq_constraint, v_solved=var_constraint, e_residual=Int[], v_residual=Int[]))
                else
                    @show s.varmask[v_residual]
                    push!(s.partitions, SystemPartition(;e_solved, v_solved, e_residual, v_residual))
                    @assert solve_equations!(obs, sys, e_solved, v_solved)
                    @views obs[e_residual] .= eqs[e_residual]
                end
            #elseif length(eq_constraint) < length(var_constraint) && order < highest_order
            #    error("ModelingToolkit cannot handle this kind of system yet. " *
            #            "Please file an issue with a reproducible example if " *
            #            "you encounter this error message.")
            #else
            #    scc_throw()
            end
        end
        println("=============== End of BLT block $bk ===============")
    end

    @info "" obs

    new_states = []
    lhss = Set{Int}()
    for (v, m) in enumerate(s.varmask); m || continue
        push!(new_states, diff2term(s.fullvars[v]))
        dv = s.varassoc[v]
        dv > 0 && push!(lhss, dv)
    end
    @info "New states" new_states

    ###
    ### Substitute reduced equations
    ###
    assign = matching(s.graph, .!s.varmask)
    @set! s.assign = assign
    @set! s.inv_assign = inverse_mapping(assign)
    var2idx = Dict{Any,Int}(v => i for (i, v) in enumerate(s.fullvars))

    sorted_eqs = sizehint!(Equation[], nsrcs(s.graph))
    new_scc = find_scc(s.graph, assign)
    new_eqs = similar(eqs, length(new_states))
    obs_dict = Dict()
    obs_eqs = Equation[]
    i = 0
    for c in new_scc
        for e in c
            #=
            vv = s.inv_assign[e]
            if vv > 0
                @info "Assigned" fullvars[vv] obs[e]
            else
                @warn "Unassigned" obs[e]
            end
            =#
            eq = obs[e]
            push!(sorted_eqs, eq)
            iv = _iszero(eq.lhs) ? -1 : var2idx[eq.lhs]
            if iv == -1 || iv in lhss # if it's algebraic or differential
                @show eq
                new_eqs[i+=1] = obs[e]
                while iv > 0 && (iv = s.inv_varassoc[iv]) > 0
                    s.inv_varassoc[iv] > 0 || break
                    nv = diff2term(s.fullvars[iv])
                    if iter == 1
                        new_eqs[i] = substitute(new_eqs[i].lhs, Dict(s.fullvars[iv] => nv)) ~ new_eqs[i].rhs
                    end
                    @show new_eqs[i+=1] = s.fullvars[iv] ~ nv
                end
            else
                obs_dict[eq.lhs] = eq.rhs
                push!(obs_eqs, eq)
            end
        end
    end
    @info "New equations" new_eqs
    @show length(new_states)
    @show count(.!s.varmask)

    @set! sys.observed = obs_eqs
    @set! sys.states = new_states
    @set! sys.eqs = map(Base.Fix2(fixpoint_sub_exclude_differential, obs_dict), new_eqs)
    @set! sys.structure = s
    return sys
end

isalgebraic(s, v) = s.varassoc[v] == 0 == s.inv_varassoc[v]

"""
$(TYPEDSIGNATURES)

Like `substitute`, but stop substitution when `Differential` is encountered.

# Example:
```julia
julia> @variables t x(t); D = Differential(t)
(::Differential) (generic function with 2 methods)

julia> StructuralTransformations.substitute_exclude_differential((D(D(x)) + D(x)).val, Dict(D(x).val => t.val))
t + Differential(t)(Differential(t)(x(t)))
```
Note that `substitute_exclude_differential` is required to make sure `D(x)` and
`D(D(x))` are treated like different variables.
"""
substitute_exclude_differential(eq::Equation, dict) = substitute_exclude_differential(eq.lhs, dict) ~ substitute_exclude_differential(eq.rhs, dict)
function substitute_exclude_differential(expr, dict)
    haskey(dict, expr) && return dict[expr]
    istree(expr) || return expr

    op = operation(expr)
    op isa Differential && return expr
    canfold = true
    args = map(arguments(expr)) do x
        x′ = substitute_exclude_differential(x, dict)
        canfold = canfold && !(x′ isa Symbolic)
        x′
    end
    canfold ? op(args...) : similarterm(expr, op, args, symtype(expr))
end
function fixpoint_sub_exclude_differential(x, dict)
    y = substitute_exclude_differential(x, dict)
    while !isequal(x, y)
        y = x
        x = substitute_exclude_differential(y, dict)
    end

    return x
end

function is_original_linear(s, islineareq, e, v)
    maxiter = 8000
    while s.inv_eqassoc[e] > 0
        e = s.inv_eqassoc[e]
        v = s.inv_varassoc[v]
        maxiter -= 1
        maxiter <= 0 && scc_throw()
    end
    return islineareq(e, v)
end

function islinearsystem(s, ts, islineareq, eqs, vars)
    @unpack inspected = ts
    resize!(inspected, length(s.fullvars))
    for v in vars
        inspected[v] = true
    end

    islinear = true
    isconstant = true
    for e in eqs, v in 𝑠neighbors(s.graph, e); inspected[v] || continue # Optimization: avoiding extra work
        islinear′, isconstant′ = is_original_linear(s, islineareq, e, v)
        islinear &= islinear′
        islinear || break
        isconstant &= isconstant′
    end

    for v in vars
        inspected[v] = false
    end
    return islinear, isconstant
end

function tearing_with_candidates!(ts, issolvable, fullvars)
    @unpack vc_candidates = ts
    foreach(x->empty!(x), vc_candidates)

    vc_candidates[1] = ts.der_v_tear

    for v in setdiff(ts.vc, ts.der_v_tear)
        # TODO
        #push!(hasmetadata(fullvars[v], VariableDefaultValue) ? vc_candidates[2] : vc_candidates[4], v)
        push!(vc_candidates[2], v)
    end

    ts.e_solved, ts.v_solved, ts.e_residual, ts.v_residual = tearEquations!(
        ts.td, issolvable, ts.ec, ts.vc_candidates;
        eSolvedFixed=ts.der_e_solved,
        vSolvedFixed=ts.der_v_solved
    )
end

@noinline scc_throw() = error("Internal error: " *
                              "Strongly connected component invariant is violated!\n" *
                              "Please file a bug report with full stack trace " *
                              "and a reproducer.")
function check_component(c, eqassoc) # return `is_highest_order`
    # Fast paths
    isempty(c) && return true
    length(c) == 1 && return eqassoc[c[1]] == 0

    # If there are more than one equation in the component, it must contain
    # highest order equations only.
    for eq in c; eqassoc[eq] == 0 && continue
        scc_throw()
    end
    return true
end

"""
$(TYPEDSIGNATURES)

The constraints are sorted in descending order.
"""
function sorted_constraints(c::Vector{Int}, s::SystemStructure)
    @unpack inv_varassoc, inv_eqassoc, inv_assign = s
    # The set of variables that this compoent solves for
    unknowns = [inv_assign[eq] for eq in c]
    eq_constraints = [c]
    var_constraints = [unknowns]
    while true
        lower_eqs = Int[]
        for eq in eq_constraints[end]
            lower_eq = inv_eqassoc[eq]
            lower_eq > 0 && push!(lower_eqs, lower_eq)
        end
        isempty(lower_eqs) && break

        push!(eq_constraints, lower_eqs)

        vars = Int[]
        for var in var_constraints[end]
            lower_var = inv_varassoc[var]
            lower_var > 0 && push!(vars, lower_var)
        end
        push!(var_constraints, vars)

        length(vars) < length(lower_eqs) && scc_throw()
    end
    return eq_constraints, var_constraints
end

function solve_equations!(obs, sys, ieqs, ivars)
    eqs = equations(sys)
    s = get_structure(sys)
    solved = true
    for (ie, iv) in zip(ieqs, ivars)
        var = s.fullvars[iv]
        eq = eqs[ie]
        sol = solve_for(eq, var; check=false)
        if sol === nothing
            solved = false
            obs[ie] = eq
        else
            obs[ie] = var ~ sol
        end
    end
    return solved
end

#=
using ModelingToolkit
@parameters t u[1:8](t)
@variables x[1:8](t)
D = Differential(t)
eqs = [
    0 ~ u[1] + x[1] - x[2]
    0 ~ u[2] + x[1] + x[2] - x[3] + D(x[6])
    0 ~ u[3] + x[1] + D(x[3]) - x[4]
    0 ~ u[4] + 2*(D^2)(x[1]) + (D^2)(x[2]) + (D^2)(x[3]) + D(x[4]) + (D^3)(x[6])
    0 ~ u[5] + 3*(D^2)(x[1]) + 2*(D^2)(x[2]) + x[5] + 0.1x[8]
    0 ~ u[6] + 2x[6] + x[7]
    0 ~ u[7] + 3x[6] + 4x[7]
    0 ~ u[8] + x[8] - sin(x[8])
]
sys = ODESystem(eqs, t)
sys = initialize_system_structure(sys); StructuralTransformations.state_selection!(sys)

# TODO: detect this case and warn that this can only be a DAEProblem
using ModelingToolkit
@parameters t u[1:8](t)
@variables x[1:8](t)
D = Differential(t)
eqs = [
    0 ~ sin(D(x[1])) + u[1]
]
sys = ODESystem(eqs, t)
sys = initialize_system_structure(sys); StructuralTransformations.state_selection!(sys)

using ModelingToolkit
@parameters t u[1:8](t)
@variables x[1:8](t) xx1(t) xx2(t) xx3(t) xx6(t) xxx6(t)
D = Differential(t)
eqs = [
    D(x[1]) ~ xx1
    D(x[2]) ~ xx2
    D(x[3]) ~ xx3
    D(x[6]) ~ xx6
    D(xx6) ~ xxx6
    0 ~ u[1] + x[1] - x[2]
    0 ~ u[2] + x[1] + x[2] - x[3] + xx6
    0 ~ u[3] + x[1] + xx3 - x[4]
    0 ~ u[4] + 2*D(xx1) + D(xx2) + D(xx3) + D(x[4]) + D(xxx6)
    0 ~ u[5] + 3*D(xx1) + 2*D(xx2) + x[5] + 0.1x[8]
    0 ~ u[6] + 2x[6] + x[7]
    0 ~ u[7] + 3x[6] + 4x[7]
    0 ~ u[8] + x[8] - sin(x[8])
]
sys = ODESystem(eqs, t)#, [x..., xx1, xx2, xx3, xx6, xxx6], u)
sys = initialize_system_structure(sys); StructuralTransformations.state_selection!(sys)
=#
