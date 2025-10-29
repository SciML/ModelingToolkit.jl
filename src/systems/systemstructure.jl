using DataStructures
using Symbolics: linear_expansion, unwrap
using SymbolicUtils: iscall, operation, arguments, Symbolic
using SymbolicUtils: quick_cancel, maketerm
using ..ModelingToolkit
import ..ModelingToolkit: isdiffeq, var_from_nested_derivative, vars!, flatten,
                          value, InvalidSystemException, isdifferential, _iszero,
                          isparameter, Connection,
                          independent_variables, SparseMatrixCLIL, AbstractSystem,
                          equations, isirreducible, input_timedomain, TimeDomain,
                          InferredTimeDomain,
                          VariableType, getvariabletype, has_equations, System
using ..BipartiteGraphs
import ..BipartiteGraphs: invview, complete
using Graphs
using UnPack
using Setfield
using SparseArrays

function quick_cancel_expr(expr)
    Rewriters.Postwalk(quick_cancel,
        similarterm = (x, f, args;
            kws...) -> maketerm(typeof(x), f, args,
            SymbolicUtils.metadata(x),
            kws...))(expr)
end

export SystemStructure, TransformationState, TearingState, mtkcompile!
export isdiffvar, isdervar, isalgvar, isdiffeq, algeqs, is_only_discrete
export dervars_range, diffvars_range, algvars_range
export DiffGraph, complete!
export get_fullvars, system_subset

struct DiffGraph <: Graphs.AbstractGraph{Int}
    primal_to_diff::Vector{Union{Int, Nothing}}
    diff_to_primal::Union{Nothing, Vector{Union{Int, Nothing}}}
end

DiffGraph(primal_to_diff::Vector{Union{Int, Nothing}}) = DiffGraph(primal_to_diff, nothing)
function DiffGraph(n::Integer, with_badj::Bool = false)
    DiffGraph(Union{Int, Nothing}[nothing for _ in 1:n],
        with_badj ? Union{Int, Nothing}[nothing for _ in 1:n] : nothing)
end

function Base.copy(dg::DiffGraph)
    DiffGraph(copy(dg.primal_to_diff),
        dg.diff_to_primal === nothing ? nothing : copy(dg.diff_to_primal))
end

@noinline function require_complete(dg::DiffGraph)
    dg.diff_to_primal === nothing &&
        error("Not complete. Run `complete` first.")
end

Graphs.is_directed(dg::DiffGraph) = true
function Graphs.edges(dg::DiffGraph)
    (i => v for (i, v) in enumerate(dg.primal_to_diff) if v !== nothing)
end
Graphs.nv(dg::DiffGraph) = length(dg.primal_to_diff)
Graphs.ne(dg::DiffGraph) = count(x -> x !== nothing, dg.primal_to_diff)
Graphs.vertices(dg::DiffGraph) = Base.OneTo(nv(dg))
function Graphs.outneighbors(dg::DiffGraph, var::Integer)
    diff = dg.primal_to_diff[var]
    return diff === nothing ? () : (diff,)
end
function Graphs.inneighbors(dg::DiffGraph, var::Integer)
    require_complete(dg)
    diff = dg.diff_to_primal[var]
    return diff === nothing ? () : (diff,)
end
function Graphs.add_vertex!(dg::DiffGraph)
    push!(dg.primal_to_diff, nothing)
    if dg.diff_to_primal !== nothing
        push!(dg.diff_to_primal, nothing)
    end
    return length(dg.primal_to_diff)
end

function Graphs.add_edge!(dg::DiffGraph, var::Integer, diff::Integer)
    dg[var] = diff
end

# Also pass through the array interface for ease of use
Base.:(==)(dg::DiffGraph, v::AbstractVector) = dg.primal_to_diff == v
Base.:(==)(dg::AbstractVector, v::DiffGraph) = v == dg.primal_to_diff
Base.eltype(::DiffGraph) = Union{Int, Nothing}
Base.size(dg::DiffGraph) = size(dg.primal_to_diff)
Base.length(dg::DiffGraph) = length(dg.primal_to_diff)
Base.getindex(dg::DiffGraph, var::Integer) = dg.primal_to_diff[var]
Base.getindex(dg::DiffGraph, a::AbstractArray) = [dg[x] for x in a]

function Base.setindex!(dg::DiffGraph, val::Union{Integer, Nothing}, var::Integer)
    if dg.diff_to_primal !== nothing
        old_pd = dg.primal_to_diff[var]
        if old_pd !== nothing
            dg.diff_to_primal[old_pd] = nothing
        end
        if val !== nothing
            #old_dp = dg.diff_to_primal[val]
            #old_dp === nothing || error("Variable already assigned.")
            dg.diff_to_primal[val] = var
        end
    end
    return dg.primal_to_diff[var] = val
end
Base.iterate(dg::DiffGraph, state...) = iterate(dg.primal_to_diff, state...)

function complete(dg::DiffGraph)
    dg.diff_to_primal !== nothing && return dg
    diff_to_primal = Union{Int, Nothing}[nothing for _ in 1:length(dg.primal_to_diff)]
    for (var, diff) in edges(dg)
        diff_to_primal[diff] = var
    end
    return DiffGraph(dg.primal_to_diff, diff_to_primal)
end

function invview(dg::DiffGraph)
    require_complete(dg)
    return DiffGraph(dg.diff_to_primal, dg.primal_to_diff)
end

struct DiffChainIterator{Descend}
    var_to_diff::DiffGraph
    v::Int
end

function Base.iterate(di::DiffChainIterator{Descend}, v = nothing) where {Descend}
    if v === nothing
        vv = di.v
        return (vv, vv)
    end
    g = Descend ? invview(di.var_to_diff) : di.var_to_diff
    v′ = g[v]
    v′ === nothing ? nothing : (v′, v′)
end

abstract type TransformationState{T} end
abstract type AbstractTearingState{T} <: TransformationState{T} end

get_fullvars(ts::TransformationState) = ts.fullvars
has_equations(::TransformationState) = true

Base.@kwdef mutable struct SystemStructure
    """Maps the index of variable x to the index of variable D(x)."""
    var_to_diff::DiffGraph
    """Maps the index of an algebraic equation to the index of the equation it is differentiated into."""
    eq_to_diff::DiffGraph
    # Can be access as
    # `graph` to automatically look at the bipartite graph
    # or as `torn` to assert that tearing has run.
    """Graph that connects equations to variables that appear in them."""
    graph::BipartiteGraph{Int, Nothing}
    """Graph that connects equations to the variable they will be solved for during simplification."""
    solvable_graph::Union{BipartiteGraph{Int, Nothing}, Nothing}
    """Variable types (brownian, variable, parameter) in the system."""
    var_types::Union{Vector{VariableType}, Nothing}
    """Whether the system is discrete."""
    only_discrete::Bool
end

function Base.copy(structure::SystemStructure)
    var_types = structure.var_types === nothing ? nothing : copy(structure.var_types)
    SystemStructure(copy(structure.var_to_diff), copy(structure.eq_to_diff),
        copy(structure.graph), copy(structure.solvable_graph),
        var_types, structure.only_discrete)
end

is_only_discrete(s::SystemStructure) = s.only_discrete
isdervar(s::SystemStructure, i) = invview(s.var_to_diff)[i] !== nothing
function isalgvar(s::SystemStructure, i)
    s.var_to_diff[i] === nothing &&
        invview(s.var_to_diff)[i] === nothing
end
function isdiffvar(s::SystemStructure, i)
    s.var_to_diff[i] !== nothing && invview(s.var_to_diff)[i] === nothing
end

function dervars_range(s::SystemStructure)
    Iterators.filter(Base.Fix1(isdervar, s), Base.OneTo(ndsts(s.graph)))
end
function diffvars_range(s::SystemStructure)
    Iterators.filter(Base.Fix1(isdiffvar, s), Base.OneTo(ndsts(s.graph)))
end
function algvars_range(s::SystemStructure)
    Iterators.filter(Base.Fix1(isalgvar, s), Base.OneTo(ndsts(s.graph)))
end

function algeqs(s::SystemStructure)
    BitSet(findall(map(1:nsrcs(s.graph)) do eq
        all(v -> !isdervar(s, v), 𝑠neighbors(s.graph, eq))
    end))
end

function complete!(s::SystemStructure)
    s.var_to_diff = complete(s.var_to_diff)
    s.eq_to_diff = complete(s.eq_to_diff)
    s.graph = complete(s.graph)
    if s.solvable_graph !== nothing
        s.solvable_graph = complete(s.solvable_graph)
    end
    s
end

mutable struct TearingState{T <: AbstractSystem} <: AbstractTearingState{T}
    """The system of equations."""
    sys::T
    """The set of variables of the system."""
    fullvars::Vector{BasicSymbolic}
    structure::SystemStructure
    extra_eqs::Vector
    param_derivative_map::Dict{BasicSymbolic, Any}
    original_eqs::Vector{Equation}
    """
    Additional user-provided observed equations. The variables calculated here
    are not used in the rest of the system.
    """
    additional_observed::Vector{Equation}
    statemachines::Vector{T}
end

TransformationState(sys::AbstractSystem) = TearingState(sys)
function system_subset(ts::TearingState, ieqs::Vector{Int}, iieqs::Vector{Int}, ivars::Vector{Int})
    eqs = equations(ts)
    initeqs = initialization_equations(ts.sys)
    @set! ts.sys.eqs = eqs[ieqs]
    @set! ts.sys.initialization_eqs = initeqs[iieqs]
    @set! ts.original_eqs = ts.original_eqs[ieqs]
    @set! ts.structure = system_subset(ts.structure, ieqs, ivars)
    if all(eq -> eq.rhs isa StateMachineOperator, get_eqs(ts.sys))
        names = Symbol[]
        for eq in get_eqs(ts.sys)
            if eq.lhs isa Transition
                push!(names, first(namespace_hierarchy(nameof(eq.rhs.from))))
                push!(names, first(namespace_hierarchy(nameof(eq.rhs.to))))
            elseif eq.lhs isa InitialState
                push!(names, first(namespace_hierarchy(nameof(eq.rhs.s))))
            else
                error("Unhandled state machine operator")
            end
        end
        @set! ts.statemachines = filter(x -> nameof(x) in names, ts.statemachines)
    else
        @set! ts.statemachines = eltype(ts.statemachines)[]
    end
    @set! ts.fullvars = ts.fullvars[ivars]
    ts
end

function system_subset(structure::SystemStructure, ieqs::Vector{Int}, ivars::Vector{Int})
    @unpack graph = structure
    fadj = Vector{Int}[]
    eq_to_diff = DiffGraph(length(ieqs))
    var_to_diff = DiffGraph(length(ivars))

    ne = 0
    old_to_new_var = zeros(Int, ndsts(graph))
    for (i, iv) in enumerate(ivars)
        old_to_new_var[iv] = i
    end
    for (i, iv) in enumerate(ivars)
        structure.var_to_diff[iv] === nothing && continue
        var_to_diff[i] = old_to_new_var[structure.var_to_diff[iv]]
    end
    for (j, eq_i) in enumerate(ieqs)
        var_adj = [old_to_new_var[i] for i in graph.fadjlist[eq_i]]
        @assert all(!iszero, var_adj)
        ne += length(var_adj)
        push!(fadj, var_adj)
        eq_to_diff[j] = structure.eq_to_diff[eq_i]
    end
    @set! structure.graph = complete(BipartiteGraph(ne, fadj, length(ivars)))
    @set! structure.eq_to_diff = eq_to_diff
    @set! structure.var_to_diff = complete(var_to_diff)
    structure
end

function Base.show(io::IO, state::TearingState)
    print(io, "TearingState of ", typeof(state.sys))
end

struct EquationsView{T} <: AbstractVector{Any}
    ts::TearingState{T}
end
equations(ts::TearingState) = EquationsView(ts)
Base.size(ev::EquationsView) = (length(equations(ev.ts.sys)) + length(ev.ts.extra_eqs),)
function Base.getindex(ev::EquationsView, i::Integer)
    eqs = equations(ev.ts.sys)
    if i > length(eqs)
        return ev.ts.extra_eqs[i - length(eqs)]
    end
    return eqs[i]
end
function Base.push!(ev::EquationsView, eq)
    push!(ev.ts.extra_eqs, eq)
end

function is_time_dependent_parameter(p, allps, iv)
    return iv !== nothing && p in allps && iscall(p) &&
           (operation(p) === getindex &&
            is_time_dependent_parameter(arguments(p)[1], allps, iv) ||
            (args = arguments(p); length(args)) == 1 && isequal(only(args), iv))
end

function symbolic_contains(var, set)
    var in set ||
        symbolic_type(var) == ArraySymbolic() &&
        Symbolics.shape(var) != Symbolics.Unknown() &&
        all(x -> x in set, Symbolics.scalarize(var))
end

"""
    $(TYPEDSIGNATURES)

Descend through the system hierarchy and look for statemachines. Remove equations from
the inner statemachine systems. Return the new `sys` and an array of top-level
statemachines.
"""
function extract_top_level_statemachines(sys::AbstractSystem)
    eqs = get_eqs(sys)

    if !isempty(eqs) && all(eq -> eq.lhs isa StateMachineOperator, eqs)
        # top-level statemachine
        with_removed = @set sys.systems = map(remove_child_equations, get_systems(sys))
        return with_removed, [sys]
    elseif !isempty(eqs) && any(eq -> eq.lhs isa StateMachineOperator, eqs)
        # error: can't mix
        error("Mixing statemachine equations and standard equations in a top-level statemachine is not allowed.")
    else
        # descend
        subsystems = get_systems(sys)
        newsubsystems = eltype(subsystems)[]
        statemachines = eltype(subsystems)[]
        for subsys in subsystems
            newsubsys, sub_statemachines = extract_top_level_statemachines(subsys)
            push!(newsubsystems, newsubsys)
            append!(statemachines, sub_statemachines)
        end
        @set! sys.systems = newsubsystems
        return sys, statemachines
    end
end

"""
    $(TYPEDSIGNATURES)

Return `sys` with all equations (including those in subsystems) removed.
"""
function remove_child_equations(sys::AbstractSystem)
    @set! sys.eqs = eltype(get_eqs(sys))[]
    @set! sys.systems = map(remove_child_equations, get_systems(sys))
    return sys
end

function TearingState(sys; quick_cancel = false, check = true, sort_eqs = true)
    # flatten system
    sys = flatten(sys)
    sys = process_parameter_equations(sys)
    ivs = independent_variables(sys)
    iv = length(ivs) == 1 ? ivs[1] : nothing
    # flatten array equations
    eqs = flatten_equations(equations(sys))
    original_eqs = copy(eqs)
    neqs = length(eqs)
    param_derivative_map = Dict{BasicSymbolic, Any}()
    # * Scalarize unknowns
    dvs = Set{BasicSymbolic}()
    fullvars = BasicSymbolic[]
    for x in unknowns(sys)
        push!(dvs, x)
        xx = Symbolics.scalarize(x)
        if xx isa AbstractArray
            union!(dvs, xx)
        end
    end
    ps = Set{Symbolic}()
    for x in full_parameters(sys)
        push!(ps, x)
        if symbolic_type(x) == ArraySymbolic() && Symbolics.shape(x) != Symbolics.Unknown()
            xx = Symbolics.scalarize(x)
            union!(ps, xx)
        end
    end
    browns = Set{BasicSymbolic}()
    for x in brownians(sys)
        push!(browns, x)
        xx = Symbolics.scalarize(x)
        if xx isa AbstractArray
            union!(browns, xx)
        end
    end
    var2idx = Dict{BasicSymbolic, Int}()
    var_types = VariableType[]
    addvar! = let fullvars = fullvars, dvs = dvs, var2idx = var2idx, var_types = var_types
        (var, vtype) -> get!(var2idx, var) do
            push!(dvs, var)
            push!(fullvars, var)
            push!(var_types, vtype)
            return length(fullvars)
        end
    end

    # build symbolic incidence
    symbolic_incidence = Vector{BasicSymbolic}[]
    varsbuf = Set()
    eqs_to_retain = trues(length(eqs))
    for (i, eq) in enumerate(eqs)
        _eq = eq
        if iscall(eq.lhs) && (op = operation(eq.lhs)) isa Differential &&
           isequal(op.x, iv) && is_time_dependent_parameter(only(arguments(eq.lhs)), ps, iv)
            # parameter derivatives are opted out by specifying `D(p) ~ missing`, but
            # we want to store `nothing` in the map because that means `fast_substitute`
            # will ignore the rule. We will this identify the presence of `eq′.lhs` in
            # the differentiated expression and error.
            param_derivative_map[eq.lhs] = coalesce(eq.rhs, nothing)
            eqs_to_retain[i] = false
            # change the equation if the RHS is `missing` so the rest of this loop works
            eq = 0.0 ~ coalesce(eq.rhs, 0.0)
        end
        is_statemachine_equation = false
        if eq.lhs isa StateMachineOperator
            is_statemachine_equation = true
            eq = eq
            rhs = eq.rhs
        elseif _iszero(eq.lhs)
            rhs = quick_cancel ? quick_cancel_expr(eq.rhs) : eq.rhs
        else
            lhs = quick_cancel ? quick_cancel_expr(eq.lhs) : eq.lhs
            rhs = quick_cancel ? quick_cancel_expr(eq.rhs) : eq.rhs
            eq = 0 ~ rhs - lhs
        end
        empty!(varsbuf)
        vars!(varsbuf, eq; op = Symbolics.Operator)
        incidence = Set{BasicSymbolic}()
        isalgeq = true
        for v in varsbuf
            # additionally track brownians in fullvars
            if v in browns
                addvar!(v, BROWNIAN)
                push!(incidence, v)
            end

            # TODO: Can we handle this without `isparameter`?
            if symbolic_contains(v, ps) ||
               getmetadata(v, SymScope, LocalScope()) isa GlobalScope && isparameter(v)
                if is_time_dependent_parameter(v, ps, iv) &&
                   !haskey(param_derivative_map, Differential(iv)(v))
                    # Parameter derivatives default to zero - they stay constant
                    # between callbacks
                    param_derivative_map[Differential(iv)(v)] = 0.0
                end
                continue
            end

            isequal(v, iv) && continue
            isdelay(v, iv) && continue

            if !symbolic_contains(v, dvs)
                isvalid = iscall(v) &&
                          (operation(v) isa Shift || isempty(arguments(v)) ||
                           is_transparent_operator(operation(v)))
                v′ = v
                while !isvalid && iscall(v′) && operation(v′) isa Union{Differential, Shift}
                    v′ = arguments(v′)[1]
                    if v′ in dvs || getmetadata(v′, SymScope, LocalScope()) isa GlobalScope
                        isvalid = true
                        break
                    end
                end
                if !isvalid
                    throw(ArgumentError("$v is present in the system but $v′ is not an unknown."))
                end

                addvar!(v, VARIABLE)
                if iscall(v) && operation(v) isa Symbolics.Operator && !isdifferential(v) &&
                   (it = input_timedomain(v)) !== nothing
                    for v′ in arguments(v)
                        addvar!(setmetadata(v′, VariableTimeDomain, it), VARIABLE)
                    end
                end
            end

            isalgeq &= !isdifferential(v)

            if symbolic_type(v) == ArraySymbolic()
                vv = collect(v)
                union!(incidence, vv)
                map(vv) do vi
                    addvar!(vi, VARIABLE)
                end
            else
                push!(incidence, v)
                addvar!(v, VARIABLE)
            end
        end
        if isalgeq || is_statemachine_equation
            eqs[i] = eq
        else
            eqs[i] = eqs[i].lhs ~ rhs
        end
        push!(symbolic_incidence, collect(incidence))
    end

    dervaridxs = OrderedSet{Int}()
    for (i, v) in enumerate(fullvars)
        while isdifferential(v)
            push!(dervaridxs, i)
            v = arguments(v)[1]
            i = addvar!(v, VARIABLE)
        end
    end
    eqs = eqs[eqs_to_retain]
    original_eqs = original_eqs[eqs_to_retain]
    neqs = length(eqs)
    symbolic_incidence = symbolic_incidence[eqs_to_retain]

    if sort_eqs
        # sort equations lexicographically to reduce simplification issues
        # depending on order due to NP-completeness of tearing.
        sortidxs = Base.sortperm(string.(eqs)) # "by = string" creates more strings
        eqs = eqs[sortidxs]
        original_eqs = original_eqs[sortidxs]
        symbolic_incidence = symbolic_incidence[sortidxs]
    end

    # Handle shifts - find lowest shift and add intermediates with derivative edges
    ### Handle discrete variables
    lowest_shift = Dict()
    for var in fullvars
        if ModelingToolkit.isoperator(var, ModelingToolkit.Shift)
            steps = operation(var).steps
            if steps > 0
                error("Only non-positive shifts allowed. Found $var with a shift of $steps")
            end
            v = arguments(var)[1]
            lowest_shift[v] = min(get(lowest_shift, v, 0), steps)
        end
    end
    for var in fullvars
        if ModelingToolkit.isoperator(var, ModelingToolkit.Shift)
            op = operation(var)
            steps = op.steps
            v = arguments(var)[1]
            lshift = lowest_shift[v]
            tt = op.t
        elseif haskey(lowest_shift, var)
            lshift = lowest_shift[var]
            steps = 0
            tt = iv
            v = var
        else
            continue
        end
        if lshift < steps
            push!(dervaridxs, var2idx[var])
        end
        for s in (steps - 1):-1:(lshift + 1)
            sf = Shift(tt, s)
            dvar = sf(v)
            idx = addvar!(dvar, VARIABLE)
            if !(idx in dervaridxs)
                push!(dervaridxs, idx)
            end
        end
    end

    # sort `fullvars` such that the mass matrix is as diagonal as possible.
    dervaridxs = collect(dervaridxs)
    sorted_fullvars = OrderedSet(fullvars[dervaridxs])
    var_to_old_var = Dict(zip(fullvars, fullvars))
    for dervaridx in dervaridxs
        dervar = fullvars[dervaridx]
        diffvar = var_to_old_var[lower_order_var(dervar, iv)]
        if !(diffvar in sorted_fullvars)
            push!(sorted_fullvars, diffvar)
        end
    end
    for v in fullvars
        if !(v in sorted_fullvars)
            push!(sorted_fullvars, v)
        end
    end
    new_fullvars = collect(sorted_fullvars)
    sortperm = indexin(new_fullvars, fullvars)
    fullvars = new_fullvars
    var_types = var_types[sortperm]
    var2idx = Dict(fullvars .=> eachindex(fullvars))
    dervaridxs = 1:length(dervaridxs)

    # build `var_to_diff`
    nvars = length(fullvars)
    diffvars = []
    var_to_diff = DiffGraph(nvars, true)
    for dervaridx in dervaridxs
        dervar = fullvars[dervaridx]
        diffvar = lower_order_var(dervar, iv)
        diffvaridx = var2idx[diffvar]
        push!(diffvars, diffvar)
        var_to_diff[diffvaridx] = dervaridx
    end

    # build incidence graph
    graph = BipartiteGraph(neqs, nvars, Val(false))
    for (ie, vars) in enumerate(symbolic_incidence), v in vars

        jv = var2idx[v]
        add_edge!(graph, ie, jv)
    end

    @set! sys.eqs = eqs

    eq_to_diff = DiffGraph(nsrcs(graph))

    ts = TearingState(sys, fullvars,
        SystemStructure(complete(var_to_diff), complete(eq_to_diff),
            complete(graph), nothing, var_types, false),
        Any[], param_derivative_map, original_eqs, Equation[], typeof(sys)[])
    return ts
end

"""
    $(TYPEDSIGNATURES)

Preemptively identify observed equations in the system and tear them. This happens before
any simplification. The equations torn by this process are ones that are already given in
an explicit form in the system and where the LHS is not present in any other equation of
the system except for other such preempitvely torn equations.
"""
function trivial_tearing!(ts::TearingState)
    @assert length(ts.original_eqs) == length(equations(ts))
    # equations that can be trivially torn an observed equations
    trivial_idxs = BitSet()
    # equations to never check
    blacklist = BitSet()
    torn_eqs = Equation[]
    # variables that have been matched to trivially torn equations
    matched_vars = BitSet()
    # variable to index in fullvars
    var_to_idx = Dict{Any, Int}(ts.fullvars .=> eachindex(ts.fullvars))
    sys_eqs = equations(ts)

    complete!(ts.structure)
    var_to_diff = ts.structure.var_to_diff
    graph = ts.structure.graph
    while true
        # track whether we added an equation to the trivial list this iteration
        added_equation = false
        for (i, eq) in enumerate(ts.original_eqs)
            # don't check already torn equations
            i in trivial_idxs && continue
            i in blacklist && continue
            # ensure it is an observed equation matched to a variable in fullvars
            vari = get(var_to_idx, eq.lhs, 0)
            iszero(vari) && continue
            # don't tear irreducible variables
            if isirreducible(eq.lhs)
                push!(blacklist, i)
                continue
            end
            # Edge case for `var ~ var` equations. They don't show up in the incidence
            # graph because `TearingState` makes them `0 ~ 0`, but they do cause `var`
            # to show up twice in `original_eqs` which fails the assertion.
            sys_eq = sys_eqs[i]
            if isequal(sys_eq.lhs, 0) && isequal(sys_eq.rhs, 0)
                continue
            end

            # if a variable was the LHS of two trivial observed equations, we wouldn't have
            # included it in the list. Error if somehow it made it through.
            @assert !(vari in matched_vars)
            # don't tear differential/shift equations (or differentiated/shifted variables)
            var_to_diff[vari] === nothing || continue
            invview(var_to_diff)[vari] === nothing || continue
            # get the equations that the candidate matched variable is present in, except
            # those equations which have already been torn as observed
            eqidxs = setdiff(𝑑neighbors(graph, vari), trivial_idxs)
            # it should only be present in this equation
            length(eqidxs) == 1 || continue
            eqi = only(eqidxs)
            @assert eqi == i

            # for every variable present in this equation, make sure it isn't _only_
            # present in trivial equations
            isvalid = true
            for v in 𝑠neighbors(graph, eqi)
                v == vari && continue
                v in matched_vars && continue
                # `> 1` and not `0` because one entry will be this equation (`eqi`)
                isvalid &= count(!in(trivial_idxs), 𝑑neighbors(graph, v)) > 1
                isvalid || break
            end
            isvalid || continue
            # skip if the LHS is present in the RHS, since then this isn't explicit
            if occursin(eq.lhs, eq.rhs)
                push!(blacklist, i)
                continue
            end

            added_equation = true
            push!(trivial_idxs, eqi)
            push!(torn_eqs, eq)
            push!(matched_vars, vari)
        end

        # if we didn't add an equation this iteration, we won't add one next iteration
        added_equation || break
    end

    deleteat!(var_to_diff.primal_to_diff, matched_vars)
    deleteat!(var_to_diff.diff_to_primal, matched_vars)
    deleteat!(ts.structure.eq_to_diff.primal_to_diff, trivial_idxs)
    deleteat!(ts.structure.eq_to_diff.diff_to_primal, trivial_idxs)
    delete_srcs!(ts.structure.graph, trivial_idxs; rm_verts = true)
    delete_dsts!(ts.structure.graph, matched_vars; rm_verts = true)
    if ts.structure.solvable_graph !== nothing
        delete_srcs!(ts.structure.solvable_graph, trivial_idxs; rm_verts = true)
        delete_dsts!(ts.structure.solvable_graph, matched_vars; rm_verts = true)
    end
    if ts.structure.var_types !== nothing
        deleteat!(ts.structure.var_types, matched_vars)
    end
    deleteat!(ts.fullvars, matched_vars)
    deleteat!(ts.original_eqs, trivial_idxs)
    ts.additional_observed = torn_eqs
    sys = ts.sys
    eqs = copy(get_eqs(sys))
    deleteat!(eqs, trivial_idxs)
    @set! sys.eqs = eqs
    ts.sys = sys
    return ts
end

function lower_order_var(dervar, t)
    if isdifferential(dervar)
        diffvar = arguments(dervar)[1]
    elseif ModelingToolkit.isoperator(dervar, ModelingToolkit.Shift)
        s = operation(dervar)
        step = s.steps - 1
        vv = arguments(dervar)[1]
        if step != 0
            diffvar = Shift(s.t, step)(vv)
        else
            diffvar = vv
        end
    else
        return Shift(t, -1)(dervar)
    end
    diffvar
end

function shift_discrete_system(ts::TearingState)
    @unpack fullvars, sys = ts
    discvars = OrderedSet()
    eqs = equations(sys)
    for eq in eqs
        vars!(discvars, eq; op = Union{Sample, Hold, Pre})
    end
    iv = get_iv(sys)

    discmap = Dict(k => StructuralTransformations.simplify_shifts(Shift(iv, 1)(k))
    for k in discvars
    if any(isequal(k), fullvars) && !isa(operation(k), Union{Sample, Hold, Pre}))

    for i in eachindex(fullvars)
        fullvars[i] = StructuralTransformations.simplify_shifts(fast_substitute(
            fullvars[i], discmap; operator = Union{Sample, Hold, Pre}))
    end
    for i in eachindex(eqs)
        eqs[i] = StructuralTransformations.simplify_shifts(fast_substitute(
            eqs[i], discmap; operator = Union{Sample, Hold, Pre}))
    end
    @set! ts.sys.eqs = eqs
    @set! ts.fullvars = fullvars
    return ts
end

using .BipartiteGraphs: Label, BipartiteAdjacencyList
struct SystemStructurePrintMatrix <:
       AbstractMatrix{Union{Label, BipartiteAdjacencyList}}
    bpg::BipartiteGraph
    highlight_graph::Union{Nothing, BipartiteGraph}
    var_to_diff::DiffGraph
    eq_to_diff::DiffGraph
    var_eq_matching::Union{Matching, Nothing}
end

"""
Create a SystemStructurePrintMatrix to display the contents
of the provided SystemStructure.
"""
function SystemStructurePrintMatrix(s::SystemStructure)
    return SystemStructurePrintMatrix(complete(s.graph),
        s.solvable_graph === nothing ? nothing :
        complete(s.solvable_graph),
        complete(s.var_to_diff),
        complete(s.eq_to_diff),
        nothing)
end
Base.size(bgpm::SystemStructurePrintMatrix) = (max(nsrcs(bgpm.bpg), ndsts(bgpm.bpg)) + 1, 9)
function compute_diff_label(diff_graph, i, symbol)
    di = i - 1 <= length(diff_graph) ? diff_graph[i - 1] : nothing
    return di === nothing ? Label("") : Label(string(di, symbol))
end
function Base.getindex(bgpm::SystemStructurePrintMatrix, i::Integer, j::Integer)
    checkbounds(bgpm, i, j)
    if i <= 1
        return (Label.(("#", "∂ₜ", " ", "  eq", "", "#", "∂ₜ", " ", "  v")))[j]
    elseif j == 5
        colors = Base.text_colors
        return Label("|", :light_black)
    elseif j == 2
        return compute_diff_label(bgpm.eq_to_diff, i, '↑')
    elseif j == 3
        return compute_diff_label(invview(bgpm.eq_to_diff), i, '↓')
    elseif j == 7
        return compute_diff_label(bgpm.var_to_diff, i, '↑')
    elseif j == 8
        return compute_diff_label(invview(bgpm.var_to_diff), i, '↓')
    elseif j == 1
        return Label((i - 1 <= length(bgpm.eq_to_diff)) ? string(i - 1) : "")
    elseif j == 6
        return Label((i - 1 <= length(bgpm.var_to_diff)) ? string(i - 1) : "")
    elseif j == 4
        return BipartiteAdjacencyList(
            i - 1 <= nsrcs(bgpm.bpg) ?
            𝑠neighbors(bgpm.bpg, i - 1) : nothing,
            bgpm.highlight_graph !== nothing &&
            i - 1 <= nsrcs(bgpm.highlight_graph) ?
            Set(𝑠neighbors(bgpm.highlight_graph, i - 1)) :
            nothing,
            bgpm.var_eq_matching !== nothing &&
            (i - 1 <= length(invview(bgpm.var_eq_matching))) ?
            invview(bgpm.var_eq_matching)[i - 1] : unassigned)
    elseif j == 9
        match = unassigned
        if bgpm.var_eq_matching !== nothing && i - 1 <= length(bgpm.var_eq_matching)
            match = bgpm.var_eq_matching[i - 1]
            isa(match, Union{Int, Unassigned}) || (match = true) # Selected Unknown
        end
        return BipartiteAdjacencyList(
            i - 1 <= ndsts(bgpm.bpg) ?
            𝑑neighbors(bgpm.bpg, i - 1) : nothing,
            bgpm.highlight_graph !== nothing &&
            i - 1 <= ndsts(bgpm.highlight_graph) ?
            Set(𝑑neighbors(bgpm.highlight_graph, i - 1)) :
            nothing, match)
    else
        @assert false
    end
end

function Base.show(io::IO, mime::MIME"text/plain", s::SystemStructure)
    @unpack graph, solvable_graph, var_to_diff, eq_to_diff = s
    if !get(io, :limit, true) || !get(io, :mtk_limit, true)
        print(io, "SystemStructure with ", length(s.graph.fadjlist), " equations and ",
            isa(s.graph.badjlist, Int) ? s.graph.badjlist : length(s.graph.badjlist),
            " variables\n")
        Base.print_matrix(io, SystemStructurePrintMatrix(s))
    else
        S = incidence_matrix(s.graph, Num(Sym{Real}(:×)))
        print(io, "Incidence matrix:")
        show(io, mime, S)
    end
end

struct MatchedSystemStructure
    structure::SystemStructure
    var_eq_matching::Matching
end

"""
Create a SystemStructurePrintMatrix to display the contents
of the provided MatchedSystemStructure.
"""
function SystemStructurePrintMatrix(ms::MatchedSystemStructure)
    return SystemStructurePrintMatrix(complete(ms.structure.graph),
        complete(ms.structure.solvable_graph),
        complete(ms.structure.var_to_diff),
        complete(ms.structure.eq_to_diff),
        complete(ms.var_eq_matching,
            nsrcs(ms.structure.graph)))
end

function Base.copy(ms::MatchedSystemStructure)
    MatchedSystemStructure(Base.copy(ms.structure), Base.copy(ms.var_eq_matching))
end

function Base.show(io::IO, mime::MIME"text/plain", ms::MatchedSystemStructure)
    s = ms.structure
    @unpack graph, solvable_graph, var_to_diff, eq_to_diff = s
    print(io, "Matched SystemStructure with ", length(graph.fadjlist), " equations and ",
        isa(graph.badjlist, Int) ? graph.badjlist : length(graph.badjlist),
        " variables\n")
    Base.print_matrix(io, SystemStructurePrintMatrix(ms))
    printstyled(io, "\n\nLegend: ")
    printstyled(io, "Solvable")
    print(io, " | ")
    printstyled(io, "(Solvable + Matched)", color = :light_yellow)
    print(io, " | ")
    printstyled(io, "Unsolvable", color = :light_black)
    print(io, " | ")
    printstyled(io, "(Unsolvable + Matched)", color = :magenta)
    print(io, " | ")
    printstyled(io, " ∫", color = :cyan)
    printstyled(io, " SelectedState")
end

function make_eqs_zero_equals!(ts::TearingState)
    neweqs = map(enumerate(get_eqs(ts.sys))) do kvp
        i, eq = kvp
        isalgeq = true
        for j in 𝑠neighbors(ts.structure.graph, i)
            isalgeq &= invview(ts.structure.var_to_diff)[j] === nothing
        end
        if isalgeq
            return 0 ~ eq.rhs - eq.lhs
        else
            return eq
        end
    end
    copyto!(get_eqs(ts.sys), neweqs)
end

function mtkcompile!(state::TearingState;
        check_consistency = true, fully_determined = true, warn_initialize_determined = true,
        inputs = Any[], outputs = Any[],
        disturbance_inputs = Any[],
        kwargs...)
    if !is_time_dependent(state.sys)
        return _mtkcompile!(state; check_consistency,
            inputs, outputs, disturbance_inputs,
            fully_determined, kwargs...)
    end
    # split_system returns one or two systems and the inputs for each
    # mod clock inference to be binary
    # if it's continuous keep going, if not then error unless given trait impl in additional passes
    ci = ModelingToolkit.ClockInference(state)
    ci = ModelingToolkit.infer_clocks!(ci)
    time_domains = merge(Dict(state.fullvars .=> ci.var_domain),
        Dict(default_toterm.(state.fullvars) .=> ci.var_domain))
    tss, clocked_inputs, continuous_id, id_to_clock = ModelingToolkit.split_system(ci)
    if !isempty(tss) && continuous_id == 0
        # do a trait check here - handle fully discrete system
        additional_passes = get(kwargs, :additional_passes, nothing)
        if !isnothing(additional_passes) && any(discrete_compile_pass, additional_passes)
            # take the first discrete compilation pass given for now
            discrete_pass_idx = findfirst(discrete_compile_pass, additional_passes)
            discrete_compile = additional_passes[discrete_pass_idx]
            deleteat!(additional_passes, discrete_pass_idx)
            return discrete_compile(tss, clocked_inputs, ci)
        end
        throw(HybridSystemNotSupportedException("""
        Discrete systems with multiple clocks are not supported with the standard \
        MTK compiler.
        """))
    end
    if length(tss) > 1
        make_eqs_zero_equals!(tss[continuous_id])
        # simplify as normal
        sys = _mtkcompile!(tss[continuous_id];
            inputs = [inputs; clocked_inputs[continuous_id]], outputs, disturbance_inputs,
            check_consistency, fully_determined,
            kwargs...)
        additional_passes = get(kwargs, :additional_passes, nothing)
        if !isnothing(additional_passes) && any(discrete_compile_pass, additional_passes)
            discrete_pass_idx = findfirst(discrete_compile_pass, additional_passes)
            discrete_compile = additional_passes[discrete_pass_idx]
            deleteat!(additional_passes, discrete_pass_idx)
            # in the case of a hybrid system, the discrete_compile pass should take the currents of sys.discrete_subsystems
            # and modifies discrete_subsystems to bea tuple of the io and anything else, while adding or manipulating the rest of sys as needed
            return discrete_compile(
                sys, tss[[i for i in eachindex(tss) if i != continuous_id]],
                clocked_inputs, ci, id_to_clock)
        end
        throw(HybridSystemNotSupportedException("""
        Hybrid continuous-discrete systems are currently not supported with \
        the standard MTK compiler. This system requires JuliaSimCompiler.jl, \
        see https://help.juliahub.com/juliasimcompiler/stable/
        """))
    end
    if get_is_discrete(state.sys) ||
       continuous_id == 1 && any(Base.Fix2(isoperator, Shift), state.fullvars)
        state.structure.only_discrete = true
        state = shift_discrete_system(state)
        sys = state.sys
        @set! sys.is_discrete = true
        state.sys = sys
    end

    sys = _mtkcompile!(state; check_consistency,
        inputs, outputs, disturbance_inputs,
        fully_determined, kwargs...)
    return sys
end

function _mtkcompile!(state::TearingState;
        check_consistency = true, fully_determined = true, warn_initialize_determined = false,
        dummy_derivative = true,
        inputs = Any[], outputs = Any[],
        disturbance_inputs = Any[],
        kwargs...)
    if fully_determined isa Bool
        check_consistency &= fully_determined
    else
        check_consistency = true
    end
    orig_inputs = Set()
    ModelingToolkit.markio!(state, orig_inputs, inputs, outputs, disturbance_inputs)
    state = ModelingToolkit.inputs_to_parameters!(state, [inputs; disturbance_inputs])
    trivial_tearing!(state)
    sys, mm = ModelingToolkit.alias_elimination!(state; fully_determined, kwargs...)
    if check_consistency
        fully_determined = ModelingToolkit.check_consistency(
            state, orig_inputs; nothrow = fully_determined === nothing)
    end
    if fully_determined && dummy_derivative
        sys = ModelingToolkit.dummy_derivative(
            sys, state; mm, check_consistency, kwargs...)
    elseif fully_determined
        var_eq_matching = pantelides!(state; finalize = false, kwargs...)
        sys = pantelides_reassemble(state, var_eq_matching)
        state = TearingState(sys)
        sys, mm = ModelingToolkit.alias_elimination!(state; fully_determined, kwargs...)
        sys = ModelingToolkit.dummy_derivative(
            sys, state; mm, check_consistency, fully_determined, kwargs...)
    else
        sys = ModelingToolkit.tearing(
            sys, state; mm, check_consistency, fully_determined, kwargs...)
    end
    fullunknowns = [observables(sys); unknowns(sys)]
    @set! sys.observed = ModelingToolkit.topsort_equations(observed(sys), fullunknowns)

    ModelingToolkit.invalidate_cache!(sys)
end

struct DifferentiatedVariableNotUnknownError <: Exception
    differentiated::Any
    undifferentiated::Any
end

function Base.showerror(io::IO, err::DifferentiatedVariableNotUnknownError)
    undiff = err.undifferentiated
    diff = err.differentiated
    print(io,
        "Variable $undiff occurs differentiated as $diff but is not an unknown of the system.")
    scope = getmetadata(undiff, SymScope, LocalScope())
    depth = expected_scope_depth(scope)
    if depth > 0
        print(io,
            "\nVariable $undiff expects $depth more levels in the hierarchy to be an unknown.")
    end
end
