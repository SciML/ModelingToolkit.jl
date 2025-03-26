using DataStructures
using Symbolics: linear_expansion, unwrap, Connection, Transition, InitialState
using SymbolicUtils: iscall, operation, arguments, Symbolic
using SymbolicUtils: quick_cancel, maketerm
using ..ModelingToolkit
import ..ModelingToolkit: isdiffeq, var_from_nested_derivative, vars!, flatten,
                          value, InvalidSystemException, isdifferential, _iszero,
                          isparameter, isconstant,
                          independent_variables, SparseMatrixCLIL, AbstractSystem,
                          equations, isirreducible, input_timedomain, TimeDomain,
                          InferredTimeDomain,
                          VariableType, getvariabletype, has_equations, ODESystem
using ..BipartiteGraphs
import ..BipartiteGraphs: invview, complete
using Graphs
using UnPack
using Setfield
using SparseArrays

function quick_cancel_expr(expr)
    Rewriters.Postwalk(quick_cancel,
        similarterm = (x, f, args; kws...) -> maketerm(typeof(x), f, args,
            SymbolicUtils.metadata(x),
            kws...))(expr)
end

export SystemStructure, TransformationState, TearingState, structural_simplify!
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
    # Maps the (index of) a variable to the (index of) the variable describing
    # its derivative.
    var_to_diff::DiffGraph
    eq_to_diff::DiffGraph
    # Can be access as
    # `graph` to automatically look at the bipartite graph
    # or as `torn` to assert that tearing has run.
    graph::BipartiteGraph{Int, Nothing}
    solvable_graph::Union{BipartiteGraph{Int, Nothing}, Nothing}
    var_types::Union{Vector{VariableType}, Nothing}
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
    sys::T
    original_eqs::Vector{Equation}
    fullvars::Vector
    structure::SystemStructure
    extra_eqs::Vector
    statemachines::Vector{T}
end

TransformationState(sys::AbstractSystem) = TearingState(sys)
function system_subset(ts::TearingState, ieqs::Vector{Int})
    eqs = equations(ts)
    @set! ts.original_eqs = ts.original_eqs[ieqs]
    @set! ts.sys.eqs = eqs[ieqs]
    @set! ts.structure = system_subset(ts.structure, ieqs)
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
    ts
end

function system_subset(structure::SystemStructure, ieqs::Vector{Int})
    @unpack graph, eq_to_diff = structure
    fadj = Vector{Int}[]
    eq_to_diff = DiffGraph(length(ieqs))
    ne = 0
    for (j, eq_i) in enumerate(ieqs)
        ivars = copy(graph.fadjlist[eq_i])
        ne += length(ivars)
        push!(fadj, ivars)
        eq_to_diff[j] = structure.eq_to_diff[eq_i]
    end
    @set! structure.graph = complete(BipartiteGraph(ne, fadj, ndsts(graph)))
    @set! structure.eq_to_diff = eq_to_diff
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

function TearingState(sys; quick_cancel = false, check = true)
    sys = flatten(sys)
    ivs = independent_variables(sys)
    iv = length(ivs) == 1 ? ivs[1] : nothing
    # scalarize array equations, without scalarizing arguments to registered functions
    original_eqs = flatten_equations(copy(equations(sys)))
    eqs = copy(original_eqs)
    neqs = length(eqs)
    dervaridxs = OrderedSet{Int}()
    var2idx = Dict{Any, Int}()
    symbolic_incidence = []
    fullvars = []
    var_counter = Ref(0)
    var_types = VariableType[]
    addvar! = let fullvars = fullvars, var_counter = var_counter, var_types = var_types
        var -> get!(var2idx, var) do
            push!(fullvars, var)
            push!(var_types, getvariabletype(var))
            var_counter[] += 1
        end
    end

    vars = OrderedSet()
    varsvec = []
    for (i, eq′) in enumerate(eqs)
        if eq′.lhs isa Connection
            check ? error("$(nameof(sys)) has unexpanded `connect` statements") :
            return nothing
        end
        is_statemachine_equation = false
        if eq′.lhs isa StateMachineOperator
            is_statemachine_equation = true
            eq = eq′
            rhs = eq.rhs
        elseif _iszero(eq′.lhs)
            rhs = quick_cancel ? quick_cancel_expr(eq′.rhs) : eq′.rhs
            eq = eq′
        else
            lhs = quick_cancel ? quick_cancel_expr(eq′.lhs) : eq′.lhs
            rhs = quick_cancel ? quick_cancel_expr(eq′.rhs) : eq′.rhs
            eq = 0 ~ rhs - lhs
        end
        vars!(vars, eq.rhs, op = Symbolics.Operator)
        for v in vars
            _var, _ = var_from_nested_derivative(v)
            any(isequal(_var), ivs) && continue
            if isparameter(_var) ||
               (iscall(_var) && isparameter(operation(_var)) || isconstant(_var))
                continue
            end
            v = scalarize(v)
            if v isa AbstractArray
                v = setmetadata.(v, VariableIrreducible, true)
                append!(varsvec, v)
            else
                push!(varsvec, v)
            end
        end
        isalgeq = true
        unknownvars = []
        for var in varsvec
            ModelingToolkit.isdelay(var, iv) && continue
            set_incidence = true
            @label ANOTHER_VAR
            _var, _ = var_from_nested_derivative(var)
            any(isequal(_var), ivs) && continue
            if isparameter(_var) ||
               (iscall(_var) && isparameter(operation(_var)) || isconstant(_var))
                continue
            end
            varidx = addvar!(var)
            set_incidence && push!(unknownvars, var)

            dvar = var
            idx = varidx
            while isdifferential(dvar)
                if !(idx in dervaridxs)
                    push!(dervaridxs, idx)
                end
                isalgeq = false
                dvar = arguments(dvar)[1]
                idx = addvar!(dvar)
            end

            dvar = var
            idx = varidx

            if iscall(var) && operation(var) isa Symbolics.Operator &&
               !isdifferential(var) && (it = input_timedomain(var)) !== nothing
                set_incidence = false
                var = only(arguments(var))
                var = setmetadata(var, VariableTimeDomain, it)
                @goto ANOTHER_VAR
            end
        end
        push!(symbolic_incidence, copy(unknownvars))
        empty!(unknownvars)
        empty!(vars)
        empty!(varsvec)
        if isalgeq || is_statemachine_equation
            eqs[i] = eq
        else
            eqs[i] = eqs[i].lhs ~ rhs
        end
    end
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
            idx = addvar!(dvar)
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

    graph = BipartiteGraph(neqs, nvars, Val(false))
    for (ie, vars) in enumerate(symbolic_incidence), v in vars
        jv = var2idx[v]
        add_edge!(graph, ie, jv)
    end

    @set! sys.eqs = eqs

    eq_to_diff = DiffGraph(nsrcs(graph))

    ts = TearingState(sys, original_eqs, fullvars,
        SystemStructure(complete(var_to_diff), complete(eq_to_diff),
            complete(graph), nothing, var_types, sys isa DiscreteSystem),
        Any[], typeof(sys)[])
    if sys isa DiscreteSystem
        ts = shift_discrete_system(ts)
    end
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
        vars!(discvars, eq; op = Union{Sample, Hold})
    end
    iv = get_iv(sys)
    discmap = Dict(k => StructuralTransformations.simplify_shifts(Shift(iv, 1)(k))
    for k in discvars
    if any(isequal(k), fullvars) && !isa(operation(k), Union{Sample, Hold}))
    for i in eachindex(fullvars)
        fullvars[i] = StructuralTransformations.simplify_shifts(fast_substitute(
            fullvars[i], discmap; operator = Union{Sample, Hold}))
    end
    for i in eachindex(eqs)
        eqs[i] = StructuralTransformations.simplify_shifts(fast_substitute(
            eqs[i], discmap; operator = Union{Sample, Hold}))
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

# TODO: clean up
function merge_io(io, inputs)
    isempty(inputs) && return io
    if io === nothing
        io = (inputs, [])
    else
        io = ([inputs; io[1]], io[2])
    end
    return io
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

function structural_simplify!(state::TearingState, io = nothing; simplify = false,
        check_consistency = true, fully_determined = true, warn_initialize_determined = true,
        kwargs...)
    if state.sys isa ODESystem
        # split_system returns one or two systems and the inputs for each
        # mod clock inference to be binary
        # if it's continous keep going, if not then error unless given trait impl in additional passes
        ci = ModelingToolkit.ClockInference(state)
        ci = ModelingToolkit.infer_clocks!(ci)
        time_domains = merge(Dict(state.fullvars .=> ci.var_domain),
            Dict(default_toterm.(state.fullvars) .=> ci.var_domain))
        tss, inputs, continuous_id, id_to_clock = ModelingToolkit.split_system(ci)
        if continuous_id == 0
            # do a trait check here - handle fully discrete system
            additional_passes = get(kwargs, :additional_passes, nothing)
            if !isnothing(additional_passes) &&
               any(discrete_compile_pass, additional_passes)
                # take the first discrete compilation pass given for now
                discrete_pass_idx = findfirst(discrete_compile_pass, additional_passes)
                discrete_compile = additional_passes[discrete_pass_idx]
                deleteat!(additional_passes, discrete_pass_idx)
                return discrete_compile(tss, inputs)
            else
                # error goes here! this is a purely discrete system
                throw(HybridSystemNotSupportedException("Discrete systems without JuliaSimCompiler are currently not supported in ODESystem."))
            end
        end
        make_eqs_zero_equals!(tss[continuous_id])
        # puts the ios passed in to the call into the continous system
        cont_io = merge_io(io, inputs[continuous_id])
        # simplify as normal
        sys, input_idxs = _structural_simplify!(tss[continuous_id], cont_io; simplify,
            check_consistency, fully_determined,
            kwargs...)
        if length(tss) > 1
            additional_passes = get(kwargs, :additional_passes, nothing)
            if !isnothing(additional_passes) &&
               any(discrete_compile_pass, additional_passes)
                discrete_pass_idx = findfirst(discrete_compile_pass, additional_passes)
                discrete_compile = additional_passes[discrete_pass_idx]
                deleteat!(additional_passes, discrete_pass_idx)
                # in the case of a hybrid system, the discrete_compile pass should take the currents of sys.discrete_subsystems
                # and modifies discrete_subsystems to bea tuple of the io and anything else, while adding or manipulating the rest of sys as needed
                sys = discrete_compile(sys, tss[2:end], inputs)
            else
                throw(HybridSystemNotSupportedException("Hybrid continuous-discrete systems are currently not supported with the standard MTK compiler. This system requires JuliaSimCompiler.jl, see https://help.juliahub.com/juliasimcompiler/stable/"))
            end
        end
        ps = [sym isa CallWithMetadata ? sym :
              setmetadata(sym, VariableTimeDomain, get(time_domains, sym, Continuous()))
              for sym in get_ps(sys)]
        @set! sys.ps = ps
    else
        sys, input_idxs = _structural_simplify!(state, io; simplify, check_consistency,
            fully_determined, kwargs...)
    end
    has_io = io !== nothing
    return has_io ? (sys, input_idxs) : sys
end

function _structural_simplify!(state::TearingState, io; simplify = false,
        check_consistency = true, fully_determined = true, warn_initialize_determined = false,
        dummy_derivative = true,
        kwargs...)
    if fully_determined isa Bool
        check_consistency &= fully_determined
    else
        check_consistency = true
    end
    has_io = io !== nothing
    orig_inputs = Set()
    if has_io
        ModelingToolkit.markio!(state, orig_inputs, io...)
    end
    if io !== nothing
        state, input_idxs = ModelingToolkit.inputs_to_parameters!(state, io)
    else
        input_idxs = 0:-1 # Empty range
    end
    sys, mm = ModelingToolkit.alias_elimination!(state; kwargs...)
    if check_consistency
        fully_determined = ModelingToolkit.check_consistency(
            state, orig_inputs; nothrow = fully_determined === nothing)
    end
    if fully_determined && dummy_derivative
        sys = ModelingToolkit.dummy_derivative(
            sys, state; simplify, mm, check_consistency, kwargs...)
    elseif fully_determined
        var_eq_matching = pantelides!(state; finalize = false, kwargs...)
        sys = pantelides_reassemble(state, var_eq_matching)
        state = TearingState(sys)
        sys, mm = ModelingToolkit.alias_elimination!(state; kwargs...)
        sys = ModelingToolkit.dummy_derivative(
            sys, state; simplify, mm, check_consistency, kwargs...)
    else
        sys = ModelingToolkit.tearing(
            sys, state; simplify, mm, check_consistency, kwargs...)
    end
    fullunknowns = [map(eq -> eq.lhs, observed(sys)); unknowns(sys)]
    @set! sys.observed = ModelingToolkit.topsort_equations(observed(sys), fullunknowns)

    ModelingToolkit.invalidate_cache!(sys), input_idxs
end
