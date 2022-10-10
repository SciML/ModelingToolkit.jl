module SystemStructures

using DataStructures
using Symbolics: linear_expansion, unwrap
using SymbolicUtils: istree, operation, arguments, Symbolic
using SymbolicUtils: quick_cancel, similarterm
using ..ModelingToolkit
import ..ModelingToolkit: isdiffeq, var_from_nested_derivative, vars!, flatten,
                          value, InvalidSystemException, isdifferential, _iszero,
                          isparameter,
                          independent_variables, SparseMatrixCLIL, AbstractSystem,
                          equations, isirreducible
using ..BipartiteGraphs
import ..BipartiteGraphs: invview, complete
using Graphs
using UnPack
using Setfield
using SparseArrays

function quick_cancel_expr(expr)
    Rewriters.Postwalk(quick_cancel,
                       similarterm = (x, f, args; kws...) -> similarterm(x, f, args,
                                                                         SymbolicUtils.symtype(x);
                                                                         metadata = SymbolicUtils.metadata(x),
                                                                         kws...))(expr)
end

export SystemStructure, TransformationState, TearingState
export initialize_system_structure, find_linear_equations
export isdiffvar, isdervar, isalgvar, isdiffeq, isalgeq, algeqs
export dervars_range, diffvars_range, algvars_range
export DiffGraph, complete!

@enum VariableType::Int8 DIFFERENTIAL_VARIABLE ALGEBRAIC_VARIABLE DERIVATIVE_VARIABLE

struct DiffGraph <: Graphs.AbstractGraph{Int}
    primal_to_diff::Vector{Union{Int, Nothing}}
    diff_to_primal::Union{Nothing, Vector{Union{Int, Nothing}}}
end
DiffGraph(primal_to_diff::Vector{Union{Int, Nothing}}) = DiffGraph(primal_to_diff, nothing)
function DiffGraph(n::Integer, with_badj::Bool = false)
    DiffGraph(Union{Int, Nothing}[nothing for _ in 1:n],
              with_badj ? Union{Int, Nothing}[nothing for _ in 1:n] : nothing)
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
end
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
end

mutable struct TearingState{T <: AbstractSystem} <: AbstractTearingState{T}
    sys::T
    fullvars::Vector
    structure::SystemStructure
    extra_eqs::Vector
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

function TearingState(sys; quick_cancel = false, check = true)
    sys = flatten(sys)
    ivs = independent_variables(sys)
    eqs = copy(equations(sys))
    neqs = length(eqs)
    algeqs = trues(neqs)
    dervaridxs = OrderedSet{Int}()
    var2idx = Dict{Any, Int}()
    symbolic_incidence = []
    fullvars = []
    var_counter = Ref(0)
    addvar! = let fullvars = fullvars, var_counter = var_counter
        var -> begin get!(var2idx, var) do
            push!(fullvars, var)
            var_counter[] += 1
        end end
    end

    vars = OrderedSet()
    for (i, eq′) in enumerate(eqs)
        if eq′.lhs isa Connection
            check ? error("$(nameof(sys)) has unexpanded `connect` statements") :
            return nothing
        end
        if _iszero(eq′.lhs)
            rhs = quick_cancel ? quick_cancel_expr(eq′.rhs) : eq′.rhs
            eq = eq′
        else
            lhs = quick_cancel ? quick_cancel_expr(eq′.lhs) : eq′.lhs
            rhs = quick_cancel ? quick_cancel_expr(eq′.rhs) : eq′.rhs
            eq = 0 ~ rhs - lhs
        end
        vars!(vars, eq.rhs)
        isalgeq = true
        statevars = []
        for var in vars
            _var, _ = var_from_nested_derivative(var)
            any(isequal(_var), ivs) && continue
            if isparameter(_var) || (istree(_var) && isparameter(operation(_var)))
                continue
            end
            varidx = addvar!(var)
            push!(statevars, var)

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
        end
        push!(symbolic_incidence, copy(statevars))
        empty!(statevars)
        empty!(vars)
        algeqs[i] = isalgeq
        if isalgeq
            eqs[i] = eq
        else
            eqs[i] = eqs[i].lhs ~ rhs
        end
    end

    # sort `fullvars` such that the mass matrix is as diagonal as possible.
    dervaridxs = collect(dervaridxs)
    sorted_fullvars = OrderedSet(fullvars[dervaridxs])
    for dervaridx in dervaridxs
        dervar = fullvars[dervaridx]
        diffvar = arguments(dervar)[1]
        if !(diffvar in sorted_fullvars)
            push!(sorted_fullvars, diffvar)
        end
    end
    for v in fullvars
        if !(v in sorted_fullvars)
            push!(sorted_fullvars, v)
        end
    end
    fullvars = collect(sorted_fullvars)
    var2idx = Dict(fullvars .=> eachindex(fullvars))
    dervaridxs = 1:length(dervaridxs)

    nvars = length(fullvars)
    diffvars = []
    var_to_diff = DiffGraph(nvars, true)
    for dervaridx in dervaridxs
        dervar = fullvars[dervaridx]
        diffvar = arguments(dervar)[1]
        diffvaridx = var2idx[diffvar]
        push!(diffvars, diffvar)
        var_to_diff[diffvaridx] = dervaridx
    end

    #=
    algvars = setdiff(states(sys), diffvars)
    for algvar in algvars
        # it could be that a variable appeared in the states, but never appeared
        # in the equations.
        algvaridx = get(var2idx, algvar, 0)
        #if algvaridx == 0
        #    check ? throw(InvalidSystemException("The system is missing an equation for $algvar.")) : return nothing
        #end
    end
    =#

    graph = BipartiteGraph(neqs, nvars, Val(false))
    for (ie, vars) in enumerate(symbolic_incidence), v in vars
        jv = var2idx[v]
        add_edge!(graph, ie, jv)
    end

    @set! sys.eqs = eqs

    eq_to_diff = DiffGraph(nsrcs(graph))

    return TearingState(sys, fullvars,
                        SystemStructure(complete(var_to_diff), complete(eq_to_diff),
                                        complete(graph), nothing), Any[])
end

using .BipartiteGraphs: Label, BipartiteAdjacencyList
struct SystemStructurePrintMatrix <:
       AbstractMatrix{Union{Label, Int, BipartiteAdjacencyList}}
    bpg::BipartiteGraph
    highlight_graph::BipartiteGraph
    var_to_diff::DiffGraph
    eq_to_diff::DiffGraph
    var_eq_matching::Union{Matching, Nothing}
end
Base.size(bgpm::SystemStructurePrintMatrix) = (max(nsrcs(bgpm.bpg), ndsts(bgpm.bpg)) + 1, 5)
function compute_diff_label(diff_graph, i)
    di = i - 1 <= length(diff_graph) ? diff_graph[i - 1] : nothing
    ii = i - 1 <= length(invview(diff_graph)) ? invview(diff_graph)[i - 1] : nothing
    return Label(string(di === nothing ? "" : string(di, '↓'),
                        di !== nothing && ii !== nothing ? " " : "",
                        ii === nothing ? "" : string(ii, '↑')))
end
function Base.getindex(bgpm::SystemStructurePrintMatrix, i::Integer, j::Integer)
    checkbounds(bgpm, i, j)
    if i <= 1
        return (Label.(("#", "∂ₜ", "eq", "∂ₜ", "v")))[j]
    elseif j == 2
        return compute_diff_label(bgpm.eq_to_diff, i)
    elseif j == 4
        return compute_diff_label(bgpm.var_to_diff, i)
    elseif j == 1
        return i - 1
    elseif j == 3
        return BipartiteAdjacencyList(i - 1 <= nsrcs(bgpm.bpg) ?
                                      𝑠neighbors(bgpm.bpg, i - 1) : nothing,
                                      bgpm.highlight_graph !== nothing &&
                                      i - 1 <= nsrcs(bgpm.highlight_graph) ?
                                      Set(𝑠neighbors(bgpm.highlight_graph, i - 1)) :
                                      nothing,
                                      bgpm.var_eq_matching !== nothing &&
                                      (i - 1 <= length(invview(bgpm.var_eq_matching))) ?
                                      invview(bgpm.var_eq_matching)[i - 1] : unassigned)
    elseif j == 5
        match = unassigned
        if bgpm.var_eq_matching !== nothing && i - 1 <= length(bgpm.var_eq_matching)
            match = bgpm.var_eq_matching[i - 1]
            isa(match, Union{Int, Unassigned}) || (match = true) # Selected State
        end
        return BipartiteAdjacencyList(i - 1 <= ndsts(bgpm.bpg) ?
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
        print(io, "SystemStructure with ", length(graph.fadjlist), " equations and ",
              isa(graph.badjlist, Int) ? graph.badjlist : length(graph.badjlist),
              " variables\n")
        Base.print_matrix(io,
                          SystemStructurePrintMatrix(complete(graph),
                                                     complete(solvable_graph),
                                                     complete(var_to_diff),
                                                     complete(eq_to_diff), nothing))
    else
        S = incidence_matrix(graph, Num(Sym{Real}(:×)))
        print(io, "Incidence matrix:")
        show(io, mime, S)
    end
end

struct MatchedSystemStructure
    structure::SystemStructure
    var_eq_matching::Matching
end

function Base.show(io::IO, mime::MIME"text/plain", ms::MatchedSystemStructure)
    s = ms.structure
    @unpack graph, solvable_graph, var_to_diff, eq_to_diff = s
    print(io, "Matched SystemStructure with ", length(graph.fadjlist), " equations and ",
          isa(graph.badjlist, Int) ? graph.badjlist : length(graph.badjlist),
          " variables\n")
    Base.print_matrix(io,
                      SystemStructurePrintMatrix(complete(graph),
                                                 complete(solvable_graph),
                                                 complete(var_to_diff),
                                                 complete(eq_to_diff),
                                                 complete(ms.var_eq_matching, nsrcs(graph))))
end

end # module
