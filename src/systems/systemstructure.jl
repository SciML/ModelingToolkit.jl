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
        similarterm = (x, f, args; kws...) -> maketerm(typeof(x), f, args,
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
end

TransformationState(sys::AbstractSystem) = TearingState(sys)
function system_subset(ts::TearingState, ieqs::Vector{Int})
    eqs = equations(ts)
    @set! ts.sys.eqs = eqs[ieqs]
    @set! ts.structure = system_subset(ts.structure, ieqs)
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

function TearingState(sys; quick_cancel = false, check = true, sort_eqs = true)
    # flatten system
    sys = flatten(sys)
    sys = process_parameter_equations(sys)
    ivs = independent_variables(sys)
    iv = length(ivs) == 1 ? ivs[1] : nothing
    # flatten array equations
    eqs = flatten_equations(equations(sys))
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
        rhs = quick_cancel ? quick_cancel_expr(eq.rhs) : eq.rhs
        if !_iszero(eq.lhs)
            lhs = quick_cancel ? quick_cancel_expr(eq.lhs) : eq.lhs
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
                isvalid = iscall(v) && operation(v) isa Union{Shift, Sample, Hold}
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
                    v′ = only(arguments(v))
                    addvar!(setmetadata(v′, VariableTimeDomain, it), VARIABLE)
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

        if isalgeq
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
    neqs = length(eqs)
    symbolic_incidence = symbolic_incidence[eqs_to_retain]

    if sort_eqs
        # sort equations lexicographically to reduce simplification issues
        # depending on order due to NP-completeness of tearing.
        sortidxs = Base.sortperm(eqs, by = string)
        eqs = eqs[sortidxs]
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
        Any[], param_derivative_map)

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

function mtkcompile!(state::TearingState; simplify = false,
        check_consistency = true, fully_determined = true, warn_initialize_determined = true,
        inputs = Any[], outputs = Any[],
        disturbance_inputs = Any[],
        kwargs...)
    ci = ModelingToolkit.ClockInference(state)
    ci = ModelingToolkit.infer_clocks!(ci)
    time_domains = merge(Dict(state.fullvars .=> ci.var_domain),
        Dict(default_toterm.(state.fullvars) .=> ci.var_domain))
    tss, clocked_inputs, continuous_id, id_to_clock = ModelingToolkit.split_system(ci)
    if length(tss) > 1
        if continuous_id == 0
            throw(HybridSystemNotSupportedException("""
            Discrete systems with multiple clocks are not supported with the standard \
            MTK compiler.
            """))
        else
            throw(HybridSystemNotSupportedException("""
            Hybrid continuous-discrete systems are currently not supported with \
            the standard MTK compiler. This system requires JuliaSimCompiler.jl, \
            see https://help.juliahub.com/juliasimcompiler/stable/
            """))
        end
    end
    if continuous_id == 1 && any(Base.Fix2(isoperator, Shift), state.fullvars)
        state.structure.only_discrete = true
    end

    sys = _mtkcompile!(state; simplify, check_consistency,
        inputs, outputs, disturbance_inputs,
        fully_determined, kwargs...)
    return sys
end

function _mtkcompile!(state::TearingState; simplify = false,
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
    has_io = !isempty(inputs) || !isempty(outputs) !== nothing ||
             !isempty(disturbance_inputs)
    orig_inputs = Set()
    if has_io
        ModelingToolkit.markio!(state, orig_inputs, inputs, outputs, disturbance_inputs)
        state = ModelingToolkit.inputs_to_parameters!(state, [inputs; disturbance_inputs])
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
