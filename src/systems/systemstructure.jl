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

#=
When we don't do subsitution, variable information is split into two different
places, i.e. `states` and the right-hand-side of `observed`.

eqs = [0 ~ z + x; 0 ~ y + z^2]
states = [y, z]
observed = [x ~ sin(y) + z]
struct Reduced
    var
    expr
    idxs
end
fullvars = [Reduced(x, sin(y) + z, [2, 3]), y, z]
active_𝑣vertices = [false, true, true]
      x   y   z
eq1:  1       1
eq2:      1   1

      x   y   z
eq1:      1   1
eq2:      1   1

for v in 𝑣vertices(graph); active_𝑣vertices[v] || continue

end
=#

export SystemStructure, TransformationState, TearingState
export initialize_system_structure, find_linear_equations, linear_subsys_adjmat
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
    vartype = fill(DIFFERENTIAL_VARIABLE, nvars)
    var_to_diff = DiffGraph(nvars, true)
    for dervaridx in dervaridxs
        vartype[dervaridx] = DERIVATIVE_VARIABLE
        dervar = fullvars[dervaridx]
        diffvar = arguments(dervar)[1]
        diffvaridx = var2idx[diffvar]
        push!(diffvars, diffvar)
        var_to_diff[diffvaridx] = dervaridx
    end

    algvars = setdiff(states(sys), diffvars)
    for algvar in algvars
        # it could be that a variable appeared in the states, but never appeared
        # in the equations.
        algvaridx = get(var2idx, algvar, 0)
        #if algvaridx == 0
        #    check ? throw(InvalidSystemException("The system is missing an equation for $algvar.")) : return nothing
        #end
        if algvaridx != 0
            vartype[algvaridx] = ALGEBRAIC_VARIABLE
        end
    end

    graph = BipartiteGraph(neqs, nvars, Val(false))
    for (ie, vars) in enumerate(symbolic_incidence), v in vars
        jv = var2idx[v]
        add_edge!(graph, ie, jv)
    end

    @set! sys.eqs = eqs

    eq_to_diff = DiffGraph(nsrcs(graph))

    return TearingState(sys, fullvars,
                        SystemStructure(var_to_diff, eq_to_diff, graph, nothing), Any[])
end

function linear_subsys_adjmat(state::TransformationState)
    fullvars = state.fullvars
    graph = state.structure.graph
    is_linear_equations = falses(nsrcs(graph))
    eqs = equations(state.sys)
    eadj = Vector{Int}[]
    cadj = Vector{Int}[]
    coeffs = Int[]
    for (i, eq) in enumerate(eqs)
        empty!(coeffs)
        linear_term = 0
        all_int_vars = true

        term = value(eq.rhs - eq.lhs)
        for j in 𝑠neighbors(graph, i)
            var = fullvars[j]
            a, b, islinear = linear_expansion(term, var)
            a = unwrap(a)
            if islinear && !(a isa Symbolic) && a isa Number && !isirreducible(var)
                if a == 1 || a == -1
                    a = convert(Integer, a)
                    linear_term += a * var
                    push!(coeffs, a)
                else
                    all_int_vars = false
                end
            end
        end

        # Check if all states in the equation is both linear and homogeneous,
        # i.e. it is in the form of
        #
        #       ``∑ c_i * v_i = 0``,
        #
        # where ``c_i`` ∈ ℤ and ``v_i`` denotes states.
        if all_int_vars && isequal(linear_term, term)
            is_linear_equations[i] = true
            push!(eadj, copy(𝑠neighbors(graph, i)))
            push!(cadj, copy(coeffs))
        else
            is_linear_equations[i] = false
        end
    end

    linear_equations = findall(is_linear_equations)
    return SparseMatrixCLIL(nsrcs(graph),
                            ndsts(graph),
                            linear_equations, eadj, cadj)
end

function Base.show(io::IO, mime::MIME"text/plain", s::SystemStructure)
    @unpack graph = s
    S = incidence_matrix(graph, Num(Sym{Real}(:×)))
    print(io, "Incidence matrix:")
    show(io, mime, S)
end

end # module
