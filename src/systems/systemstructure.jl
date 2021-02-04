module SystemStructures

using ..ModelingToolkit
using ..ModelingToolkit: isdiffeq, var_from_nested_derivative, vars!
using SymbolicUtils: arguments
using ..BipartiteGraphs
using UnPack
using Setfield
using SparseArrays

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

export SystemStructure, initialize_system_structure
export diffvars_range, dervars_range, algvars_range
export isdiffvars, isdervars, isalgvars
export DIFFERENTIAL_VARIABLE, ALGEBRAIC_VARIABLE, DERIVATIVE_VARIABLE
export DIFFERENTIAL_EQUATION, ALGEBRAIC_EQUATION
export vartype, eqtype

struct SystemStructure
    dxvar_offset::Int
    fullvars::Vector # [xvar; dxvars; algvars]
    varassoc::Vector{Int}
    algeqs::BitVector
    graph::BipartiteGraph{Int}
    solvable_graph::BipartiteGraph{Int}
    assign::Vector{Int}
    inv_assign::Vector{Int}
    scc::Vector{Vector{Int}}
    partitions::Vector{NTuple{4, Vector{Int}}}
end

diffvars_range(s::SystemStructure) = 1:s.dxvar_offset
dervars_range(s::SystemStructure) = s.dxvar_offset+1:2s.dxvar_offset
algvars_range(s::SystemStructure) = 2s.dxvar_offset+1:length(s.fullvars)

isdiffvar(s::SystemStructure, var::Integer) = var in diffvars_range(s)
isdervar(s::SystemStructure, var::Integer) = var in dervars_range(s)
isalgvar(s::SystemStructure, var::Integer) = var in algvars_range(s)

@enum VariableType DIFFERENTIAL_VARIABLE ALGEBRAIC_VARIABLE DERIVATIVE_VARIABLE

function vartype(s::SystemStructure, var::Integer)::VariableType
    isdiffvar(s, var) ? DIFFERENTIAL_VARIABLE :
    isdervar(s, var)  ? DERIVATIVE_VARIABLE :
    isalgvar(s, var)  ? ALGEBRAIC_VARIABLE : error("Variable $var out of bounds")
end

@enum EquationType DIFFERENTIAL_EQUATION ALGEBRAIC_EQUATION

eqtype(s::SystemStructure, eq::Integer)::EquationType = s.algeqs[eq] ? ALGEBRAIC_EQUATION : DIFFERENTIAL_EQUATION

function initialize_system_structure(sys)
    sys, dxvar_offset, fullvars, varassoc, algeqs, graph, solvable_graph = init_graph(sys)
    @set sys.structure = SystemStructure(
                                         dxvar_offset,
                                         fullvars,
                                         varassoc,
                                         algeqs,
                                         graph,
                                         solvable_graph,
                                         Int[],
                                         Int[],
                                         Vector{Int}[],
                                         NTuple{4, Vector{Int}}[]
                                        )
end

function Base.show(io::IO, s::SystemStructure)
    @unpack fullvars, dxvar_offset, solvable_graph, graph = s
    algvar_offset = 2dxvar_offset
    print(io, "xvars: ")
    print(io, fullvars[1:dxvar_offset])
    print(io, "\ndxvars: ")
    print(io, fullvars[dxvar_offset+1:algvar_offset])
    print(io, "\nalgvars: ")
    print(io, fullvars[algvar_offset+1:end], '\n')

    S = incidence_matrix(graph, Num(Sym{Real}(:×)))
    print(io, "Incidence matrix:")
    show(io, S)
end

# V-nodes `[x_1, x_2, x_3, ..., dx_1, dx_2, ..., y_1, y_2, ...]` where `x`s are
# differential variables and `y`s are algebraic variables.
function collect_variables(sys)
    dxvars = []
    eqs = equations(sys)
    algeqs = falses(length(eqs))
    for (i, eq) in enumerate(eqs)
        if isdiffeq(eq)
            algeqs[i] = true
            lhs = eq.lhs
            # Make sure that the LHS is a first order derivative of a var.
            @assert !(arguments(lhs)[1] isa Differential) "The equation $eq is not first order"

            push!(dxvars, lhs)
        end
    end

    xvars = (first ∘ var_from_nested_derivative).(dxvars)
    algvars  = setdiff(states(sys), xvars)
    return xvars, dxvars, algvars, algeqs
end

function init_graph(sys)
    xvars, dxvars, algvars, algeqs = collect_variables(sys)
    dxvar_offset = length(xvars)
    algvar_offset = 2dxvar_offset

    fullvars = [xvars; dxvars; algvars]
    eqs = equations(sys)
    idxmap = Dict(fullvars .=> 1:length(fullvars))
    graph = BipartiteGraph(length(eqs), length(fullvars))
    solvable_graph = BipartiteGraph(length(eqs), length(fullvars))

    vs = Set()
    for (i, eq) in enumerate(eqs)
        # TODO: custom vars that handles D(x)
        # TODO: add checks here
        lhs = eq.lhs
        if isdiffeq(eq)
            v = lhs
            haskey(idxmap, v) && add_edge!(graph, i, idxmap[v])
        else
            vars!(vs, lhs)
        end
        vars!(vs, eq.rhs)
        for v in vs
            haskey(idxmap, v) && add_edge!(graph, i, idxmap[v])
        end
        empty!(vs)
    end

    varassoc = Int[(1:dxvar_offset) .+ dxvar_offset; zeros(Int, length(fullvars) - dxvar_offset)] # variable association list
    sys, dxvar_offset, fullvars, varassoc, algeqs, graph, solvable_graph
end

end # module
