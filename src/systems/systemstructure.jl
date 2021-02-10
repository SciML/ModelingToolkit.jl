module SystemStructures

using DataStructures
using SymbolicUtils: istree, operation, arguments
using ..ModelingToolkit
import ..ModelingToolkit: isdiffeq, var_from_nested_derivative, vars!, flatten, value
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
export isdiffvar, isdervar, isalgvar, isdiffeq, isalgeq
export DIFFERENTIAL_VARIABLE, ALGEBRAIC_VARIABLE, DERIVATIVE_VARIABLE
export DIFFERENTIAL_EQUATION, ALGEBRAIC_EQUATION
export vartype, eqtype

struct SystemStructure
    dxvar_offset::Int
    fullvars::Vector # [diffvars; dervars; algvars]
    varassoc::Vector{Int}
    algeqs::BitVector
    graph::BipartiteGraph{Int,Nothing}
    solvable_graph::BipartiteGraph{Int,Vector{Vector{Int}}}
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

isalgeq(s::SystemStructure, eq::Integer) = s.algeqs[eq]
isdiffeq(s::SystemStructure, eq::Integer) = !isalgeq(s, eq)
eqtype(s::SystemStructure, eq::Integer)::EquationType = isalgeq(s, eq) ? ALGEBRAIC_EQUATION : DIFFERENTIAL_EQUATION

function initialize_system_structure(sys)
    sys, dxvar_offset, fullvars, varassoc, algeqs, graph, solvable_graph = init_graph(flatten(sys))
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

function init_graph(sys)
    iv = independent_variable(sys)
    eqs = equations(sys)
    neqs = length(eqs)
    algeqs = trues(neqs)
    varsadj = Vector{Any}(undef, neqs)
    dervars = OrderedSet()
    diffvars = OrderedSet()

    for (i, eq) in enumerate(eqs)
        vars = OrderedSet()
        vars!(vars, eq)
        isalgeq = true
        for var in vars
            if istree(var) && operation(var) isa Differential
                isalgeq = false
                diffvar = arguments(var)[1]
                @assert !(diffvar isa Differential) "The equation [ $eq ] is not first order"
                push!(dervars, var)
                push!(diffvars, diffvar)
            end
        end
        algeqs[i] = isalgeq
        varsadj[i] = vars
    end

    algvars  = setdiff(states(sys), diffvars)
    fullvars = [collect(diffvars); collect(dervars); algvars]

    dxvar_offset = length(diffvars)
    algvar_offset = 2dxvar_offset

    nvars = length(fullvars)
    idxmap = Dict(fullvars .=> 1:nvars)
    graph = BipartiteGraph(neqs, nvars)
    solvable_graph = BipartiteGraph(neqs, nvars, metadata=map(_->Int[], 1:neqs))

    for (i, vs) in enumerate(varsadj)
        eq = eqs[i]
        for v in vs
            j = get(idxmap, v, nothing)
            if j !== nothing
                add_edge!(graph, i, idxmap[v])
                j > algvar_offset || continue
                D = Differential(fullvars[j])
                c = value(expand_derivatives(D(eq.rhs - eq.lhs), false))
                if c isa Number && c != 0
                    # 0 here is a sentinel value for non-integer coefficients
                    coeff = c isa Integer ? c : 0
                    add_edge!(solvable_graph, i, j, coeff)
                end
            end
        end
    end

    varassoc = Int[(1:dxvar_offset) .+ dxvar_offset; zeros(Int, length(fullvars) - dxvar_offset)] # variable association list
    sys, dxvar_offset, fullvars, varassoc, algeqs, graph, solvable_graph
end

end # module
