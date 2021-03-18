module SystemStructures

using DataStructures
using SymbolicUtils: istree, operation, arguments, Symbolic
using ..ModelingToolkit
import ..ModelingToolkit: isdiffeq, var_from_nested_derivative, vars!, flatten,
    value, InvalidSystemException, isdifferential
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

export SystemStructure, initialize_system_structure, find_linear_equations
export isdiffvar, isdervar, isalgvar, isdiffeq, isalgeq

@enum VariableType::Int8 DIFFERENTIAL_VARIABLE ALGEBRAIC_VARIABLE DERIVATIVE_VARIABLE

Base.@kwdef struct SystemStructure
    fullvars::Vector
    vartype::Vector{VariableType}
    inv_varassoc::Vector{Int}
    varassoc::Vector{Int}
    algeqs::BitVector
    graph::BipartiteGraph{Int,Nothing}
    solvable_graph::BipartiteGraph{Int,Nothing}
    assign::Vector{Int}
    inv_assign::Vector{Int}
    scc::Vector{Vector{Int}}
    partitions::Vector{NTuple{4, Vector{Int}}}
end

isdervar(s::SystemStructure, var::Integer) = s.vartype[var] === DERIVATIVE_VARIABLE
isdiffvar(s::SystemStructure, var::Integer) = s.vartype[var] === DIFFERENTIAL_VARIABLE
isalgvar(s::SystemStructure, var::Integer) = s.vartype[var] === ALGEBRAIC_VARIABLE

dervars_range(s::SystemStructure) = Iterators.filter(Base.Fix1(s, isdervar), eachindex(s.vartype))
diffvars_range(s::SystemStructure) = Iterators.filter(Base.Fix1(s, isdiffvar), eachindex(s.vartype))
algvars_range(s::SystemStructure) = Iterators.filter(Base.Fix1(s, isalgeq), eachindex(s.vartype))

isalgeq(s::SystemStructure, eq::Integer) = s.algeqs[eq]
isdiffeq(s::SystemStructure, eq::Integer) = !isalgeq(s, eq)

function initialize_system_structure(sys)
    iv = independent_variable(sys)
    eqs = equations(sys)
    neqs = length(eqs)
    algeqs = trues(neqs)
    dervaridxs = Int[]
    var2idx = Dict{Any,Int}()
    symbolic_incidence = []
    fullvars = []
    var_counter = 0

    for (i, eq) in enumerate(eqs)
        vars = OrderedSet()
        vars!(vars, eq)
        push!(symbolic_incidence, copy(vars))
        isalgeq = true
        for var in vars
            varidx = get(var2idx, var, 0)
            if varidx == 0 # new var
                var_counter += 1
                push!(fullvars, var)
            end

            if isdifferential(var)
                isalgeq = false
                diffvar = arguments(var)[1]
                if diffvar isa Differential
                    throw(InvalidSystemException("The equation [ $eq ] is not first order"))
                end
                push!(dervaridxs, varidx)
            end
        end
        algeqs[i] = isalgeq
    end

    diffvars = []
    varassoc = zeros(Int, length(fullvars))
    inv_varassoc = zeros(Int, length(fullvars))
    for dervaridx in dervaridxs
        dervar = fullvars[dervaridx]
        diffvar = arguments(dervar)[1]
        diffvaridx = get(var2idx, diffvar, 0)
        if diffvaridx != 0
            push!(diffvars, diffvar)
            varassoc[diffvaridx] = dervaridx
            inv_varassoc[dervaridx] = diffvaridx
        end
    end

    algvars = setdiff(states(sys), diffvars)
    for algvar in algvars
        # it could be that a variable appeared in the states, but never appeared
        # in the equations.
        algvaridx = get(var2idx, algvar, 0)
        if algvaridx != 0
            varassoc[algvaridx] = -1
        end
    end

    neqs, nvars = length(eqs), length(fullvars)
    graph = BipartiteGraph(neqs, nvars)
    for (ie, vars) in enumerate(symbolic_incidence), v in vars
        jv = var2idx[v]
        add_edge!(graph, ie, jv)
    end

    SystemStructure(
        fullvars = fullvars,
        varassoc = varassoc,
        inv_varassoc = inv_varassoc,
        algeqs = algeqs,
        graph = graph,
        solvable_graph = BipartiteGraph(nsrcs(graph), ndsts(graph)),
        assign = Int[],
        inv_assign = Int[],
        scc = Vector{Int}[],
        partitions = NTuple{4, Vector{Int}}[],
    )
end

function find_linear_equations(sys)
    s = structure(sys)
    @unpack fullvars, graph = s
    is_linear_equations = falses(ndsts(graph))
    eqs = equations(sys)
    eadj = Vector{Int}[]
    cadj = Vector{Int}[]
    coeffs = Int[]
    for (i, eq) in enumerate(eqs); isdiffeq(eq) && continue
        empty!(coeffs)
        linear_term = 0
        all_int_algvars = true

        term = value(eq.rhs - eq.lhs)
        for j in 𝑠neighbors(graph, i)
            if !isalgvar(s, j)
                all_int_algvars = false
                continue
            end
            var = fullvars[j]
            c = expand_derivatives(Differential(var)(term), false)
            # test if `var` is linear in `eq`.
            if !(c isa Symbolic) && c isa Number
                if isinteger(c) && !iszero(c)
                    c = convert(Integer, c)
                    linear_term += c * var
                    push!(coeffs, c)
                else
                    all_int_algvars = false
                end
            end
        end

        # Check if there are only algebraic variables and the equation is both
        # linear and homogeneous, i.e. it is in the form of
        #
        #       ``∑ c_i * a_i = 0``,
        #
        # where ``c_i`` ∈ ℤ and ``a_i`` denotes algebraic variables.
        if all_int_algvars && isequal(linear_term, term)
            is_linear_equations[i] = true
            push!(eadj, copy(𝑠neighbors(graph, i)))
            push!(cadj, copy(coeffs))
        else
            is_linear_equations[i] = false
        end
    end
    sys, dxvar_offset, fullvars, varassoc, algeqs, graph = init_graph(flatten(sys))

    solvable_graph = BipartiteGraph(nsrcs(graph), ndsts(graph))

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
    return is_linear_equations, eadj, cadj
end

end # module
