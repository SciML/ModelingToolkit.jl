module SystemStructures

using DataStructures
using SymbolicUtils: istree, operation, arguments, Symbolic
using ..ModelingToolkit
import ..ModelingToolkit: isdiffeq, var_from_nested_derivative, vars!, flatten,
    value, InvalidSystemException, isdifferential, _iszero, isparameter
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
export dervars_range, diffvars_range, algvars_range

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

dervars_range(s::SystemStructure) = Iterators.filter(Base.Fix1(isdervar, s), eachindex(s.vartype))
diffvars_range(s::SystemStructure) = Iterators.filter(Base.Fix1(isdiffvar, s), eachindex(s.vartype))
algvars_range(s::SystemStructure) = Iterators.filter(Base.Fix1(isalgvar, s), eachindex(s.vartype))

isalgeq(s::SystemStructure, eq::Integer) = s.algeqs[eq]
isdiffeq(s::SystemStructure, eq::Integer) = !isalgeq(s, eq)

function initialize_system_structure(sys)
    sys = flatten(sys)

    iv = independent_variable(sys)
    eqs = copy(equations(sys))
    neqs = length(eqs)
    algeqs = trues(neqs)
    dervaridxs = Int[]
    var2idx = Dict{Any,Int}()
    symbolic_incidence = []
    fullvars = []
    var_counter = Ref(0)
    addvar! = let fullvars=fullvars, var_counter=var_counter
        var -> begin
            get!(var2idx, var) do
                push!(fullvars, var)
                var_counter[] += 1
            end
        end
    end

    for (i, eq) in enumerate(eqs)
        vars = OrderedSet()
        vars!(vars, eq)
        isalgeq = true
        statevars = []
        for var in vars
            isequal(var, iv) && continue
            if isparameter(var) || (istree(var) && isparameter(operation(var)))
                continue
            end
            varidx = addvar!(var)
            push!(statevars, var)

            dvar = var
            idx = varidx
            while isdifferential(dvar)
                push!(dervaridxs, idx)
                isalgeq = false
                dvar = arguments(dvar)[1]
                idx = addvar!(dvar)
            end
        end
        push!(symbolic_incidence, copy(statevars))
        empty!(statevars)
        algeqs[i] = isalgeq
        if isalgeq && !_iszero(eq.lhs)
            eqs[i] = 0 ~ eq.rhs - eq.lhs
        end
    end

    nvars = length(fullvars)
    diffvars = []
    vartype = fill(DIFFERENTIAL_VARIABLE, nvars)
    varassoc = zeros(Int, nvars)
    inv_varassoc = zeros(Int, nvars)
    for dervaridx in dervaridxs
        vartype[dervaridx] = DERIVATIVE_VARIABLE
        dervar = fullvars[dervaridx]
        diffvar = arguments(dervar)[1]
        diffvaridx = var2idx[diffvar]
        push!(diffvars, diffvar)
        varassoc[diffvaridx] = dervaridx
        inv_varassoc[dervaridx] = diffvaridx
    end

    algvars = setdiff(states(sys), diffvars)
    for algvar in algvars
        # it could be that a variable appeared in the states, but never appeared
        # in the equations.
        algvaridx = get(var2idx, algvar, 0)
        vartype[algvaridx] = ALGEBRAIC_VARIABLE
    end

    graph = BipartiteGraph(neqs, nvars)
    for (ie, vars) in enumerate(symbolic_incidence), v in vars
        jv = var2idx[v]
        add_edge!(graph, ie, jv)
    end

    @set! sys.eqs = eqs
    @set! sys.structure = SystemStructure(
        fullvars = fullvars,
        vartype = vartype,
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
    return sys
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

    return is_linear_equations, eadj, cadj
end

function Base.show(io::IO, s::SystemStructure)
    @unpack graph = s
    S = incidence_matrix(graph, Num(Sym{Real}(:×)))
    print(io, "Incidence matrix:")
    show(io, S)
end

end # module
