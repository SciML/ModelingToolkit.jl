using SparseArrays

const SHOW_EQUATIONS = Ref(false)
struct SystemStructure
    sys
    dxvar_offset::Int
    fullvars::Vector # [xvar; dxvars; algvars]
    varassoc::Vector{Int}
    graph::BipartiteGraph{Int}
    solvable_graph::BipartiteGraph{Int}
end
function SystemStructure(sys)
    sys = ModelingToolkit.flatten(sys)
    sys, dxvar_offset, fullvars, varassoc, graph, solvable_graph = init_graph(sys)
    SystemStructure(sys, dxvar_offset, fullvars, varassoc, graph, solvable_graph)
end

ModelingToolkit.equations(s::SystemStructure) = equations(s.sys)

function Base.show(io::IO, s::SystemStructure)
    @unpack fullvars, dxvar_offset, solvable_graph, graph = s
    algvar_offset = 2dxvar_offset
    print(io, "xvars: ")
    print(io, fullvars[1:dxvar_offset])
    print(io, "\ndxvars: ")
    print(io, fullvars[dxvar_offset+1:algvar_offset])
    print(io, "\nalgvars: ")
    print(io, fullvars[algvar_offset+1:end], '\n')

    if SHOW_EQUATIONS[]
        println(io, "Edges:")
        eqs = equations(s)
        for ev in 𝑠vertices(graph)
            print(io, "  $(eqs[ev])\n    -> ")
            vars = 𝑠neighbors(graph, ev)
            solvars = 𝑠neighbors(solvable_graph, ev)
            solvable = intersect(vars, solvars)
            notsolvable = setdiff(vars, solvars)

            print(io, join(string.(fullvars[notsolvable]), ", "))
            for ii in solvable
                print(io, ", ")
                var = fullvars[ii]
                Base.printstyled(io, string(fullvars[ii]), color=:cyan)
            end
            println(io)
        end
    end

    S = incidence_matrix(graph, Num(Sym{Real}(:×)))
    print(io, "Incidence matrix:")
    show(io, S)
end

# V-nodes `[x_1, x_2, x_3, ..., dx_1, dx_2, ..., y_1, y_2, ...]` where `x`s are
# differential variables and `y`s are algebraic variables.
function get_vnodes(sys)
    dxvars = []
    eqs = equations(sys)
    for (i, eq) in enumerate(eqs)
        if eq.lhs isa Symbolic
            # Make sure that the LHS is a first order derivative of a var.
            @assert operation(eq.lhs) isa Differential "The equation $eq is not in the form of `D(...) ~ ...`"
            @assert !(arguments(eq.lhs)[1] isa Differential) "The equation $eq is not first order"

            push!(dxvars, eq.lhs)
        end
    end

    xvars = (first ∘ var_from_nested_derivative).(dxvars)
    algvars  = setdiff(states(sys), xvars)
    return xvars, dxvars, algvars
end

function init_graph(sys)
    xvars, dxvars, algvars = get_vnodes(sys)
    dxvar_offset = length(xvars)
    algvar_offset = 2dxvar_offset

    fullvars = [xvars; dxvars; algvars]
    sys = reordersys(sys, dxvar_offset, fullvars)
    eqs = equations(sys)
    idxmap = Dict(fullvars .=> 1:length(fullvars))
    graph = BipartiteGraph(length(eqs), length(fullvars))
    solvable_graph = BipartiteGraph(length(eqs), length(fullvars))

    for (i, eq) in enumerate(eqs)
        if isdiffeq(eq)
            v = eq.lhs
            haskey(idxmap, v) && add_edge!(graph, i, idxmap[v])
        end
        # TODO: custom vars that handles D(x)
        vs = vars(eq.rhs)
        for v in vs
            haskey(idxmap, v) && add_edge!(graph, i, idxmap[v])
        end
    end

    varassoc = Int[(1:dxvar_offset) .+ dxvar_offset; zeros(Int, length(fullvars) - dxvar_offset)] # variable association list
    sys, dxvar_offset, fullvars, varassoc, graph, solvable_graph
end

function reordersys(sys, dxvar_offset, fullvars)
    eqs = equations(sys)
    neweqs = Vector{Equation}(undef, length(eqs))
    eqidxmap = Dict(@view(fullvars[dxvar_offset+1:2dxvar_offset]) .=> (1:dxvar_offset))
    varidxmap = Dict([@view(fullvars[1:dxvar_offset]); @view(fullvars[2dxvar_offset+1:end])] .=> (1:length(fullvars)-dxvar_offset))
    algidx = dxvar_offset
    for eq in eqs
        if isdiffeq(eq)
            neweqs[eqidxmap[eq.lhs]] = eq
        else
            neweqs[algidx+=1] = eq
        end
    end
    sts = states(sys)
    @set! sys.eqs = neweqs
    @set! sys.states = sts[map(s->varidxmap[s], sts)]
end
