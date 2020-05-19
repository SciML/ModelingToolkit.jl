"""
```julia
equation_dependencies(sys::AbstractSystem; variables=states(sys))
```

Given an `AbstractSystem` calculate for each equation the variables it depends on. 

Notes:
- Variables that are not in `variables` are filtered out.
- `get_variables!` is used to determine the variables within a given equation. 
- returns a `Vector{Vector{Variable}}()` mapping the index of an equation to the `variables` it depends on.

Example:
```julia
using ModelingToolkit
@parameters β γ κ η t
@variables S(t) I(t) R(t)

# use a reaction system to easily generate ODE and jump systems
rxs = [Reaction(β, [S,I], [I], [1,1], [2]),
       Reaction(γ, [I], [R]),
       Reaction(κ+η, [R], [S])]
rs = ReactionSystem(rxs, t, [S,I,R], [β,γ,κ,η])

# ODEs:
os = convert(ODESystem, rs)

# dependency of each ODE on state variables
equation_dependencies(os)    

# dependency of each ODE on parameters
equation_dependencies(os, variables=parameters(os))

# Jumps
js = convert(JumpSystem, rs)

# dependency of each jump rate function on state variables
equation_dependencies(js)    

# dependency of each jump rate function on parameters
equation_dependencies(js, variables=parameters(js))    
```
"""
function equation_dependencies(sys::AbstractSystem; variables=states(sys))
    eqs  = equations(sys)
    deps = Set{Operation}()
    depeqs_to_vars = Vector{Vector{Variable}}(undef,length(eqs))

    for (i,eq) in enumerate(eqs)      
        get_variables!(deps, eq, variables)
        depeqs_to_vars[i] = [convert(Variable,v) for v in deps]
        empty!(deps)
    end

    depeqs_to_vars
end

"""
$(TYPEDEF)

A bipartite graph representation between two, possibly distinct, sets of vertices 
(source and dependencies). Maps source vertices, labelled `1:N₁`, to vertices 
on which they depend (labelled `1:N₂`).

# Fields
$(FIELDS)

# Example
```julia
using ModelingToolkit

ne = 4
srcverts = 1:4
depverts = 1:2

# six source vertices
fadjlist = [[1],[1],[2],[2],[1],[1,2]]

# two vertices they depend on 
badjlist = [[1,2,5,6],[3,4,6]]

bg = BipartiteGraph(7, fadjlist, badjlist)
```
"""
mutable struct BipartiteGraph{T <: Integer}
    """Number of edges from source vertices to vertices they depend on."""
    ne::Int
    """Forward adjacency list mapping index of source vertices to the vertices they depend on."""
    fadjlist::Vector{Vector{T}}  # fadjlist[src] = [dest1,dest2,...]
    """Backwrad adjacency list mapping index of vertices that are dependencies to the source vertices that depend on them."""
    badjlist::Vector{Vector{T}}  # badjlist[dst] = [src1,src2,...]
end

function Base.isequal(bg1::BipartiteGraph{T}, bg2::BipartiteGraph{T}) where {T<:Integer} 
    iseq = (bg1.ne == bg2.ne)
    iseq &= (bg1.fadjlist == bg2.fadjlist) 
    iseq &= (bg1.badjlist == bg2.badjlist) 
    iseq
end

"""
```julia
asgraph(eqdeps, vtois)    
```

Convert a collection of equation dependencies, for example as returned by 
`equation_dependencies`, to a `BipartiteGraph`.

Notes:
- `vtois` should provide `Dict` like mapping from variable dependency in `eqdeps`
  to the integer idx of the variable to use in the graph.

Example:
Continuing the example started in [`equation_dependencies`](@ref)
```julia
digr = asgraph(equation_dependencies(os), Dict(s => i for (i,s) in enumerate(states(os))))
```
"""
function asgraph(eqdeps, vtois)    
    fadjlist = Vector{Vector{Int}}(undef, length(eqdeps))
    for (i,dep) in enumerate(eqdeps)
        fadjlist[i] = sort!([vtois[var] for var in dep])
    end

    badjlist = [Vector{Int}() for i = 1:length(vtois)]
    ne = 0
    for (eqidx,vidxs) in enumerate(fadjlist)
        foreach(vidx -> push!(badjlist[vidx], eqidx), vidxs)
        ne += length(vidxs)
    end

    BipartiteGraph(ne, fadjlist, badjlist)
end


# could be made to directly generate graph and save memory
"""
```julia
asgraph(sys::AbstractSystem; variables=states(sys), 
                                      variablestoids=Dict(convert(Variable, v) => i for (i,v) in enumerate(variables)))
```

Convert an `AbstractSystem` to a `BipartiteGraph` mapping equations
to variables they depend on.

Notes:
- Defaults for kwargs creating a mapping from `equations(sys)` to `states(sys)`
  they depend on.
- `variables` should provide the list of variables to use for generating 
  the dependency graph.
- `variablestoids` should provide `Dict` like mapping from a variable to its 
  integer index within `variables`.

Example:
Continuing the example started in [`equation_dependencies`](@ref)
```julia
digr = asgraph(os)
```
"""
function asgraph(sys::AbstractSystem; variables=states(sys), 
                                      variablestoids=Dict(convert(Variable, v) => i for (i,v) in enumerate(variables)))
    asgraph(equation_dependencies(sys, variables=variables), variablestoids)
end
# function asgraph(sys::AbstractSystem; variables=nothing, variablestoids=nothing)
#     vs     = isnothing(variables) ? states(sys) : variables    
#     eqdeps = equation_dependencies(sys, variables=vs)
#     vtois  = isnothing(variablestoids) ? Dict(convert(Variable, v) => i for (i,v) in enumerate(vs)) : variablestoids
#     asgraph(eqdeps, vtois)
# end

# for each variable determine the equations that modify it
function variable_dependencies(sys::AbstractSystem; variables=states(sys), variablestoids=nothing)
    eqs   = equations(sys)
    vtois = isnothing(variablestoids) ? Dict(convert(Variable, v) => i for (i,v) in enumerate(variables)) : variablestoids

    deps = Set{Operation}()
    badjlist = Vector{Vector{Int}}(undef, length(eqs))    
    for (eidx,eq) in enumerate(eqs)
        modified_states!(deps, eq, variables)
        badjlist[eidx] = sort!([vtois[convert(Variable,var)] for var in deps])
        empty!(deps)
    end

    fadjlist = [Vector{Int}() for i = 1:length(variables)]
    ne = 0
    for (eqidx,vidxs) in enumerate(badjlist)
        foreach(vidx -> push!(fadjlist[vidx], eqidx), vidxs)
        ne += length(vidxs)
    end

    BipartiteGraph(ne, fadjlist, badjlist)
end

"""
- The resulting `SimpleDiGraph` unifies the two sets of vertices (equations 
  and then states in the case `eqdeps` comes from `equation_dependencies`), producing
  one ordered set of integer vertices (as `SimpleDiGraph` does not support two distinct
  collections of nodes.
"""
# convert BipartiteGraph to LightGraph.SimpleDiGraph
function asdigraph(g::BipartiteGraph, sys::AbstractSystem; variables = states(sys), equationsfirst = true)
    neqs     = length(equations(sys))
    nvars    = length(variables)
    fadjlist = deepcopy(g.fadjlist)
    badjlist = deepcopy(g.badjlist)

    # offset is for determining indices for the second set of vertices
    offset = equationsfirst ? neqs : nvars
    for i = 1:offset
        fadjlist[i] .+= offset
    end

    # add empty rows for vertices without connections
    append!(fadjlist, [Vector{Int}() for i=1:(equationsfirst ? nvars : neqs)])
    prepend!(badjlist, [Vector{Int}() for i=1:(equationsfirst ? neqs : nvars)])

    SimpleDiGraph(g.ne, fadjlist, badjlist)
end

# maps the i'th eq to equations that depend on it
function eqeq_dependencies(eqdeps::BipartiteGraph{T}, vardeps::BipartiteGraph{T}) where {T <: Integer}
    g = SimpleDiGraph{T}(length(eqdeps.fadjlist))
    
    for (eqidx,sidxs) in enumerate(vardeps.badjlist)
        # states modified by eqidx
        for sidx in sidxs
            # equations depending on sidx
            foreach(v -> add_edge!(g, eqidx, v), eqdeps.badjlist[sidx])
        end
    end

    g
end

# maps the i'th variable to variables that depend on it
varvar_dependencies(eqdeps::BipartiteGraph{T}, vardeps::BipartiteGraph{T}) where {T <: Integer} = eqeq_dependencies(vardeps, eqdeps)
