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
odesys = convert(ODESystem, rs)

# dependency of each ODE on state variables
equation_dependencies(odesys)

# dependency of each ODE on parameters
equation_dependencies(odesys, variables=parameters(odesys))

# Jumps
jumpsys = convert(JumpSystem, rs)

# dependency of each jump rate function on state variables
equation_dependencies(jumpsys)

# dependency of each jump rate function on parameters
equation_dependencies(jumpsys, variables=parameters(jumpsys))
```
"""
function equation_dependencies(sys::AbstractSystem; variables=states(sys))
    eqs  = equations(sys)
    deps = Set()
    depeqs_to_vars = Vector{Vector}(undef,length(eqs))

    for (i,eq) in enumerate(eqs)
        get_variables!(deps, eq, variables)
        depeqs_to_vars[i] = [value(v) for v in deps]
        empty!(deps)
    end

    depeqs_to_vars
end

"""
```julia
asgraph(eqdeps, vtois)
```

Convert a collection of equation dependencies, for example as returned by
`equation_dependencies`, to a [`BipartiteGraph`](@ref).

Notes:
- `vtois` should provide a `Dict` like mapping from each `Variable` dependency in
  `eqdeps` to the integer idx of the variable to use in the graph.

Example:
Continuing the example started in [`equation_dependencies`](@ref)
```julia
digr = asgraph(equation_dependencies(odesys), Dict(s => i for (i,s) in enumerate(states(odesys))))
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

Convert an `AbstractSystem` to a [`BipartiteGraph`](@ref) mapping the index of equations
to the indices of variables they depend on.

Notes:
- Defaults for kwargs creating a mapping from `equations(sys)` to `states(sys)`
  they depend on.
- `variables` should provide the list of variables to use for generating
  the dependency graph.
- `variablestoids` should provide `Dict` like mapping from a `Variable` to its
  `Int` index within `variables`.

Example:
Continuing the example started in [`equation_dependencies`](@ref)
```julia
digr = asgraph(odesys)
```
"""
function asgraph(sys::AbstractSystem; variables=states(sys),
                                      variablestoids=Dict(v => i for (i,v) in enumerate(variables)))
    asgraph(equation_dependencies(sys, variables=variables), variablestoids)
end

"""
```julia
variable_dependencies(sys::AbstractSystem; variables=states(sys), variablestoids=nothing)
```

For each variable determine the equations that modify it and return as a [`BipartiteGraph`](@ref).

Notes:
- Dependencies are returned as a [`BipartiteGraph`](@ref) mapping variable
  indices to the indices of equations that modify them.
- `variables` denotes the list of variables to determine dependencies for.
- `variablestoids` denotes a `Dict` mapping `Variable`s to their `Int` index in `variables`.

Example:
Continuing the example of [`equation_dependencies`](@ref)
```julia
variable_dependencies(odesys)
```
"""
function variable_dependencies(sys::AbstractSystem; variables=states(sys), variablestoids=nothing)
    eqs   = equations(sys)
    vtois = isnothing(variablestoids) ? Dict(v => i for (i,v) in enumerate(variables)) : variablestoids

    deps = Set()
    badjlist = Vector{Vector{Int}}(undef, length(eqs))
    for (eidx,eq) in enumerate(eqs)
        modified_states!(deps, eq, variables)
        badjlist[eidx] = sort!([vtois[var] for var in deps])
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
```julia
asdigraph(g::BipartiteGraph, sys::AbstractSystem; variables = states(sys), equationsfirst = true)
```

Convert a [`BipartiteGraph`](@ref) to a `LightGraph.SimpleDiGraph`.

Notes:
- The resulting `SimpleDiGraph` unifies the two sets of vertices (equations
  and then states in the case it comes from [`asgraph`](@ref)), producing one
  ordered set of integer vertices (`SimpleDiGraph` does not support two distinct
  collections of vertices so they must be merged).
- `variables` gives the variables that `g` is associated with (usually the
  `states` of a system).
- `equationsfirst` (default is `true`) gives whether the [`BipartiteGraph`](@ref)
  gives a mapping from equations to variables they depend on (`true`), as calculated
  by [`asgraph`](@ref), or whether it gives a mapping from variables to the equations
  that modify them, as calculated by [`variable_dependencies`](@ref).

Example:
Continuing the example in [`asgraph`](@ref)
```julia
dg = asdigraph(digr)
```
"""
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

"""
```julia
eqeq_dependencies(eqdeps::BipartiteGraph{T}, vardeps::BipartiteGraph{T}) where {T <: Integer}
```

Calculate a `LightGraph.SimpleDiGraph` that maps each equation to equations they depend on.

Notes:
- The `fadjlist` of the `SimpleDiGraph` maps from an equation to the equations that
  modify variables it depends on.
- The `badjlist` of the `SimpleDiGraph` maps from an equation to equations that
  depend on variables it modifies.

Example:
Continuing the example of `equation_dependencies`
```julia
eqeqdep = eqeq_dependencies(asgraph(odesys), variable_dependencies(odesys))
```
"""
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

"""
```julia
varvar_dependencies(eqdeps::BipartiteGraph{T}, vardeps::BipartiteGraph{T}) where {T <: Integer} = eqeq_dependencies(vardeps, eqdeps)
```

Calculate a `LightGraph.SimpleDiGraph` that maps each variable to variables they depend on.

Notes:
- The `fadjlist` of the `SimpleDiGraph` maps from a variable to the variables that
  depend on it.
- The `badjlist` of the `SimpleDiGraph` maps from a variable to variables on which
  it depends.

Example:
Continuing the example of `equation_dependencies`
```julia
varvardep = varvar_dependencies(asgraph(odesys), variable_dependencies(odesys))
```
"""
varvar_dependencies(eqdeps::BipartiteGraph{T}, vardeps::BipartiteGraph{T}) where {T <: Integer} = eqeq_dependencies(vardeps, eqdeps)
