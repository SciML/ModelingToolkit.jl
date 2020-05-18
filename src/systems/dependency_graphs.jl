# variables equations depend on as a vector of vectors of variables
# each system type should define extract_variables! for a single equation
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

# modeled on LightGraphs SimpleGraph
mutable struct BipartiteGraph{T <: Integer}
    ne::Int
    fadjlist::Vector{Vector{T}}  # fadjlist[src] = [dest1,dest2,...]
    badjlist::Vector{Vector{T}}  # badjlist[dst] = [src1,src2,...]
end

# convert equation-variable dependencies to a bipartite graph
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

function Base.isequal(bg1::BipartiteGraph{T}, bg2::BipartiteGraph{T}) where {T<:Integer} 
    iseq = (bg1.ne == bg2.ne)
    iseq &= (bg1.fadjlist == bg2.fadjlist) 
    iseq &= (bg1.badjlist == bg2.badjlist) 
    iseq
end

# could be made to directly generate graph and save memory
function asgraph(sys::AbstractSystem; variables=nothing, variablestoids=nothing)
    vs     = isnothing(variables) ? states(sys) : variables    
    eqdeps = equation_dependencies(sys, variables=vs)
    vtois  = isnothing(variablestoids) ? Dict(convert(Variable, v) => i for (i,v) in enumerate(vs)) : variablestoids
    asgraph(eqdeps, vtois)
end

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
