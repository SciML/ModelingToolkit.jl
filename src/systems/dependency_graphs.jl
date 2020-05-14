# variables equations depend on as a vector of vectors of variables
# each system type should define extract_variables! for a single equation
function equation_dependencies(sys::AbstractSystem; variables=states(sys))
    eqs  = equations(sys)
    deps = Set{Variable}()
    depeqs_to_vars = Vector{Vector{Variable}}(undef,length(eqs))

    for (i,eq) in enumerate(eqs)      
        depeqs_to_vars[i] = collect(get_variables!(deps, eq, variables))
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

    deps = Set{Variable}()
    badjlist = Vector{Vector{Int}}(undef, length(eqs))    
    for (eidx,eq) in enumerate(eqs)
        modified_states!(deps, eq, variables)
        badjlist[eidx] = sort!([vtois[var] for var in deps])
        empty!(deps)
    end

    fadjlist = [Vector{Int}() for i = 1:length(variables)]
    ne = 0
    for (eqidx,vidxs) in enumerate(badjlist)
        println(vidxs)
        foreach(vidx -> push!(fadjlist[vidx], eqidx), vidxs)
        ne += length(vidxs)
    end

    BipartiteGraph(ne, fadjlist, badjlist)
end