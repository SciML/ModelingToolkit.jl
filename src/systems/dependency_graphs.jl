# each system type should define extract_variables! for a single equation
function extract_variables(eqs, vars)
    deps = Set{Variable}()
    depeqs_to_vars = Vector{Vector{Variable}}(undef,length(eqs))

    for (i,eq) in enumerate(eqs)      
        depeqs_to_vars[i] = collect(extract_variables!(deps, eq, vars))
        empty!(deps)
    end

    depeqs_to_vars
end

# variables equations depend on as a vector of vectors of variables
equation_dependencies(sys::AbstractSystem; variables=states(sys)) = extract_variables(equations(sys), variables)

# modeled on LightGraphs SimpleGraph
mutable struct BiPartiteGraph{T <: Integer}
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
        ne += 1
    end

    BiPartiteGraph(ne, fadjlist, badjlist)
end

# could be made to directly generate graph and save memory
function asgraph(sys::AbstractSystem; variables=nothing, variablestoids=nothing)
    vs   = isnothing(variables) ? states(sys) : variables    
    eqdeps = extract_variables(equations(sys), vs)
    vtois  = isnothing(variablestoids) ? Dict(convert(Variable, v) => i for (i,v) in enumerate(vs)) : variablestoids
    asgraph(eqdeps, vtois)
end

# map each variable to the eqs depending on it 
function variables_to_depeqs(sys::AbstractSystem; equationdeps = nothing, statestoids = nothing)
    sts   = states(sys)
    stoi  = isnothing(statestoids) ? Dict(convert(Variable,state) => i for (i,state) in enumerate(sts)) : statestoids

    # map from eqs to states they depend on
    eqdeps = isnothing(equationdeps) ? equation_dependencies(sys) : equationdeps

    # reverse map and switch to integer indices of states
    dg = [Vector{Int}() for i = 1:length(sts)]
    for (k,dep) in enumerate(eqdeps)
        for state in dep
            push!(dg[stoi[state]],k)
        end
    end
    foreach(dep -> sort!(dep), dg)

    dg
end