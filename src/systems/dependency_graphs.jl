
# each system type should define extract_variables! for a single equation
function extract_variables(eqs, vars)
    deps = Set{Variable}()
    depeqs_to_vars = Vector{Vector{Variable}}(undef,length(eqs))

    for (i,eq) in enumerate(eqs)        
        depeqs_to_vars[i] = collect(extract_variables!(deps, eq, vars))
        empty!(deps)
    end
end

equation_dependencies(sys) = extract_variables(equations(sys), states(sys))

# map each variable to the eqs depending on it 
function variables_to_depeqs(sys::AbstractSystem; equationdeps = nothing, statestoids = nothing)
    sts   = states(sys)
    stoi  = isnothing(statestoids) ? Dict(convert(Variable,state) => i for (i,state) in enumerate(sts)) : statestoids

    # map from eqs to states they depend on
    eqdeps = isnothing(equationdeps) ? equation_dependencies(equations(sys), sts) : equationdeps

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