struct IONotFoundError <: Exception
    variant::String
    sysname::Symbol
    not_found::Any
end
function Base.showerror(io::IO, err::IONotFoundError)
end
function markio!(state, orig_inputs::Set{SymbolicT},
                 inputs::Vector{SymbolicT}, outputs::Vector{SymbolicT},
                 disturbances::Vector{SymbolicT}; check = true)
    fullvars = get_fullvars(state)
    inputset = Dict{SymbolicT, Bool}()
    for i in inputs
        inputset[i] = false
    end
    outputset = Dict{SymbolicT, Bool}()
    for o in outputs
        outputset[o] = false
    end
    disturbanceset = Dict{SymbolicT, Bool}()
    for d in disturbances
        disturbanceset[d] = false
    end
    for (i, v) in enumerate(fullvars)
        if haskey(inputset, v)
            if v in keys(outputset)
                v = setio(v, true, true)
                outputset[v] = true
            else
                v = setio(v, true, false)
            end
            inputset[v] = true
            fullvars[i] = v
        elseif haskey(outputset, v)
            v = setio(v, false, true)
            outputset[v] = true
            fullvars[i] = v
        else
            if isinput(v)
                push!(orig_inputs, v)
            end
            v = setio(v, false, false)
            fullvars[i] = v
        end
        if haskey(disturbanceset, v)
            v = setio(v, true, false)
            v = setdisturbance(v, true)
            disturbanceset[v] = true
            fullvars[i] = v
        end
    end
    if check
        ikeys = keys(filter(!last, inputset))
        if !isempty(ikeys)
            throw(IONotFoundError("inputs", nameof(state.sys), ikeys))
        end
        dkeys = keys(filter(!last, disturbanceset))
        if !isempty(dkeys)
            throw(IONotFoundError("disturbance inputs", nameof(state.sys), ikeys))
        end
        okeys = keys(filter(!last, outputset))
        if !isempty(okeys)
            throw(IONotFoundError("outputs", nameof(state.sys), okeys))
        end
    end
    state, orig_inputs
end
