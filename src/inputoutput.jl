using Symbolics: get_variables
""""""
inputs(sys) = collect(get_inputs(sys))
""""""
function outputs(sys)
    return collect(get_outputs(sys))
end
""""""
bound_inputs(sys) = filter(x -> is_bound(sys, x), inputs(sys))
""""""
unbound_inputs(sys) = filter(x -> !is_bound(sys, x), inputs(sys))
""""""
bound_outputs(sys) = filter(x -> is_bound(sys, x), outputs(sys))
""""""
unbound_outputs(sys) = filter(x -> !is_bound(sys, x), outputs(sys))
""""""
function inputs_to_parameters!(state::TransformationState, inputsyms::Vector{SymbolicT})
    check_bound = inputsyms === nothing
    @unpack structure, fullvars, sys = state
    @unpack var_to_diff, graph, solvable_graph = structure
    @assert solvable_graph === nothing
    inputs = BitSet()
    var_reidx = zeros(Int, length(fullvars))
    ninputs = 0
    nvar = 0
    new_parameters = SymbolicT[]
    input_to_parameters = Dict{SymbolicT, SymbolicT}()
    new_fullvars = SymbolicT[]
    for (i, v) in enumerate(fullvars)
        if isinput(v) && !(check_bound && is_bound(sys, v))
            if var_to_diff[i] !== nothing
                error("Input $(fullvars[i]) is differentiated!")
            end
            push!(inputs, i)
            ninputs += 1
            var_reidx[i] = -1
            p = toparam(v)
            push!(new_parameters, p)
            input_to_parameters[v] = p
        else
            nvar += 1
            var_reidx[i] = nvar
            push!(new_fullvars, v)
        end
    end
    if ninputs == 0
        @set! sys.inputs = OrderedSet{SymbolicT}()
        @set! sys.outputs = OrderedSet{SymbolicT}(filter(isoutput, fullvars))
        state.sys = sys
        return state
    end
    nvars = ndsts(graph) - ninputs
    new_graph = BipartiteGraph(nsrcs(graph), nvars, Val(false))
    for ie in 1:nsrcs(graph)
        for iv in 𝑠neighbors(graph, ie)
            iv = var_reidx[iv]
            iv > 0 || continue
            add_edge!(new_graph, ie, iv)
        end
    end
    new_var_to_diff = DiffGraph(nvars, true)
    for (i, v) in enumerate(var_to_diff)
        new_i = var_reidx[i]
        (new_i < 1 || v === nothing) && continue
        new_v = var_reidx[v]
        @assert new_v > 0
        new_var_to_diff[new_i] = new_v
    end
    @set! structure.var_to_diff = complete(new_var_to_diff)
    @set! structure.graph = complete(new_graph)
    @set! sys.eqs = isempty(input_to_parameters) ? equations(sys) :
                    substitute(equations(sys), input_to_parameters)
    @set! sys.unknowns = setdiff(unknowns(sys), keys(input_to_parameters))
    ps = parameters(sys)
    @set! sys.ps = [ps; new_parameters]
    @set! sys.inputs = OrderedSet{SymbolicT}(new_parameters)
    @set! sys.outputs = OrderedSet{SymbolicT}(filter(isoutput, fullvars))
    @set! state.sys = sys
    @set! state.fullvars = Vector{SymbolicT}(new_fullvars)
    @set! state.structure = structure
    return state
end
