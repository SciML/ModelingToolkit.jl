function mtkcompile(
        sys::System; additional_passes = (), simplify = false, split = true,
        kwargs...)
    newsys = __mtkcompile(sys; simplify,
        kwargs...)
    @set! newsys.parent = complete(sys; split = false, flatten = false)
    newsys = complete(newsys; split)
    return newsys
end
function __mtkcompile(sys::AbstractSystem; simplify = false,
        sort_eqs = true,
        kwargs...)
    if isempty(equations(sys)) && !is_time_dependent(sys) && !_iszero(cost(sys))
        return simplify_optimization_system(sys; kwargs..., sort_eqs, simplify)::System
    end
    sys, statemachines = extract_top_level_statemachines(sys)
    sys = flatten(sys)
    state = TearingState(sys)
    append!(state.statemachines, statemachines)
    @unpack structure, fullvars = state
    @unpack graph, var_to_diff, var_types = structure
    return state.sys
end
function simplify_optimization_system(sys::System; split = true, kwargs...)
    sys = flatten(sys)
    econs = constraints(sys)
    dvs = SymbolicT[]
    for var in unknowns(sys)
        sh = SU.shape(var)::SU.ShapeVecT
        if isempty(sh)
            push!(dvs, var)
        else
            append!(dvs, vec(collect(var)::Array{SymbolicT})::Vector{SymbolicT})
        end
    end
    nlsys = System(econs, dvs, parameters(sys); name = :___tmp_nlsystem)
    snlsys = mtkcompile(nlsys; kwargs..., fully_determined = false)::System
    return sys
end
