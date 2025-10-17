function mtkcompile(
        sys::System; additional_passes = (), simplify = false, split = true,
        kwargs...)
    newsys = __mtkcompile(sys; simplify,
        kwargs...)
    @set! newsys.parent = complete(sys; split = false, flatten = false)
    newsys = complete(newsys; split)
end
function __mtkcompile(sys::AbstractSystem; simplify = false,
        kwargs...)
    if isempty(equations(sys)) && !is_time_dependent(sys) && !_iszero(cost(sys))
        return simplify_optimization_system(sys; kwargs..., sort_eqs, simplify)::System
    end
    state = TearingState(sys)
    return state.sys
end
function simplify_optimization_system(sys::System; split = true, kwargs...)
    for var in unknowns(sys)
        if isempty(sh)
        end
    end
    snlsys = mtkcompile(nlsys; kwargs..., fully_determined = false)::System
end
