const ObsSubberT = typeof(SU.IRSubstituter{false}(IRStructure{VartypeT}(), Dict{SymbolicT, SymbolicT}()))

struct IRInfo
    bound_ps_idxs::Vector{Int}
    bindings_idxs::Vector{Int}
    obs_vars_idxs::Vector{Int}
    obs_idxs::Vector{Int}
    dervars_idxs::Vector{Int}
    dervals_idxs::Vector{Int}
    eqs_idxs::Vector{Int}
    obs_subber::ObsSubberT
end

function should_invalidate_mutable_cache_entry(::Type{IRInfo}, patch::NamedTuple)
    return haskey(patch, :eqs) || haskey(patch, :bindings) || haskey(patch, :parameter_bindings_graph) ||
        haskey(patch, :observed) || haskey(patch, :schedule) || haskey(patch, :irstructure)
end

function get_ir_info(sys::System)
    info = check_mutable_cache(sys, IRInfo, IRInfo, nothing)
    if info isa IRInfo
        return info
    end
    ir = get_irstructure(sys)
    sizehint!(ir, 50 * length(observed(sys)) + 10 * length(equations(sys)))

    subrules = Dict{SymbolicT, SymbolicT}()
    subber = SU.IRSubstituter{false}(ir, subrules)

    pbg = get_parameter_bindings_graph(sys)
    if pbg === nothing
        pbg = ParameterBindingsGraph(sys)
    end
    bound_ps = pbg.bound_ps
    binds = bindings(sys)
    bound_ps_idxs = Int[]
    bindings_idxs = Int[]
    sizehint!(bound_ps_idxs, length(bound_ps))
    sizehint!(bindings_idxs, length(bound_ps))
    for par in bound_ps
        push!(bound_ps_idxs, populate_ir!(ir, par))
        newsym = subber(binds[par])
        push!(bindings_idxs, ir[newsym])
        subrules[par] = newsym
    end

    obs = observed(sys)
    obs_vars_idxs = Int[]
    obs_idxs = Int[]
    sizehint!(obs_vars_idxs, length(obs))
    sizehint!(obs_idxs, length(obs))
    for eq in obs
        push!(obs_vars_idxs, populate_ir!(ir, eq.lhs))
        newsym = subber(eq.rhs)
        push!(obs_idxs, ir[newsym])
        subrules[eq.lhs] = newsym
    end

    eqs = equations(sys)
    eqs_idxs = Int[]
    sizehint!(eqs_idxs, length(eqs))
    for eq in eqs
        newsym = subber(eq.rhs)
        push!(eqs_idxs, ir[newsym])
    end

    dervars_idxs = Int[]
    dervals_idxs = Int[]
    sizehint!(dervars_idxs, length(eqs))
    sizehint!(dervals_idxs, length(eqs))
    if isscheduled(sys)
        sched::Schedule = get_schedule(sys)
        for (k, v) in sched.dummy_sub
            ttk = default_toterm(k)
            isequal(ttk, v) && continue

            push!(dervars_idxs, populate_ir!(ir, ttk))
            newsym = subber(v)
            push!(dervals_idxs, ir[newsym])
            subrules[ttk] = newsym
        end
    elseif is_time_dependent(sys)
        for (i, eq) in enumerate(eqs)
            isdiffeq(eq) || continue
            ttk = default_toterm(eq.lhs)
            isequal(ttk, eq.rhs) && continue

            push!(dervars_idxs, populate_ir!(ir, ttk))
            push!(dervals_idxs, eqs_idxs[i])
            subrules[ttk] = ir[eqs_idxs[i]]
        end
    end

    irinfo = IRInfo(
        bound_ps_idxs, bindings_idxs, obs_vars_idxs, obs_idxs, dervars_idxs,
        dervals_idxs, eqs_idxs, subber
    )
    store_to_mutable_cache!(sys, IRInfo, irinfo)

    return irinfo
end
