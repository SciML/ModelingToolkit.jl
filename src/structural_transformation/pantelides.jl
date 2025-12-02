###
### Reassemble: structural information -> system
###

const NOTHING_EQ = nothing ~ nothing

function pantelides_reassemble(state::TearingState, var_eq_matching)
    fullvars = state.fullvars
    @unpack var_to_diff, eq_to_diff = state.structure
    sys = state.sys
    # Step 1: write derivative equations
    in_eqs = equations(sys)
    out_eqs = Vector{Equation}(undef, nv(eq_to_diff))
    fill!(out_eqs, NOTHING_EQ)
    out_eqs[1:length(in_eqs)] .= in_eqs

    out_vars = Vector{SymbolicT}(undef, nv(var_to_diff))
    fill!(out_vars, ModelingToolkit.COMMON_NOTHING)
    out_vars[1:length(fullvars)] .= fullvars

    iv = get_iv(sys)
    D = Differential(iv)

    for (varidx, diff) in edges(var_to_diff)
        # fullvars[diff] = D(fullvars[var])
        vi = out_vars[varidx]
        @assert vi!==ModelingToolkit.COMMON_NOTHING "Something went wrong on reconstructing unknowns from variable association list"
        # `fullvars[i]` needs to be not a `D(...)`, because we want the DAE to be
        # first-order.
        if isdifferential(vi)
            vi = out_vars[varidx] = diff2term_with_unit(vi, iv)
        end
        out_vars[diff] = D(vi)
    end

    d_dict = Dict{SymbolicT, Int}(zip(fullvars, 1:length(fullvars)))
    for (eqidx, diff) in edges(eq_to_diff)
        # LHS variable is looked up from var_to_diff
        # the var_to_diff[i]-th variable is the differentiated version of var at i
        eq = out_eqs[eqidx]
        lhs = if SU.isconst(eq.lhs)
            Symbolics.COMMON_ZERO
        elseif isdiffeq(eq)
            # look up the variable that represents D(lhs)
            lhsarg1 = arguments(eq.lhs)[1]
            @assert !(lhsarg1 isa Differential) "The equation $eq is not first order"
            i = get(d_dict, lhsarg1, nothing)
            if i === nothing
                D(eq.lhs)
            else
                # remove clashing equations
                lhs = ModelingToolkit.COMMON_NOTHING
            end
        else
            D(eq.lhs)
        end
        rhs = ModelingToolkit.expand_derivatives(D(eq.rhs))
        rhs = substitute(rhs, state.param_derivative_map)
        substitution_dict = Dict(x.lhs => x.rhs
        for x in out_eqs if x !== NOTHING_EQ && !SU.isconst(x.lhs))
        sub_rhs = substitute(rhs, substitution_dict)
        out_eqs[diff] = lhs ~ sub_rhs
    end

    final_vars = unique(filter(x -> !(operation(x) isa Differential), fullvars))
    final_eqs = map(identity,
        filter(x -> x.lhs !== ModelingToolkit.COMMON_NOTHING,
            out_eqs[sort(filter(x -> x !== unassigned, var_eq_matching))]))

    @set! sys.eqs = final_eqs
    @set! sys.unknowns = final_vars
    return sys
end

"""
    dae_index_lowering(sys::System; kwargs...) -> System

Perform the Pantelides algorithm to transform a higher index DAE to an index 1
DAE. `kwargs` are forwarded to [`pantelides!`](@ref). End users are encouraged to call [`mtkcompile`](@ref)
instead, which calls this function internally.
"""
function dae_index_lowering(sys::System; kwargs...)
    state = TearingState(sys)
    var_eq_matching = StateSelection.pantelides!(state; finalize = false, kwargs...)
    return invalidate_cache!(pantelides_reassemble(state, var_eq_matching))
end
