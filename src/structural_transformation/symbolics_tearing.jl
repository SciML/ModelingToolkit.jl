function tearing(state::TearingState;
                 tearing_alg::StateSelection.TearingAlgorithm = StateSelection.DummyDerivativeTearing(),
                 kwargs...)
    state.structure.solvable_graph === nothing && StateSelection.find_solvables!(state; kwargs...)
    StateSelection.complete!(state.structure)
    tearing_alg(state.structure)
end

"""
    tearing(sys)

Tear the nonlinear equations in system. When `simplify=true`, we simplify the
new residual equations after tearing. End users are encouraged to call [`mtkcompile`](@ref)
instead, which calls this function internally.
"""
function tearing(sys::AbstractSystem, state = TearingState(sys); mm = nothing,
        reassemble_alg::ReassembleAlgorithm = DefaultReassembleAlgorithm(),
        fully_determined = true, kwargs...)
    tearing_result, extras = tearing(state; kwargs...)
    invalidate_cache!(reassemble_alg(state, tearing_result, mm; fully_determined))
end

"""
    dummy_derivative(sys)

Perform index reduction and use the dummy derivative technique to ensure that
the system is balanced.
"""
function dummy_derivative(sys, state = TearingState(sys);
        reassemble_alg::ReassembleAlgorithm = DefaultReassembleAlgorithm(),
        mm = nothing, fully_determined = true, kwargs...)
    jac = let state = state
        (eqs, vars) -> begin
            symeqs = equations(state)[eqs]
            _J = Symbolics.jacobian((x -> x.rhs).(symeqs), state.fullvars[vars])
            J = similar(_J, Int)
            for i in eachindex(_J)
                el = _J[i]
                Moshi.Match.@match el begin
                    BSImpl.Const(; val) && if val isa Number end => begin
                        isinteger(val)::Bool || return nothing
                        val = Int(val)
                        typemin(Int) <= val <= typemax(Int) || return nothing
                        J[i] = val
                    end
                    _ => return nothing
                end
            end
            return J
        end
    end
    state_priority = let state = state
        var -> begin
            p = 0.0
            var_to_diff = state.structure.var_to_diff
            diff_to_var = invview(var_to_diff)
            while var_to_diff[var] !== nothing
                var = var_to_diff[var]
            end
            while true
                p = max(p, ModelingToolkit.state_priority(state.fullvars[var]))
                (var = diff_to_var[var]) === nothing && break
            end
            p
        end
    end
    tearing_result, extras = StateSelection.dummy_derivative_graph!(
        state, jac; state_priority, kwargs...)
    reassemble_alg(state, tearing_result, mm; fully_determined)
end
