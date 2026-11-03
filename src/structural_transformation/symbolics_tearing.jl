"""
    $TYPEDEF

Tearing algorithm for stationary (no independent variable) systems. Behaves like the default
`DummyDerivativeTearing`, except that every strongly connected component with more than
`max_linear_scc_tearing` equations, all of which are linear in the component's own unknowns, is
left untorn. Tearing such a block (for example the interior of a discretized Laplace
equation) selects a few iteration variables and expresses the rest by a recursion whose
conditioning degrades exponentially with the block size, whereas the untorn block is a linear
system the nonlinear solver handles directly.
"""
struct LinearBlockTearing{S <: TearingState} <: StateSelection.TearingAlgorithm
    state::S
    max_linear_scc_tearing::Int
end

function (alg::LinearBlockTearing)(structure, dummy_derivatives = ())
    (; state, max_linear_scc_tearing) = alg
    mm = state.mm
    big_linear_vars, big_linear_eqs = large_linear_sccs(state, structure, max_linear_scc_tearing)
    if !isempty(big_linear_vars) && mm !== nothing
        # Without its rows in `mm`, the tearing algorithm does not triangularize these blocks.
        keep = findall(!in(big_linear_eqs), mm.nzrows)
        mm = typeof(mm)(
            mm.nparentrows, mm.ncols, mm.nzrows[keep], mm.row_cols[keep], mm.row_vals[keep]
        )
    end
    tearing_result, extra = StateSelection.DummyDerivativeTearing()(structure, dummy_derivatives; mm)
    # A block with no matched variable is solved all at once by the nonlinear solver.
    for v in big_linear_vars
        tearing_result.var_eq_matching[v] = unassigned
    end
    return tearing_result, extra
end

"""
    $(TYPEDSIGNATURES)

Return the variables and equations of the strongly connected components of more than
`max_size` equations that are all linear in the component's variables. The components are
those the tearing algorithm sees, namely of a maximal matching of the structural graph.
"""
function large_linear_sccs(state::TearingState, structure, max_size::Int)
    (; graph) = structure
    big_linear_vars = BitSet()
    big_linear_eqs = BitSet()
    matching = complete(maximal_matching(graph), ndsts(graph))
    eqs = equations(state)
    for scc in StateSelection.find_var_sccs(graph, matching)
        length(scc) > max_size || continue
        scc_vars = BitSet(scc)
        all(v -> matching[v] isa Int && is_linear_in(state, graph, matching[v], scc_vars), scc) || continue
        union!(big_linear_vars, scc)
        union!(big_linear_eqs, matching[v] for v in scc)
    end
    return big_linear_vars, big_linear_eqs
end

"""
    $(TYPEDSIGNATURES)

Check whether equation `ieq` is linear in the variables `scc_vars` it is incident on, with
coefficients that do not depend on any of them. This follows the check for inlining linear
SCCs in reassembly: incidence comes from `graph`, and each variable is expanded with the
cached `get_linear_expander_for!`, continuing with the remainder.
"""
function is_linear_in(state::TearingState, graph, ieq::Int, scc_vars::BitSet)
    eq = equations(state)[ieq]
    expr = SU._iszero(eq.lhs) ? eq.rhs : eq.rhs - eq.lhs
    vars = [v for v in 𝑠neighbors(graph, ieq) if v in scc_vars]
    for (i, v) in enumerate(vars)
        coeff, expr, islinear = MTKBase.get_linear_expander_for!(state.sys, state.fullvars[v], true)(expr)
        islinear || return false
        SU.isconst(coeff) && continue
        # A coefficient containing another unknown of the block makes the equation nonlinear.
        for w in @view vars[(i + 1):end]
            dcoeff, _, coeff_linear = MTKBase.get_linear_expander_for!(state.sys, state.fullvars[w], true)(coeff)
            coeff_linear && SU._iszero(dcoeff) || return false
        end
    end
    return true
end

function tearing(
        state::TearingState;
        tearing_alg::StateSelection.TearingAlgorithm = StateSelection.DummyDerivativeTearing(),
        kwargs...
    )
    state.structure.solvable_graph === nothing && StateSelection.find_solvables!(state; kwargs...)
    StateSelection.complete!(state.structure)
    return tearing_alg(state.structure)
end

"""
    tearing(sys)

Tear the nonlinear equations in system. When `simplify=true`, we simplify the
new residual equations after tearing. End users are encouraged to call [`ModelingToolkitBase.mtkcompile`](@ref)
instead, which calls this function internally.
"""
function tearing(
        sys::AbstractSystem, state = TearingState(sys);
        reassemble_alg::ReassembleAlgorithm = DefaultReassembleAlgorithm(),
        fully_determined = true,
        # Not forwarded: StateSelection and MTKTearing do not accept `verbose`.
        verbose::MTKVerbosity = MTKVerbosity(),
        kwargs...
    )
    tearing_result, extras = tearing(state; kwargs...)
    return invalidate_cache!(reassemble_alg(state, tearing_result, state.mm; fully_determined, kwargs...))
end

function safe_isinteger(@nospecialize(x::Number))
    if x isa Int64
        return typemin(Int) <= x <= typemax(Int)
    elseif x isa Int32
        return true
    elseif x isa BigInt
        return typemin(Int) <= x <= typemax(Int)
    elseif x isa Float64
        return isinteger(x) && typemin(Int) <= x <= typemax(Int)
    elseif x isa Float32
        return isinteger(x) && typemin(Int) <= x <= typemax(Int)
    elseif x isa BigFloat
        return isinteger(x) && typemin(Int) <= x <= typemax(Int)
    elseif x isa Rational{Int64}
        return isinteger(x) && typemin(Int) <= x <= typemax(Int)
    elseif x isa Rational{Int32}
        return isinteger(x) && typemin(Int) <= x <= typemax(Int)
    elseif x isa Rational{BigInt}
        return isinteger(x) && typemin(Int) <= x <= typemax(Int)
    else
        return isinteger(x)::Bool && (typemin(Int) <= x)::Bool && (x <= typemax(Int))::Bool
    end
end

"""
    dummy_derivative(sys)

Perform index reduction and use the dummy derivative technique to ensure that
the system is balanced.

For systems without an independent variable, strongly connected components of more than
`max_linear_scc_tearing` (default 16) equations that are all linear in their own unknowns
are not torn, see `LinearBlockTearing`. Pass `max_linear_scc_tearing = typemax(Int)` to
tear them as before.
"""
function dummy_derivative(
        sys, state = TearingState(sys);
        reassemble_alg::ReassembleAlgorithm = DefaultReassembleAlgorithm(),
        fully_determined = true, max_linear_scc_tearing::Int = 16,
        # Not forwarded: StateSelection and MTKTearing do not accept `verbose`.
        verbose::MTKVerbosity = MTKVerbosity(),
        kwargs...
    )
    jac = let state = state
        (eqs, vars) -> begin
            symeqs = equations(state)[eqs]
            _J = Symbolics.jacobian((x -> x.rhs).(symeqs), state.fullvars[vars])
            J = similar(_J, Int)
            for i in eachindex(_J)
                el = _J[i]
                Moshi.Match.@match el begin
                    BSImpl.Const(; val) && if val isa Number end => begin
                        safe_isinteger(val)|| return nothing
                        J[i] = convert(Int, val)::Int
                    end
                    _ => return nothing
                end
            end
            return J
        end
    end
    state_priority = Base.Fix1(getindex, state.structure.state_priorities)
    tearing_result, extras = if get_iv(sys) === nothing
        # `dummy_derivative_graph!(state, jac)` does not forward a tearing algorithm, so
        # replicate its preamble.
        state.structure.solvable_graph === nothing && find_solvables!(state; kwargs...)
        StateSelection.complete!(state.structure)
        var_eq_matching = complete(StateSelection.pantelides!(state; kwargs...))
        StateSelection.dummy_derivative_graph!(
            state.structure, var_eq_matching, jac, state_priority, Val(false);
            tearing_alg = LinearBlockTearing(state, max_linear_scc_tearing),
            mm = StateSelection.get_mm(state)
        )
    else
        StateSelection.dummy_derivative_graph!(state, jac; state_priority, kwargs...)
    end
    return reassemble_alg(state, tearing_result, state.mm; fully_determined, kwargs...)
end
