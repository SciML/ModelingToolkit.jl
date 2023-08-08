function System(eqs::AbstractVector{<:Equation}, iv = nothing, args...; name = nothing,
    kw...)
    ODESystem(eqs, iv, args...; name, checks = false)
end

"""
$(SIGNATURES)

Structurally simplify algebraic equations in a system and compute the
topological sort of the observed equations. When `simplify=true`, the `simplify`
function will be applied during the tearing process. It also takes kwargs
`allow_symbolic=false` and `allow_parameter=true` which limits the coefficient
types during tearing.

The optional argument `io` may take a tuple `(inputs, outputs)`.
This will convert all `inputs` to parameters and allow them to be unconnected, i.e.,
simplification will allow models where `n_states = n_equations - n_inputs`.
"""
function structural_simplify(sys::AbstractSystem, io = nothing; simplify = false,
    kwargs...)
    sys = expand_connections(sys)
    sys isa DiscreteSystem && return sys
    state = TearingState(sys)

    @unpack structure, fullvars = state
    @unpack graph, var_to_diff, var_types = structure
    eqs = equations(state)
    brown_vars = Int[]
    new_idxs = zeros(Int, length(var_types))
    idx = 0
    for (i, vt) in enumerate(var_types)
        if vt === BROWNIAN
            push!(brown_vars, i)
        else
            new_idxs[i] = (idx += 1)
        end
    end
    if isempty(brown_vars)
        return structural_simplify!(state, io; simplify, kwargs...)
    else
        Is = Int[]
        Js = Int[]
        vals = Num[]
        new_eqs = copy(eqs)
        dvar2eq = Dict{Any, Int}()
        for (v, dv) in enumerate(var_to_diff)
            dv === nothing && continue
            deqs = 𝑑neighbors(graph, dv)
            if length(deqs) != 1
                error("$(eqs[deqs]) is not handled.")
            end
            dvar2eq[fullvars[dv]] = only(deqs)
        end
        for (j, bj) in enumerate(brown_vars), i in 𝑑neighbors(graph, bj)
            push!(Is, i)
            push!(Js, j)
            eq = new_eqs[i]
            brown = fullvars[bj]
            (coeff, residual, islinear) = Symbolics.linear_expansion(eq, brown)
            islinear || error("$brown isn't linear in $eq")
            new_eqs[i] = 0 ~ residual
            push!(vals, coeff)
        end
        g = Matrix(sparse(Is, Js, vals))
        sys = state.sys
        @set! sys.eqs = new_eqs
        @set! sys.states = [v
                            for (i, v) in enumerate(fullvars)
                                if !iszero(new_idxs[i]) &&
                                   invview(var_to_diff)[i] === nothing]
        # TODO: IO is not handled.
        ode_sys = structural_simplify(sys, io; simplify, kwargs...)
        eqs = equations(ode_sys)
        sorted_g_rows = zeros(Num, length(eqs), size(g, 2))
        for (i, eq) in enumerate(eqs)
            dvar = eq.lhs
            # differential equations always precede algebraic equations
            _iszero(dvar) && break
            g_row = get(dvar2eq, dvar, 0)
            iszero(g_row) && error("$dvar isn't handled.")
            g_row > size(g, 1) && continue
            @views copyto!(sorted_g_rows[i, :], g[g_row, :])
        end

        return SDESystem(full_equations(ode_sys), sorted_g_rows,
            get_iv(ode_sys), states(ode_sys), parameters(ode_sys);
            name = nameof(ode_sys))
    end
end
