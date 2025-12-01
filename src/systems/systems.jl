@doc """
    function mtkcompile(sys::System; kwargs...)

Compile the given system into a form that ModelingToolkit can generate code for. Also
performs a variety of symbolic-numeric enhancements. For ODEs, this includes processes
such as order reduction, index reduction, alias elimination and tearing. A subset of the
unknowns of the system may be eliminated as observables, eliminating the need for the
numerical solver to solve for these variables.

Does not rely on metadata to identify variables/parameters/brownians. Instead, queries
the system for which symbolic quantites belong to which category. Any variables not
present in the equations of the system will be removed in this process.

# Keyword Arguments

+ When `simplify=true`, the `simplify` function will be applied during the tearing process.
+ `allow_symbolic=false`, `allow_parameter=true`, and `conservative=false` limit the coefficient types during tearing. In particular, `conservative=true` limits tearing to only solve for trivial linear systems where the coefficient has the absolute value of ``1``.
+ `fully_determined=true` controls whether or not an error will be thrown if the number of equations don't match the number of inputs, outputs, and equations.
+ `inputs`, `outputs` and `disturbance_inputs` are passed as keyword arguments.` All inputs` get converted to parameters and are allowed to be unconnected, allowing models where `n_unknowns = n_equations - n_inputs`.
+ `sort_eqs=true` controls whether equations are sorted lexicographically before simplification or not.
""" mtkcompile

function MTKBase.__mtkcompile(sys::System;
        inputs::OrderedSet{SymbolicT} = OrderedSet{SymbolicT}(),
        outputs::OrderedSet{SymbolicT} = OrderedSet{SymbolicT}(),
        disturbance_inputs::OrderedSet{SymbolicT} = OrderedSet{SymbolicT}(),
        sort_eqs = true,
        kwargs...)
    sys, statemachines = extract_top_level_statemachines(sys)
    sys = expand_connections(sys)
    state = TearingState(sys; sort_eqs)
    append!(state.statemachines, statemachines)

    @unpack structure, fullvars = state
    @unpack graph, var_to_diff, var_types = structure
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
        return mtkcompile!(
            state; inputs, outputs, disturbance_inputs, kwargs...)
    else
        Is = Int[]
        Js = Int[]
        vals = SymbolicT[]
        make_eqs_zero_equals!(state)
        new_eqs = copy(equations(state))
        dvar2eq = Dict{SymbolicT, Int}()
        eqs = equations(state)
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
            new_eqs[i] = COMMON_ZERO ~ residual
            push!(vals, coeff)
        end
        g = Matrix(sparse(Is, Js, vals))
        sys = state.sys
        @set! sys.eqs = new_eqs
        @set! sys.unknowns = [v
                              for (i, v) in enumerate(fullvars)
                              if !iszero(new_idxs[i]) &&
                                 invview(var_to_diff)[i] === nothing]
        ode_sys = mtkcompile(
            sys; inputs, outputs, disturbance_inputs, kwargs...)
        eqs = equations(ode_sys)
        sorted_g_rows = fill(COMMON_ZERO, length(eqs), size(g, 2))
        for (i, eq) in enumerate(eqs)
            dvar = eq.lhs
            # differential equations always precede algebraic equations
            _iszero(dvar) && break
            g_row = get(dvar2eq, dvar, 0)
            iszero(g_row) && error("$dvar isn't handled.")
            g_row > size(g, 1) && continue
            @views copyto!(sorted_g_rows[i, :], g[g_row, :])
        end
        # Fix for https://github.com/SciML/ModelingToolkit.jl/issues/2490
        if sorted_g_rows isa AbstractMatrix && size(sorted_g_rows, 2) == 1
            # If there's only one brownian variable referenced across all the equations,
            # we get a Nx1 matrix of noise equations, which is a special case known as scalar noise
            noise_eqs = reshape(sorted_g_rows[:, 1], (:, 1))
            is_scalar_noise = true
        elseif __num_isdiag_noise(sorted_g_rows)
            # If each column of the noise matrix has either 0 or 1 non-zero entry, then this is "diagonal noise".
            # In this case, the solver just takes a vector column of equations and it interprets that to
            # mean that each noise process is independent
            noise_eqs = __get_num_diag_noise(sorted_g_rows)
            is_scalar_noise = false
        else
            noise_eqs = sorted_g_rows
            is_scalar_noise = false
        end

        noise_eqs = substitute_observed(ode_sys, noise_eqs)
        ssys = System(Vector{Equation}(full_equations(ode_sys)),
            get_iv(ode_sys), unknowns(ode_sys),
            [parameters(ode_sys); collect(bound_parameters(ode_sys))]; noise_eqs,
            name = nameof(ode_sys), observed = observed(ode_sys), bindings = bindings(sys),
            initial_conditions = initial_conditions(sys),
            assertions = assertions(sys),
            guesses = guesses(sys), initialization_eqs = initialization_equations(sys),
            continuous_events = continuous_events(sys),
            discrete_events = discrete_events(sys),
            gui_metadata = get_gui_metadata(sys))
        return ssys
    end
end

function MTKBase.simplify_sde_system(sys::System; kwargs...)
    __mtkcompile(sys; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Given a system that has been simplified via `mtkcompile`, return a `Dict` mapping
variables of the system to equations that are used to solve for them. This includes
observed variables.

# Keyword Arguments

- `rename_dummy_derivatives`: Whether to rename dummy derivative variable keys into their
  `Differential` forms. For example, this would turn the key `yˍt(t)` into
  `Differential(t)(y(t))`.
"""
function map_variables_to_equations(sys::AbstractSystem; rename_dummy_derivatives = true)
    if !has_tearing_state(sys)
        throw(ArgumentError("$(typeof(sys)) is not supported."))
    end
    ts = get_tearing_state(sys)
    if ts === nothing
        throw(ArgumentError("`map_variables_to_equations` requires a simplified system. Call `mtkcompile` on the system before calling this function."))
    end

    dummy_sub = Dict()
    if rename_dummy_derivatives && has_schedule(sys) && (sc = get_schedule(sys)) !== nothing
        dummy_sub = Dict(v => k for (k, v) in sc.dummy_sub if isequal(default_toterm(k), v))
    end

    mapping = Dict{Union{Num, BasicSymbolic}, Equation}()
    eqs = equations(sys)
    for eq in eqs
        isdifferential(eq.lhs) || continue
        var = arguments(eq.lhs)[1]
        var = get(dummy_sub, var, var)
        mapping[var] = eq
    end

    graph = ts.structure.graph
    algvars = BitSet(findall(
        Base.Fix1(StructuralTransformations.isalgvar, ts.structure), 1:ndsts(graph)))
    algeqs = BitSet(findall(1:nsrcs(graph)) do eq
        all(!Base.Fix1(isdervar, ts.structure), 𝑠neighbors(graph, eq))
    end)
    alge_var_eq_matching = complete(maximal_matching(graph, in(algeqs), in(algvars)))
    for (i, eq) in enumerate(alge_var_eq_matching)
        eq isa Unassigned && continue
        mapping[get(dummy_sub, ts.fullvars[i], ts.fullvars[i])] = eqs[eq]
    end
    for eq in observed(sys)
        mapping[get(dummy_sub, eq.lhs, eq.lhs)] = eq
    end

    return mapping
end

"""
Mark whether an extra pass `p` can support compiling discrete systems.
"""
discrete_compile_pass(p) = false
