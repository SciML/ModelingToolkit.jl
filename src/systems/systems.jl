function System(eqs::AbstractVector{<:Equation}, iv, args...; name = nothing,
        kw...)
    ODESystem(eqs, iv, args...; name, kw..., checks = false)
end

const REPEATED_SIMPLIFICATION_MESSAGE = "Structural simplification cannot be applied to a completed system. Double simplification is not allowed."

struct RepeatedStructuralSimplificationError <: Exception end

function Base.showerror(io::IO, e::RepeatedStructuralSimplificationError)
    print(io, REPEATED_SIMPLIFICATION_MESSAGE)
end

"""
$(SIGNATURES)

Structurally simplify algebraic equations in a system and compute the
topological sort of the observed equations in `sys`.

### Optional Arguments:
+ optional argument `io` may take a tuple `(inputs, outputs)`. This will convert all `inputs` to parameters and allow them to be unconnected, i.e., simplification will allow models where `n_unknowns = n_equations - n_inputs`.

### Optional Keyword Arguments:
+ When `simplify=true`, the `simplify` function will be applied during the tearing process.
+ `allow_symbolic=false`, `allow_parameter=true`, and `conservative=false` limit the coefficient types during tearing. In particular, `conservative=true` limits tearing to only solve for trivial linear systems where the coefficient has the absolute value of ``1``.
+ `fully_determined=true` controls whether or not an error will be thrown if the number of equations don't match the number of inputs, outputs, and equations.
"""
function structural_simplify(
        sys::AbstractSystem, io = nothing; additional_passes = [], simplify = false, split = true,
        allow_symbolic = false, allow_parameter = true, conservative = false, fully_determined = true,
        kwargs...)
    isscheduled(sys) && throw(RepeatedStructuralSimplificationError())
    newsys′ = __structural_simplify(sys, io; simplify,
        allow_symbolic, allow_parameter, conservative, fully_determined,
        kwargs...)
    if newsys′ isa Tuple
        @assert length(newsys′) == 2
        newsys = newsys′[1]
    else
        newsys = newsys′
    end
    if newsys isa DiscreteSystem &&
       any(eq -> symbolic_type(eq.lhs) == NotSymbolic(), equations(newsys))
        error("""
            Encountered algebraic equations when simplifying discrete system. Please construct \
            an ImplicitDiscreteSystem instead.
        """)
    end
    for pass in additional_passes
        newsys = pass(newsys)
    end
    if newsys isa ODESystem || has_parent(newsys)
        @set! newsys.parent = complete(sys; split = false, flatten = false)
    end
    newsys = complete(newsys; split)
    if newsys′ isa Tuple
        idxs = [parameter_index(newsys, i) for i in io[1]]
        return newsys, idxs
    else
        return newsys
    end
end

function __structural_simplify(sys::JumpSystem, args...; kwargs...)
    return sys
end

function __structural_simplify(sys::SDESystem, args...; kwargs...)
    return __structural_simplify(ODESystem(sys), args...; kwargs...)
end

function __structural_simplify(sys::AbstractSystem, io = nothing; simplify = false,
        kwargs...)
    sys = expand_connections(sys)
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
        @set! sys.unknowns = [v
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
        # Fix for https://github.com/SciML/ModelingToolkit.jl/issues/2490
        if sorted_g_rows isa AbstractMatrix && size(sorted_g_rows, 2) == 1
            # If there's only one brownian variable referenced across all the equations,
            # we get a Nx1 matrix of noise equations, which is a special case known as scalar noise
            noise_eqs = sorted_g_rows[:, 1]
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

        noise_eqs = StructuralTransformations.tearing_substitute_expr(ode_sys, noise_eqs)
        ssys = SDESystem(Vector{Equation}(full_equations(ode_sys)), noise_eqs,
            get_iv(ode_sys), unknowns(ode_sys), parameters(ode_sys);
            name = nameof(ode_sys), is_scalar_noise, observed = observed(ode_sys), defaults = defaults(sys),
            parameter_dependencies = parameter_dependencies(sys), assertions = assertions(sys),
            guesses = guesses(sys), initialization_eqs = initialization_equations(sys))
        @set! ssys.tearing_state = get_tearing_state(ode_sys)
        return ssys
    end
end

"""
    $(TYPEDSIGNATURES)

Given a system that has been simplified via `structural_simplify`, return a `Dict` mapping
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
        throw(ArgumentError("`map_variables_to_equations` requires a simplified system. Call `structural_simplify` on the system before calling this function."))
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
