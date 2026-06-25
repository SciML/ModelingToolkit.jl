MTKBase.singular_check(ts::TearingState) = StateSelection.singular_check(ts)

function MTKBase.get_initialization_problem_type(
        sys::System, isys::System;
        warn_initialize_determined = true,
        use_scc = true, kwargs...
    )
    neqs = length(equations(isys))
    nunknown = length(unknowns(isys))
    ts = get_tearing_state(isys)::TearingState

    if use_scc
        scc_message = """
        `SCCNonlinearProblem` can only be used for initialization of fully determined \
        systems and hence will not be used here.
        """
    else
        scc_message = ""
    end

    if warn_initialize_determined && neqs > nunknown
        @warn overdetermined_initialization_message(neqs, nunknown, scc_message)
    end
    if warn_initialize_determined && neqs < nunknown
        @warn underdetermined_initialization_message(neqs, nunknown, scc_message)
    end

    unassigned_vars = MTKBase.singular_check(ts)
    return if nunknown > 0 && nunknown <= neqs && calculate_A_b(isys; throw = false) !== nothing
        LinearInitializationProblem
    elseif neqs == nunknown && isempty(unassigned_vars)
        if use_scc && neqs > 0
            if is_split(isys)
                SCCNonlinearProblem
            else
                @warn """
                `SCCNonlinearProblem` can only be used with `split = true` systems. \
                Simplify your `System` with `split = true` or pass `use_scc = false` to \
                disable this warning
                """
                NonlinearProblem
            end
        else
            NonlinearProblem
        end
    else
        NonlinearLeastSquaresProblem
    end
end

"""
    report_initialization_nullspace(prob; rtol = 1e-8, atol = 0.0, var_threshold = 1e-3, verbose = true)

Diagnose rank deficiency of a system's initialization problem by inspecting the null
space of its residual Jacobian, and report which unknowns span that null space.

When the initialization system is underdetermined or otherwise rank deficient, the
initialization solve can converge to different solutions depending on the solver and on
the (often nondeterministic) ordering of the assembled unknowns/equations. This is a
common source of intermittent initialization failures: some orderings land on a singular
realization and throw, others succeed. This utility makes the offending degrees of freedom
explicit.

It evaluates the Jacobian of the initialization residual at the initial guess `u0` (via
`ForwardDiff`), computes its singular value decomposition, and reports the unknowns that
participate in the null space — i.e. the directions along which the initialization is free
to move. Pinning those unknowns (or adding equations that constrain them) makes the
initialization well-posed and deterministic.

`prob` may be a problem that carries initialization data (e.g. an `ODEProblem`/`DAEProblem`
built from a `System`), or an initialization `NonlinearProblem` /
`NonlinearLeastSquaresProblem` directly.

# Keyword arguments
  - `rtol`: a singular value `σ` is treated as zero when `σ ≤ max(atol, rtol * σmax)`,
    where `σmax` is the largest singular value. Increase it to also surface
    near-singular directions.
  - `atol`: absolute singular-value threshold (see `rtol`).
  - `var_threshold`: only unknowns whose null-space participation exceeds this value are
    reported.
  - `verbose`: print a human-readable report.

# Returns

A `NamedTuple` `(; jacobian, singular_values, nullity, variables)`. `variables` is a vector
of `unknown => weight` pairs sorted by decreasing `weight`. The weight, in `[0, 1]`, is the
squared norm of the unknown's row in an orthonormal basis of the null space (the diagonal
of the null-space projector): how strongly that unknown participates in the
underdetermined directions, independent of the arbitrary basis chosen within the null
space. A `nullity` of `0` means the initialization Jacobian has full column rank at `u0`
and the initialization is locally well-posed.

!!! note
    The Jacobian is evaluated at a single point (`u0`), so this reports the *local* rank
    structure there. A structurally well-posed initialization can still be numerically
    rank deficient at a particular operating point, and vice versa.
"""
function report_initialization_nullspace(
        prob; rtol = 1e-8, atol = 0.0, var_threshold = 1e-3, verbose = true)
    empty_result = (; jacobian = nothing, singular_values = Float64[],
        nullity = 0, variables = Pair[])
    iprob = _initialization_problem(prob)
    if iprob === nothing
        verbose &&
            @info "No initialization problem to analyze: the system is fully determined by its initial conditions."
        return empty_result
    end
    u0 = state_values(iprob)
    if u0 === nothing || isempty(u0)
        verbose && @info "The initialization problem has no unknowns to solve for."
        return empty_result
    end
    u0 = collect(float.(u0))
    p = parameter_values(iprob)
    f = iprob.f
    nresid = f.resid_prototype !== nothing ? length(f.resid_prototype) : length(u0)
    residual = if SciMLBase.isinplace(f)
        u -> (du = similar(u, nresid); f(du, u, p); du)
    else
        u -> f(u, p)
    end
    J = ForwardDiff.jacobian(residual, u0)
    n = length(u0)
    fact = svd(J; full = true)
    S = fact.S
    σmax = isempty(S) ? zero(eltype(S)) : maximum(S)
    threshold = max(atol, rtol * σmax)
    # A right-singular vector is a null direction when its singular value is below the
    # threshold; columns of `V` beyond `length(S)` (present when n > nresid) have an
    # implied singular value of zero and are always null directions.
    is_null = [(j <= length(S) ? S[j] : zero(σmax)) <= threshold for j in 1:n]
    nullity = count(is_null)
    syms = variable_symbols(iprob)
    # Participation = diagonal of the null-space projector V_null * V_null', i.e. the
    # squared row norm of each unknown over an orthonormal null-space basis. This is
    # invariant to the arbitrary choice of basis within the null space.
    weights = nullity == 0 ? zeros(n) : vec(sum(abs2, @view(fact.V[:, is_null]); dims = 2))
    variables = sort(
        [syms[i] => weights[i] for i in 1:n if weights[i] > var_threshold];
        by = last, rev = true)
    if verbose
        println("Initialization Jacobian null-space report")
        println("  residual Jacobian: $(nresid)×$(n), rank ≈ $(n - nullity), nullity ≈ $nullity")
        if !isempty(S)
            k = min(5, length(S))
            println("  smallest $k singular value(s): ",
                round.(S[(end - k + 1):end], sigdigits = 3))
        end
        if nullity == 0
            println("  Full column rank at u0 — initialization is locally well-posed.")
        else
            println("  Unknowns spanning the null space (participation weight ∈ [0, 1]):")
            for (s, w) in variables
                println("    ", lpad(string(round(w, digits = 4)), 8), "  ", s)
            end
            isempty(variables) &&
                println("  (no individual unknown exceeds var_threshold = $var_threshold)")
        end
    end
    return (; jacobian = J, singular_values = S, nullity, variables)
end

# Return the initialization `NonlinearProblem`/`NonlinearLeastSquaresProblem` carried by
# `prob`, or `prob` itself if it is already a nonlinear problem, or `nothing` if there is
# no initialization problem to analyze.
function _initialization_problem(prob)
    prob isa SciMLBase.AbstractNonlinearProblem && return prob
    f = prob.f
    if hasproperty(f, :initialization_data) && f.initialization_data !== nothing
        return f.initialization_data.initializeprob
    end
    return nothing
end
