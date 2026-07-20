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
                # Modelica `homotopy(actual, simplified)` nodes are handled per-SCC inside
                # the `SCCNonlinearProblem` constructor: a block whose equations carry them
                # is built as a `HomotopyProblem` and solved by continuation, while the
                # remaining blocks keep their plain Newton solves.
                SCCNonlinearProblem
            else
                @warn """
                `SCCNonlinearProblem` can only be used with `split = true` systems. \
                Simplify your `System` with `split = true` or pass `use_scc = false` to \
                disable this warning
                """
                MTKBase.get_nonlinear_problem_type(isys)
            end
        else
            # No SCC decomposition to hang per-block continuation off, so a
            # `homotopy`-carrying system is solved as a single `HomotopyProblem`.
            MTKBase.get_nonlinear_problem_type(isys)
        end
    else
        NonlinearLeastSquaresProblem
    end
end

"""
    analyze_initialization_jacobian(prob; rtol = 1e-8, atol = 0.0, threshold = 1e-3, verbose = true, autodiff = nothing)

Diagnose rank deficiency of a system's initialization problem by inspecting the singular
value decomposition of its residual Jacobian, and report both the **unknowns** that span
its (right) null space and the **equations** that are locally redundant (its left null
space).

When the initialization system is rank deficient the solve can converge to different
solutions, or fail, depending on the solver and on the (often nondeterministic) ordering of
the assembled unknowns/equations — a common source of intermittent initialization failures:
some orderings land on a singular realization and throw, others succeed. This utility makes
the offending degrees of freedom explicit:

  - **Underdetermined unknowns** (right null space): directions along which the
    initialization is free to move. Pinning these unknowns (or adding equations that
    constrain them) makes the initialization well-posed and deterministic.
  - **Redundant equations** (left null space): combinations of initialization equations
    whose Jacobian rows are linearly dependent at `u0`, i.e. equations that do not locally
    constrain any additional degree of freedom. These explain why an initialization can
    have more equations than unknowns yet still be underdetermined.

It evaluates the Jacobian of the initialization residual at the initial guess `u0` and
computes its SVD. By default the Jacobian is formed by central finite differences, which
reuses the already-compiled residual function; pass an ADTypes backend via `autodiff`
(e.g. `AutoForwardDiff()`) to differentiate through `DifferentiationInterface` instead,
at the cost of compiling the residual for the backend's number types. `prob` may be a problem that carries initialization
data (e.g. an `ODEProblem`/`DAEProblem` built from a `System`), or an initialization
`NonlinearProblem`/`NonlinearLeastSquaresProblem` directly.

# Keyword arguments
  - `rtol`: a singular value `σ` is treated as zero when `σ ≤ max(atol, rtol * σmax)`,
    where `σmax` is the largest singular value. Increase it to also surface near-singular
    directions.
  - `atol`: absolute singular-value threshold (see `rtol`).
  - `threshold`: only unknowns/equations whose participation exceeds this value are
    reported.
  - `verbose`: print a human-readable report.
  - `autodiff`: `nothing` (default) computes the Jacobian by central finite differences
    with step `cbrt(eps)`, whose truncation error is far below the default rank tolerance
    and whose first call avoids recompiling the generated residual (the dominant cost of
    dual-number AD on large systems). Alternatively an ADTypes backend evaluated through
    `DifferentiationInterface`.

# Returns

A `NamedTuple` `(; jacobian, singular_values, rank, nullity, redundancy, underdetermined_unknowns, redundant_equations)`.

  - `nullity = ncols - rank` is the number of underdetermined directions and `redundancy =
    nrows - rank` the number of redundant equations.
  - `underdetermined_unknowns` is a vector of `unknown => weight` pairs and
    `redundant_equations` a vector of `equation => weight` pairs, each sorted by decreasing
    `weight ∈ [0, 1]`. The weight is the diagonal of the corresponding null-space projector
    (the squared row norm over an orthonormal basis of the right/left null space): how
    strongly that unknown/equation participates in the underdetermined/redundant directions,
    independent of the arbitrary basis chosen within the null space.

A `nullity` of `0` means the Jacobian has full column rank at `u0` (no underdetermined
unknowns); a `redundancy` of `0` means full row rank (no redundant equations).

!!! note
    The Jacobian is evaluated at a single point (`u0`), so this reports the *local* rank
    structure there. A structurally well-posed initialization can still be numerically
    rank deficient at a particular operating point, and vice versa.
"""
function analyze_initialization_jacobian(
        prob; rtol = 1.0e-8, atol = 0.0, threshold = 1.0e-3, verbose = true, autodiff = nothing
    )
    empty_result = (;
        jacobian = nothing, singular_values = Float64[], rank = 0,
        nullity = 0, redundancy = 0, underdetermined_unknowns = Pair[],
        redundant_equations = Pair[],
    )
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
    J = autodiff === nothing ? _central_difference_jacobian(residual, u0) :
        DI.jacobian(residual, autodiff, u0)
    nrows, ncols = size(J)
    fact = svd(J; full = true)
    S = fact.S
    σmax = isempty(S) ? zero(eltype(S)) : maximum(S)
    tol = max(atol, rtol * σmax)
    rank = count(>(tol), S)
    # A singular vector is a null direction when its singular value is below the
    # tolerance. Right vectors (columns of `V`) beyond `length(S)` (present when
    # ncols > nrows) and left vectors (columns of `U`) beyond `length(S)` (present when
    # nrows > ncols) have an implied singular value of zero and are always null.
    col_isnull = [(j <= length(S) ? S[j] : zero(σmax)) <= tol for j in 1:ncols]
    row_isnull = [(i <= length(S) ? S[i] : zero(σmax)) <= tol for i in 1:nrows]
    nullity = count(col_isnull)
    redundancy = count(row_isnull)
    # Participation = diagonal of the null-space projector (squared row norm over an
    # orthonormal null-space basis), invariant to the arbitrary basis within the null space.
    col_w = nullity == 0 ? zeros(ncols) :
        vec(sum(abs2, @view(fact.V[:, col_isnull]); dims = 2))
    row_w = redundancy == 0 ? zeros(nrows) :
        vec(sum(abs2, @view(fact.U[:, row_isnull]); dims = 2))
    syms = variable_symbols(iprob)
    eqs = _initialization_equations(iprob)
    underdetermined_unknowns = sort(
        [syms[i] => col_w[i] for i in 1:ncols if col_w[i] > threshold];
        by = last, rev = true
    )
    redundant_equations = sort(
        [(eqs === nothing ? i : eqs[i]) => row_w[i] for i in 1:nrows if row_w[i] > threshold];
        by = last, rev = true
    )
    if verbose
        println("Initialization Jacobian rank analysis")
        println("  residual Jacobian: $(nrows)×$(ncols), rank ≈ $rank, nullity ≈ $nullity, redundancy ≈ $redundancy")
        if !isempty(S)
            k = min(5, length(S))
            println(
                "  smallest $k singular value(s): ",
                round.(S[(end - k + 1):end], sigdigits = 3)
            )
        end
        if nullity == 0
            println("  Full column rank at u0 — no underdetermined unknowns.")
        else
            println("  Underdetermined unknowns (participation ∈ [0, 1]):")
            for (s, w) in underdetermined_unknowns
                println("    ", lpad(string(round(w, digits = 4)), 8), "  ", s)
            end
        end
        if redundancy == 0
            println("  Full row rank at u0 — no redundant equations.")
        else
            println("  Redundant equations (participation ∈ [0, 1]):")
            for (e, w) in redundant_equations
                println("    ", lpad(string(round(w, digits = 4)), 8), "  ", e)
            end
        end
    end
    return (;
        jacobian = J, singular_values = S, rank, nullity, redundancy,
        underdetermined_unknowns, redundant_equations,
    )
end

# Symbolic equations of an initialization problem, in residual order, or `nothing` if not
# available (e.g. a problem not built from a `System`).
function _initialization_equations(iprob)
    f = iprob.f
    hasproperty(f, :sys) || return nothing
    sys = f.sys
    sys === nothing && return nothing
    return equations(sys)
end

# Return the initialization `NonlinearProblem`/`NonlinearLeastSquaresProblem` carried by
# `prob`, or `prob` itself if it is already a nonlinear problem, or `nothing` if there is
# no initialization problem to analyze.
# Dense Jacobian by central finite differences. Reuses the residual function already
# compiled for the element type of `u0`; `ForwardDiff.jacobian` instead triggers a
# dual-number recompilation of the generated residual, which dominates the runtime on
# large systems. The O(h^2) truncation error (h = cbrt(eps)) is far below the default
# singular-value tolerance of the rank analysis.
function _central_difference_jacobian(residual, u0)
    T = eltype(u0)
    up = copy(u0)
    um = copy(u0)
    J = nothing
    for i in eachindex(u0)
        h = cbrt(eps(T)) * max(one(T), abs(u0[i]))
        up[i] = u0[i] + h
        um[i] = u0[i] - h
        col = (residual(up) .- residual(um)) ./ (2h)
        J === nothing && (J = similar(col, length(col), length(u0)))
        J[:, i] = col
        up[i] = u0[i]
        um[i] = u0[i]
    end
    return J === nothing ? zeros(T, 0, length(u0)) : J
end

function _initialization_problem(prob)
    prob isa SciMLBase.AbstractNonlinearProblem && return prob
    f = prob.f
    if hasproperty(f, :initialization_data) && f.initialization_data !== nothing
        return f.initialization_data.initializeprob
    end
    return nothing
end
