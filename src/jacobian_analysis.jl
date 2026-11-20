# Numerical Jacobian analysis utilities shared by `analyze_initialization_jacobian` (rank
# structure of the initialization residual at a point) and `analyze_solution_jacobian`
# (local modes / stiffness of the dynamics along a solution).

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

# Dense Jacobian of `residual` at `u0`: central finite differences when `autodiff ===
# nothing`, otherwise the given ADTypes backend through DifferentiationInterface.
function _dense_jacobian(residual, u0, autodiff)
    return autodiff === nothing ? _central_difference_jacobian(residual, u0) :
        DI.jacobian(residual, autodiff, u0)
end

# `names[i] => w[i]` for every `w[i] > threshold`, sorted by decreasing weight.
function _weighted_pairs(names, w, threshold)
    return sort(
        [names[i] => w[i] for i in eachindex(w) if w[i] > threshold];
        by = last, rev = true
    )
end

"""
    analyze_solution_jacobian(sol; ts = nothing, threshold = 0.05, nmodes = 3, verbose = true, autodiff = nothing, check_state = true)

Analyze the local linearized dynamics along an ODE solution: at each analyzed time point,
evaluate the Jacobian `∂f/∂u` of the right-hand side at `(sol(t), p, t)`, compute its
eigenvalues (the local *modes*), their time constants, and the **participation factors**
that attribute each mode to the unknowns of the system. This answers questions like

  - *which unknowns are responsible for the stiffness of my model?* — the unknowns that
    participate in the fastest modes;
  - *what are the fast and slow time constants of the system, and how do they change
    along the trajectory?*
  - *is the model locally unstable somewhere?* — modes with a positive real part.
  - *which oscillatory modes are lightly damped, and which unknowns do they involve?* —
    the relative damping of the complex eigenvalue pairs.

For a mass-matrix DAE `M u' = f(u, p, t)` with a diagonal mass matrix, the algebraic
unknowns are eliminated by a Schur complement, so the modes are those of the underlying
ODE on the differential unknowns, `J_dd - J_da J_aa⁻¹ J_ad` (scaled by the inverse of the
nonzero diagonal of `M`). A nonsingular non-diagonal mass matrix is handled as `M⁻¹ J`.

# Participation factors

With the eigendecomposition `J = V Λ V⁻¹`, the participation of unknown `i` in mode `k` is
`|V[i, k] (V⁻¹)[k, i]|`, normalized so that the participations of each mode sum to one.
Unlike the eigenvector entries themselves, participation factors are independent of the
scaling of the unknowns and of the eigenvectors. A participation close to one means the
mode is (locally) essentially the dynamics of that single unknown. For a defective or
nearly defective eigenvalue (eigenvalue condition number above `0.01 / √δ`, with `δ` the
relative accuracy of the Jacobian), the eigenvectors are nearly parallel and do not
resolve the mode into contributions of individual unknowns; its participation is `NaN`.

# Arguments

  - `sol`: an `ODESolution` of an `ODEProblem` (including mass-matrix DAEs), typically
    built from a `System`.

The states must be real-valued; solutions with complex states are rejected with an
`ArgumentError`.

# Keywords

  - `ts`: the time points to analyze. Defaults to the saved time points `sol.t` (using
    the saved states `sol.u`); otherwise `sol` is interpolated at `ts`.
  - `threshold`: only unknowns whose participation in a mode exceeds this value are
    reported in `fast_unknowns`/`slow_unknowns` and in the printed report.
  - `nmodes`: the number of fastest and slowest modes listed in the printed report.
  - `verbose`: print a human-readable report. Accepts a `Bool`, a `SciMLLogging` preset,
    or an [`MTKVerbosity`](@ref); the report is controlled by its
    `solution_jacobian_analysis` toggle.
  - `autodiff`: `nothing` (default) uses the analytic Jacobian of the problem if one was
    generated (`ODEProblem(sys, ...; jac = true)`) and `AutoForwardDiff()` otherwise.
    Alternatively an ADTypes backend evaluated through `DifferentiationInterface`. The
    analytic Jacobian and automatic differentiation are taken to be accurate to the
    working precision `ϵ`, a finite-difference backend to `√ϵ`; this accuracy determines
    which eigenvalues are resolved (see below).
  - `check_state`: for a solution with dense output, check that the first saved state is
    consistent with the right-hand side (`M u' ≈ f(u, p, t)`). See the note below on
    `save_idxs`. Set to `false` for solutions whose parameters were modified by callbacks.

# Returns

A `NamedTuple` with fields

  - `t`: the analyzed time points.
  - `unknowns`: the unknowns the modes are expressed in (the differential unknowns for a
    mass-matrix DAE), in the order of the rows of `participation`.
  - `jacobians`: the (reduced) Jacobian at each time point.
  - `eigenvalues`: the eigenvalues at each time point, sorted from fastest to slowest
    (by decreasing `|real(λ)|`).
  - `time_constants`: `1 / |real(λ)|` for each eigenvalue (`Inf` for zero eigenvalues and
    for modes whose real part is zero to the accuracy of the computed spectrum, e.g.
    undamped oscillatory modes).
  - `natural_frequencies`: `|λ|` for each eigenvalue.
  - `damping`: the relative damping `ζ = -real(λ) / |λ|` of each eigenvalue (`1` for a
    decaying real mode, `-1` for a growing real mode, `NaN` for a zero eigenvalue).
  - `zero_eigenvalues`: at each time point, a mask of the eigenvalues that are zero to the
    accuracy of the computed spectrum, e.g. integrators and conservation laws (see below).
  - `eigenvalue_conditions`: the condition number `κ = ‖v‖ ‖w‖ / |w' v|` of each
    eigenvalue, with right and left eigenvectors `v` and `w`.
  - `participation`: at each time point, a matrix whose column `k` holds the participation
    of every unknown in mode `k` (columns ordered as `eigenvalues`). The two eigenvalues of
    a complex-conjugate pair describe one oscillatory mode in a two-dimensional real
    subspace; their participation columns are identical and are typically distributed over
    at least two unknowns (e.g. a position and a velocity). The printed report lists such a
    pair once, as `real(λ) ± imag(λ)im`.
  - `stiffness_ratio`: at each time point, `max |real(λ)| / min |real(λ)|` over the
    decaying modes (`real(λ) < 0`), or `NaN` if there are none.
  - `stiffest_index`: the index of the time point with the largest stiffness ratio.
  - `fast_unknowns` / `slow_unknowns`: `unknown => participation` pairs of the fastest /
    slowest decaying mode at the stiffest time point, sorted by decreasing participation.
  - `max_unstable`: the largest number of modes with a positive real part at any analyzed
    time point.
  - `least_damped`: the oscillatory mode (`imag(λ) ≠ 0`) with the smallest relative damping
    over all analyzed time points, as a `NamedTuple` `(; index, t, eigenvalue, damping,
    natural_frequency, unknowns)`, where `index` is the time-point index and `unknowns`
    holds the `unknown => participation` pairs of that mode, or `nothing` if no analyzed
    time point has an oscillatory mode.

# Examples

```julia
sol = solve(prob, Rodas5P())
res = analyze_solution_jacobian(sol)
res.fast_unknowns          # which unknowns make the model stiff
extrema(res.stiffness_ratio)
res.least_damped           # the least-damped oscillatory mode and its unknowns
```

# Zero eigenvalues and numerical resolution

To first order, a perturbation `E` of the Jacobian moves an eigenvalue `λ` by at most
`κ ‖E‖`, where `κ` is its condition number. With Jacobian entries of relative accuracy `δ`
(see `autodiff`), `‖E‖ ≈ δ S`, where `S` is `‖J‖` for an ODE and, for a mass-matrix DAE,
the scale of the terms entering the Schur complement (so that a reduced Jacobian that is
small because of cancellation is not resolved beyond the accuracy of those terms). The
error bound of `λ` is taken as `100 δ S κ`. This first-order bound is valid only while it
does not reach another eigenvalue; otherwise, the `m` eigenvalues within it form a (nearly)
defective cluster, and the bound is limited to `100 δ^(1/m) S`. For a defective cluster
perturbed into split eigenvalues, `κ` is large and the bound has the order of the
splitting, so that, e.g., the split eigenvalues of a double or quadruple integrator
expressed in a non-modal basis are recognized as zero eigenvalues, while a nearly defective
pair of eigenvalues near `-1` is resolved. An
eigenvalue within its bound of the origin is reported as a zero eigenvalue (integrator,
conservation law), and a real part within its bound of zero is treated as zero; the sign
of such a real part is not resolved, and the mode is neither counted as unstable nor
included in the stiffness ratio. Zero eigenvalues are also excluded from `least_damped`.
Well-conditioned eigenvalues are resolved down to about `100 δ S`.

!!! note
    The analysis is local: eigenvalues of the Jacobian along a nonlinear (or
    time-varying) trajectory describe the linearized dynamics at each point, which is
    what determines stiffness and step-size restrictions of explicit integrators, but
    not necessarily the global behaviour. Parameters are taken from `sol.prob`, so
    parameter changes made by callbacks during the solve are not reflected. The solution
    must contain the complete state vector in problem order, i.e. it must not be solved
    with `save_idxs`. For a problem with a symbolic system (e.g. built from a `System`),
    the solution records `save_idxs`, and such solutions are rejected. Without a symbolic
    system, a `save_idxs` that permutes all unknowns is not recorded; it is detected only
    by the consistency check of `check_state`, which requires dense output and does not
    detect a permutation that leaves the first saved state consistent (e.g. equal initial
    values of the permuted unknowns).
"""
function analyze_solution_jacobian(
        sol::SciMLBase.AbstractODESolution; ts = nothing, threshold = 0.05,
        nmodes = 3, verbose = true, autodiff = nothing, check_state = true
    )
    verbose = MTKVerbosity(verbose)
    prob = sol.prob
    prob isa SciMLBase.AbstractODEProblem || throw(
        ArgumentError(
            "`analyze_solution_jacobian` requires the solution of an `ODEProblem`, got a solution of a `$(nameof(typeof(prob)))`."
        )
    )
    f = prob.f
    p = parameter_values(prob)
    ts_ = ts === nothing ? collect(sol.t) : collect(ts)
    u0 = state_values(prob)
    eltype(u0) <: Real || throw(
        ArgumentError(
            "`analyze_solution_jacobian` supports real-valued states only; the state has element type $(eltype(u0))."
        )
    )
    n = length(u0)
    if hasproperty(sol, :saved_subsystem) && sol.saved_subsystem !== nothing
        throw(
            ArgumentError(
                "`analyze_solution_jacobian` requires the complete state vector; solve without `save_idxs`."
            )
        )
    end
    us = ts === nothing ? sol.u : [sol(t) for t in ts_]
    if !isempty(sol.u) && length(first(sol.u)) != n
        throw(
            ArgumentError(
                "The saved state vector has length $(length(first(sol.u))) but the problem has state dimension $n; solve without `save_idxs`."
            )
        )
    end
    M = f.mass_matrix
    check_state && _check_saved_state(sol, f, u0, p, M)
    syms = _state_names(sol, n)
    diffidxs = _differential_indices(M, n)

    jacobians = Matrix{Float64}[]
    eigenvalues = Vector{ComplexF64}[]
    eigenvalue_conditions = Vector{Float64}[]
    time_constants = Vector{Float64}[]
    natural_frequencies = Vector{Float64}[]
    damping = Vector{Float64}[]
    participation = Matrix{Float64}[]
    stiffness_ratio = Float64[]
    zero_eigenvalues = BitVector[]
    real_parts = Vector{Float64}[]
    nunstable = Int[]
    for (u, t) in zip(us, ts_)
        J, δ = _rhs_jacobian(f, u0, u, p, t, autodiff)
        size(J) == (n, n) || throw(
            DimensionMismatch(
                "The Jacobian at t = $t has size $(size(J)), expected ($n, $n)."
            )
        )
        Jr, S = _reduced_jacobian(J, M, diffidxs, t, δ)
        λ, P, κ = _modes(Jr, δ)
        radius = _eigenvalue_error_bounds(λ, κ, δ, S)
        zm = BitVector(abs.(λ) .<= radius)
        σ = [zm[k] || abs(real(λ[k])) <= radius[k] ? 0.0 : real(λ[k]) for k in eachindex(λ)]
        push!(jacobians, Jr)
        push!(eigenvalues, λ)
        push!(eigenvalue_conditions, κ)
        push!(zero_eigenvalues, zm)
        push!(real_parts, σ)
        push!(time_constants, [iszero(s) ? Inf : 1 / abs(s) for s in σ])
        push!(natural_frequencies, abs.(λ))
        push!(damping, [zm[k] ? NaN : -σ[k] / abs(λ[k]) for k in eachindex(λ)])
        push!(participation, P)
        push!(stiffness_ratio, _stiffness_ratio(σ))
        push!(nunstable, count(>(0), σ))
    end

    names = syms[diffidxs]
    stiffest_index = isempty(stiffness_ratio) ? 0 :
        something(_argmax_skipnan(stiffness_ratio), 1)
    fast_unknowns = Pair[]
    slow_unknowns = Pair[]
    if stiffest_index > 0
        P = participation[stiffest_index]
        decaying = findall(<(0), real_parts[stiffest_index])
        if !isempty(decaying)
            fast_unknowns = _weighted_pairs(names, P[:, first(decaying)], threshold)
            slow_unknowns = _weighted_pairs(names, P[:, last(decaying)], threshold)
        end
    end
    max_unstable = isempty(nunstable) ? 0 : maximum(nunstable)
    least_damped = _least_damped_mode(
        ts_, names, eigenvalues, zero_eigenvalues, damping, participation, threshold
    )

    @SciMLMessage(verbose, :solution_jacobian_analysis) do
        _solution_jacobian_report(
            ts_, names, eigenvalues, time_constants, damping, participation,
            stiffness_ratio, stiffest_index, nunstable, least_damped, threshold, nmodes
        )
    end
    return (;
        t = ts_, unknowns = names, jacobians, eigenvalues, eigenvalue_conditions,
        time_constants, natural_frequencies, damping, zero_eigenvalues, participation,
        stiffness_ratio, stiffest_index, fast_unknowns, slow_unknowns, max_unstable,
        least_damped,
    )
end

# Check that the first saved state of `sol` is consistent with the right-hand side, i.e.
# that `M u'` from the interpolation agrees with `f(u, p, t)`. Neither the solution nor its
# interpolation records a numeric `save_idxs`, so a permutation of all unknowns can only be
# detected this way. The tolerance is loose, since the interpolated derivative has the
# accuracy of the solver. The check requires dense output: without it (e.g. with `saveat`),
# the interpolated derivative is a secant slope and is skipped.
function _check_saved_state(sol, f, u0, p, M)
    isempty(sol.t) && return nothing
    (hasproperty(sol, :dense) && sol.dense) || return nothing
    t = first(sol.t)
    u = first(sol.u)
    du = try
        sol(t, Val{1})
    catch
        return nothing
    end
    rhs = _rhs_function(f, u0, p, t)
    x = u isa Number ? [float(u)] : collect(float.(vec(u)))
    fx = rhs(x)
    dv = du isa Number ? [du] : vec(du)
    Mdu = M isa UniformScaling ? M.λ .* dv : Matrix(M) * dv
    mismatch = LinearAlgebra.norm(Mdu .- fx)
    scale = LinearAlgebra.norm(fx) + LinearAlgebra.norm(Mdu)
    if mismatch > 0.1 * scale + sqrt(eps(Float64)) * (1 + LinearAlgebra.norm(x))
        throw(
            ArgumentError(
                "The saved state at t = $t is inconsistent with the right-hand side of the problem (relative mismatch $(round(mismatch / scale, sigdigits = 3)) between `M u'` and `f(u, p, t)`). This occurs when the solution was computed with `save_idxs` that reorder the unknowns, or when callbacks modified the parameters. Solve without `save_idxs`, or pass `check_state = false` to skip this check."
            )
        )
    end
    return nothing
end

# The oscillatory mode with the smallest relative damping over all time points, or
# `nothing` if there is no oscillatory mode.
function _least_damped_mode(
        ts, names, eigenvalues, zero_eigenvalues, damping, participation, threshold
    )
    best = nothing
    for i in eachindex(eigenvalues)
        λ = eigenvalues[i]
        for k in eachindex(λ)
            (iszero(imag(λ[k])) || zero_eigenvalues[i][k]) && continue
            if best === nothing || damping[i][k] < damping[best[1]][best[2]]
                best = (i, k)
            end
        end
    end
    best === nothing && return nothing
    i, k = best
    return (;
        index = i, t = ts[i], eigenvalue = eigenvalues[i][k], damping = damping[i][k],
        natural_frequency = abs(eigenvalues[i][k]),
        unknowns = _weighted_pairs(names, view(participation[i], :, k), threshold),
    )
end

# Symbolic names of the unknowns of `sol`, or their indices if it was not built from a
# symbolic system.
function _state_names(sol, n)
    syms = try
        collect(variable_symbols(sol))
    catch
        nothing
    end
    return syms === nothing || length(syms) != n ? collect(1:n) : syms
end

# A callable for the right-hand side that accepts dual numbers. With the default
# `AutoSpecialize`, the `f` field of an `ODEFunction` is a function wrapper restricted to the
# argument types of the solver; only then is the unwrapped user function returned. Any other
# SciML function (e.g. an operator-backed `ODEFunction` or a `DynamicalODEFunction`) is
# returned as is, so that its calling convention is preserved.
function _raw_function(g)
    g isa SciMLBase.ODEFunction || return g
    h = SciMLBase.unwrapped_f(g.f)
    return h === g.f ? g : h
end

# The right-hand side as a function of a flat real vector `x`, evaluated with `x` converted
# to the representation of the problem's initial state `u0` (a scalar, or an array of the
# same type and shape with the element type of `x`). For a `SplitFunction` it is `f1 + f2`.
function _rhs_function(f, u0, p, t)
    iip = SciMLBase.isinplace(f)
    eval_rhs = if f isa SciMLBase.SplitFunction
        f1, f2 = _raw_function(f.f1), _raw_function(f.f2)
        iip ? (du, u) -> (f1(du, u, p, t); tmp = similar(du); f2(tmp, u, p, t); du .+= tmp) :
            u -> f1(u, p, t) .+ f2(u, p, t)
    else
        g = _raw_function(f)
        iip ? (du, u) -> g(du, u, p, t) : u -> g(u, p, t)
    end
    return if u0 isa Number
        x -> [eval_rhs(only(x))]
    elseif iip
        x -> (u = _restructure(u0, x); du = similar(u); eval_rhs(du, u); vec(du))
    else
        x -> vec(eval_rhs(_restructure(u0, x)))
    end
end

# Jacobian of the ODE right-hand side `f(u, p, t)` with respect to `u`, and an estimate of
# its relative accuracy `δ`. Without an explicit `autodiff` backend, the analytic Jacobian
# of the problem is used if it has one (except for a `SplitFunction`, whose analytic
# Jacobian covers only `f1`), and ForwardDiff otherwise; both are accurate to the working
# precision. A finite-difference backend is accurate only to about `√ϵ`.
function _rhs_jacobian(f, u0, u, p, t, autodiff)
    T = float(real(eltype(u0)))
    x0 = u isa Number ? [float(u)] : collect(float.(vec(u)))
    precision = max(eps(T), eps(Float64))
    if autodiff === nothing && SciMLBase.has_jac(f) && !(f isa SciMLBase.SplitFunction) &&
            !(u0 isa Number)
        J = if SciMLBase.isinplace(f)
            Jbuf = f.jac_prototype === nothing ? zeros(eltype(u), length(u), length(u)) :
                similar(f.jac_prototype, eltype(u))
            f.jac(Jbuf, u, p, t)
            Jbuf
        else
            f.jac(u, p, t)
        end
        return Matrix{Float64}(J), precision
    end
    backend = autodiff === nothing ? ADTypes.AutoForwardDiff() : autodiff
    J = DI.jacobian(_rhs_function(f, u0, p, t), backend, x0)
    fd = backend isa Union{ADTypes.AutoFiniteDiff, ADTypes.AutoFiniteDifferences}
    δ = fd ? sqrt(precision) : precision
    return Matrix{Float64}(J), δ
end

# `x` in the array type and shape of `template`, with the element type of `x`.
function _restructure(template::StaticArrays.StaticArray, x)
    return StaticArrays.similar_type(template, eltype(x))(x)
end
function _restructure(template, x)
    y = similar(template, eltype(x))
    copyto!(y, x)
    return y
end

# Indices of the differential unknowns: those with a nonzero diagonal mass-matrix entry
# (all unknowns for an identity or nonsingular non-diagonal mass matrix).
function _differential_indices(M, n)
    M isa UniformScaling && return collect(1:n)
    Md = Matrix(M)
    LinearAlgebra.isdiag(Md) || return collect(1:n)
    return findall(!iszero, LinearAlgebra.diag(Md))
end

# The Jacobian `Jr` of the underlying ODE on the differential unknowns, and the scale `S` of
# its absolute error: if the entries of each row of `J` have the relative accuracy `δ`, the
# entries of `Jr` are accurate to about `δ S`. `Jr` is `M⁻¹ J` for an identity/nonsingular
# mass matrix, or the Schur complement `J_dd - J_da J_aa⁻¹ J_ad` (scaled by the inverse of
# the nonzero diagonal of `M`) for a singular diagonal mass matrix. For the Schur
# complement, the error scale is derived from the terms entering it rather than from the
# result, which may be small because of cancellation. The algebraic rows are scaled to
# unit norm, so that the rank test and the error scale do not depend on the scaling of the
# algebraic equations; the algebraic block is rejected when its smallest singular value is
# not resolved by the accuracy of the scaled rows.
function _reduced_jacobian(J, M, diffidxs, t, δ)
    n = size(J, 1)
    norm = LinearAlgebra.norm
    if M isa UniformScaling
        return J ./ M.λ, norm(J) / abs(M.λ)
    end
    Md = Matrix{Float64}(M)
    if !LinearAlgebra.isdiag(Md)
        c = LinearAlgebra.cond(Md)
        c * eps(Float64) < 1.0e-2 || throw(
            ArgumentError(
                "`analyze_solution_jacobian` supports singular mass matrices only when they are diagonal."
            )
        )
        Jr = Md \ J
        # The uncertainty of `J` propagates through `M⁻¹`; the rounding of the solve adds
        # an error of about `cond(M) ϵ ‖Jr‖`, expressed relative to `δ`.
        S = LinearAlgebra.opnorm(inv(Md)) * norm(J) + c * (eps(Float64) / δ) * norm(Jr)
        return Jr, max(S, norm(Jr))
    end
    m = LinearAlgebra.diag(Md)
    d = diffidxs
    length(d) == n && return J ./ m, norm(J ./ m)
    a = setdiff(1:n, d)
    rownorms = [norm(@view J[i, :]) for i in a]
    σmin = any(iszero, rownorms) ? 0.0 :
        minimum(LinearAlgebra.svdvals(J[a, a] ./ rownorms))
    σmin > 100 * δ * sqrt(length(a)) || throw(
        ArgumentError(
            "The Jacobian of the algebraic equations with respect to the algebraic unknowns is singular to the accuracy of the Jacobian at t = $t (smallest singular value $(round(σmin, sigdigits = 3)) after scaling the algebraic equations to unit norm), so the algebraic unknowns cannot be eliminated. The DAE may be of index greater than one there."
        )
    )
    isempty(d) && return zeros(0, 0), 0.0
    Ã = J[a, a] ./ rownorms
    X = Ã \ (J[a, d] ./ rownorms)
    Y = J[d, a] / Ã
    Jr = (J[d, d] .- J[d, a] * X) ./ m[d]
    S = (1 + norm(X)) * (norm(J[d, :]) + norm(Y) * sqrt(length(a))) / minimum(abs, m[d])
    return Jr, max(S, norm(Jr))
end

# Eigenvalues of `J`, sorted from fastest to slowest (decreasing `|real(λ)|`, ties broken
# by decreasing `|λ|`), the matrix of normalized participation factors (column `k` for
# mode `k`) and the eigenvalue condition numbers `κ`. The two eigenvalues of a
# complex-conjugate pair are kept adjacent, the one with positive imaginary part first. The
# participation of a mode whose eigenvalue condition number exceeds `0.01 / √δ` is `NaN`:
# its eigenvector is nearly parallel to others (a defective or nearly defective
# eigenvalue), so the eigenbasis does not resolve that mode into contributions of
# individual unknowns.
function _modes(J, δ)
    n = size(J, 1)
    n == 0 && return ComplexF64[], zeros(0, 0), Float64[]
    # Unsorted, LAPACK returns each conjugate pair of a real matrix consecutively, as exact
    # conjugates with conjugate eigenvectors. Julia's default sorting by (real, imag) would
    # separate repeated pairs, e.g. [a - bi, a - bi, a + bi, a + bi].
    E = LinearAlgebra.eigen(J; sortby = nothing)
    λ = ComplexF64.(E.values)
    V = ComplexF64.(E.vectors)
    W = try
        inv(V)
    catch err
        err isa LinearAlgebra.SingularException || rethrow()
        fill(complex(NaN), n, n)
    end
    κ = Float64[]
    for k in 1:n
        v = @view V[:, k]
        w = @view W[k, :]
        κk = LinearAlgebra.norm(v) * LinearAlgebra.norm(w) / abs(LinearAlgebra.dot(conj(w), v))
        push!(κ, isfinite(κk) ? κk : Inf)
    end
    P = abs.(V .* transpose(W))
    for k in 1:n
        if κ[k] > 1.0e-2 / sqrt(δ)
            P[:, k] .= NaN
        else
            s = sum(@view P[:, k])
            s > 0 && (P[:, k] ./= s)
        end
    end
    groups = _conjugate_groups(λ)
    order = sortperm(
        [(abs(real(λ[first(g)])), abs(λ[first(g)])) for g in groups]; rev = true
    )
    perm = reduce(vcat, groups[order])
    return λ[perm], P[:, perm], κ[perm]
end

# Bounds on the error of the computed eigenvalues `λ` of a matrix whose entries have the
# absolute error `δ S`, with eigenvalue condition numbers `κ`. To first order, a
# perturbation `E` moves the eigenvalue `λ_k` by at most `κ_k ‖E‖`, and a safety factor of
# 100 is applied. The first-order bound is valid only while it does not reach another
# eigenvalue. Otherwise, the eigenvalues within the bound form a (nearly) defective cluster
# of size `m`, whose perturbation is of order `δ^(1/m) S`, and the bound is limited to
# `100 δ^(1/m) S`. For a split defective cluster, the computed condition numbers are of order
# `δ^(-(m - 1) / m)`, so both bounds have the order of the splitting; for an exactly or
# nearly defective cluster computed without splitting, the condition numbers are
# meaningless and the limit applies. An eigenvalue within its bound of the origin is a zero
# eigenvalue; a real part within its bound of zero is not resolved.
function _eigenvalue_error_bounds(λ, κ, δ, S)
    δ = max(δ, eps(Float64))
    bounds = similar(κ)
    for k in eachindex(λ)
        first_order = 100 * δ * S * κ[k]
        m = count(l -> abs(l - λ[k]) <= first_order, λ)
        bounds[k] = m == 1 ? first_order : min(first_order, 100 * δ^(1 / m) * S)
    end
    return bounds
end

# `max |σ| / min |σ|` over the negative (resolved) real parts `σ`.
function _stiffness_ratio(σ)
    rates = [-s for s in σ if s < 0]
    isempty(rates) && return NaN
    return maximum(rates) / minimum(rates)
end

function _argmax_skipnan(x)
    best = nothing
    for i in eachindex(x)
        isnan(x[i]) && continue
        (best === nothing || x[i] > x[best]) && (best = i)
    end
    return best
end

# The human-readable report of `analyze_solution_jacobian`.
function _solution_jacobian_report(
        ts, names, eigenvalues, time_constants, damping, participation,
        stiffness_ratio, stiffest_index, nunstable, least_damped, threshold, nmodes
    )
    io = IOBuffer()
    println(io, "Solution Jacobian mode analysis")
    print(io, "  state dimension $(length(names)), $(length(ts)) time point(s)")
    isempty(ts) || print(io, " in [$(first(ts)), $(last(ts))]")
    println(io)
    finite = filter(!isnan, stiffness_ratio)
    if isempty(finite)
        print(io, "  No decaying modes found.")
    else
        lo, hi = extrema(finite)
        print(
            io, "  stiffness ratio ∈ [", round(lo, sigdigits = 3), ", ",
            round(hi, sigdigits = 3), "], largest at t = ", ts[stiffest_index]
        )
    end
    nu = isempty(nunstable) ? 0 : maximum(nunstable)
    if nu > 0
        i = findfirst(>(0), nunstable)
        print(io, "\n", "  Up to $nu locally unstable mode(s) (real(λ) > 0), first at t = $(ts[i]).")
    end
    if least_damped !== nothing
        ld = least_damped
        print(
            io, "\n  Least-damped oscillatory mode: ζ = ", round(ld.damping, sigdigits = 3),
            ", ωn = ", round(ld.natural_frequency, sigdigits = 4), " at t = ", ld.t
        )
        for (s, w) in ld.unknowns
            print(io, "\n      ", lpad(string(round(w, digits = 3)), 6), "  ", s)
        end
    end
    stiffest_index > 0 || return String(take!(io))
    λ = eigenvalues[stiffest_index]
    τ = time_constants[stiffest_index]
    ζ = damping[stiffest_index]
    P = participation[stiffest_index]
    groups = _conjugate_groups(λ)
    ng = length(groups)
    tstr = ts[stiffest_index]
    if ng <= 2nmodes
        print(io, "\n  Modes at t = $tstr, fastest first:")
        _print_modes(io, groups, λ, τ, ζ, P, names, threshold)
    else
        print(io, "\n  Fastest $nmodes mode(s) at t = $tstr:")
        _print_modes(io, groups[1:nmodes], λ, τ, ζ, P, names, threshold)
        print(io, "\n  Slowest $nmodes mode(s) at t = $tstr:")
        _print_modes(io, groups[(ng - nmodes + 1):ng], λ, τ, ζ, P, names, threshold)
    end
    return String(take!(io))
end

# Indices of the eigenvalues grouped into modes: a complex-conjugate pair forms one
# oscillatory mode, represented by the eigenvalue with positive imaginary part. The
# eigenvalues of a real matrix are either exactly real or come in exactly conjugate pairs
# that are adjacent, both in the unsorted LAPACK output and after the sorting in `_modes`,
# so no tolerance is involved.
function _conjugate_groups(λ)
    groups = Vector{Int}[]
    k = 1
    while k <= length(λ)
        l = λ[k]
        if k < length(λ) && !iszero(imag(l)) && λ[k + 1] == conj(l)
            push!(groups, imag(l) > 0 ? [k, k + 1] : [k + 1, k])
            k += 2
        else
            push!(groups, [k])
            k += 1
        end
    end
    return groups
end

function _print_modes(io, groups, λ, τ, ζ, P, names, threshold)
    for g in groups
        j = first(g)
        l = λ[j]
        lstr = length(g) == 1 ?
            (iszero(imag(l)) ? string(round(real(l), sigdigits = 4)) : string(round(l, sigdigits = 4))) :
            string(round(real(l), sigdigits = 4), " ± ", round(abs(imag(l)), sigdigits = 4), "im")
        print(io, "\n    λ = ", lstr, ", τ = ", round(τ[j], sigdigits = 3))
        iszero(imag(l)) || print(
            io, ", ωn = ", round(abs(l), sigdigits = 4), ", ζ = ", round(ζ[j], sigdigits = 3)
        )
        for (s, w) in _weighted_pairs(names, view(P, :, j), threshold)
            print(io, "\n      ", lpad(string(round(w, digits = 3)), 6), "  ", s)
        end
    end
    return nothing
end
