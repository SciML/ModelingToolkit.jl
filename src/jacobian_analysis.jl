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
    analyze_solution_jacobian(sol; ts = nothing, threshold = 0.05, nmodes = 3, verbose = true, autodiff = nothing)

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
mode is (locally) essentially the dynamics of that single unknown.

# Arguments

  - `sol`: an `ODESolution` of an `ODEProblem` (including mass-matrix DAEs), typically
    built from a `System`.

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
    generated (`ODEProblem(sys, ...; jac = true)`) and central finite differences
    otherwise. Alternatively an ADTypes backend (e.g. `AutoForwardDiff()`) evaluated
    through `DifferentiationInterface`.

# Returns

A `NamedTuple` with fields

  - `t`: the analyzed time points.
  - `unknowns`: the unknowns the modes are expressed in (the differential unknowns for a
    mass-matrix DAE), in the order of the rows of `participation`.
  - `jacobians`: the (reduced) Jacobian at each time point.
  - `eigenvalues`: the eigenvalues at each time point, sorted from fastest to slowest
    (by decreasing `|real(λ)|`).
  - `time_constants`: `1 / |real(λ)|` for each eigenvalue (`Inf` for purely oscillatory
    or constant modes, e.g. those of conservation laws, whose real part is zero to
    rounding).
  - `natural_frequencies`: `|λ|` for each eigenvalue.
  - `damping`: the relative damping `ζ = -real(λ) / |λ|` of each eigenvalue (`1` for a
    decaying real mode, `-1` for a growing real mode, `NaN` for `λ = 0`).
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

!!! note
    The analysis is local: eigenvalues of the Jacobian along a nonlinear (or
    time-varying) trajectory describe the linearized dynamics at each point, which is
    what determines stiffness and step-size restrictions of explicit integrators, but
    not necessarily the global behaviour. Parameters are taken from `sol.prob`, so
    parameter changes made by callbacks during the solve are not reflected. The solution
    must contain the complete state vector in problem order, i.e. it must not be solved
    with `save_idxs`.
"""
function analyze_solution_jacobian(
        sol::SciMLBase.AbstractODESolution; ts = nothing, threshold = 0.05,
        nmodes = 3, verbose = true, autodiff = nothing
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
    n = length(u0)
    if hasproperty(sol, :saved_subsystem) && sol.saved_subsystem !== nothing
        throw(
            ArgumentError(
                "`analyze_solution_jacobian` requires the complete state vector; solve without `save_idxs`."
            )
        )
    end
    us = ts === nothing ? sol.u : [sol(t) for t in ts_]
    if !isempty(us) && length(first(us)) != n
        throw(
            ArgumentError(
                "The saved state vector has length $(length(first(us))) but the problem has state dimension $n; solve without `save_idxs`."
            )
        )
    end
    syms = _state_names(sol, n)
    M = f.mass_matrix
    diffidxs = _differential_indices(M, n)

    jacobians = Matrix{Float64}[]
    eigenvalues = Vector{ComplexF64}[]
    time_constants = Vector{Float64}[]
    natural_frequencies = Vector{Float64}[]
    damping = Vector{Float64}[]
    participation = Matrix{Float64}[]
    stiffness_ratio = Float64[]
    nunstable = Int[]
    for (u, t) in zip(us, ts_)
        J = _rhs_jacobian(f, u0, u, p, t, autodiff)
        Jr = _reduced_jacobian(J, M, diffidxs, t)
        λ, P = _modes(Jr)
        push!(jacobians, Jr)
        push!(eigenvalues, λ)
        tol = _zero_tol(λ)
        push!(time_constants, [abs(real(l)) <= tol ? Inf : 1 / abs(real(l)) for l in λ])
        push!(natural_frequencies, abs.(λ))
        push!(damping, [iszero(l) ? NaN : -real(l) / abs(l) for l in λ])
        push!(participation, P)
        push!(stiffness_ratio, _stiffness_ratio(λ))
        push!(nunstable, count(l -> real(l) > tol, λ))
    end

    names = syms[diffidxs]
    stiffest_index = isempty(stiffness_ratio) ? 0 :
        something(_argmax_skipnan(stiffness_ratio), 1)
    fast_unknowns = Pair[]
    slow_unknowns = Pair[]
    if stiffest_index > 0
        λ = eigenvalues[stiffest_index]
        P = participation[stiffest_index]
        decaying = findall(l -> real(l) < -_zero_tol(λ), λ)
        if !isempty(decaying)
            fast_unknowns = _weighted_pairs(names, P[:, first(decaying)], threshold)
            slow_unknowns = _weighted_pairs(names, P[:, last(decaying)], threshold)
        end
    end
    max_unstable = isempty(nunstable) ? 0 : maximum(nunstable)
    least_damped = _least_damped_mode(
        ts_, names, eigenvalues, damping, participation, threshold
    )

    @SciMLMessage(verbose, :solution_jacobian_analysis) do
        _solution_jacobian_report(
            ts_, names, eigenvalues, time_constants, damping, participation,
            stiffness_ratio, stiffest_index, nunstable, least_damped, threshold, nmodes
        )
    end
    return (;
        t = ts_, unknowns = names, jacobians, eigenvalues, time_constants,
        natural_frequencies, damping, participation, stiffness_ratio, stiffest_index,
        fast_unknowns, slow_unknowns, max_unstable, least_damped,
    )
end

# The oscillatory mode with the smallest relative damping over all time points, or
# `nothing` if there is no oscillatory mode.
function _least_damped_mode(ts, names, eigenvalues, damping, participation, threshold)
    best = nothing
    for i in eachindex(eigenvalues)
        λ = eigenvalues[i]
        for k in eachindex(λ)
            iszero(imag(λ[k])) && continue
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

# Jacobian of the ODE right-hand side `f(u, p, t)` with respect to `u`, differentiated with
# respect to a flat vector `x` and evaluated with `x` converted back to the representation
# of the problem's initial state `u0` (a scalar, or an array of the same type and shape).
# For a `SplitFunction` the analytic Jacobian covers only its first component, so the
# complete right-hand side `f1 + f2` is differentiated instead.
function _rhs_jacobian(f, u0, u, p, t, autodiff)
    x0 = u isa Number ? [float(u)] : collect(float.(vec(u)))
    if autodiff === nothing && SciMLBase.has_jac(f) && !(f isa SciMLBase.SplitFunction) &&
            !(u0 isa Number)
        if SciMLBase.isinplace(f)
            J = f.jac_prototype === nothing ? zeros(eltype(u), length(u), length(u)) :
                similar(f.jac_prototype, eltype(u))
            f.jac(J, u, p, t)
            return Matrix{Float64}(J)
        else
            return Matrix{Float64}(f.jac(u, p, t))
        end
    end
    rhs = if u0 isa Number
        x -> [f(only(x), p, t)]
    elseif SciMLBase.isinplace(f)
        x -> (xs = _restructure(u0, x); du = similar(xs); f(du, xs, p, t); vec(du))
    else
        x -> vec(f(_restructure(u0, x), p, t))
    end
    return Matrix{Float64}(_dense_jacobian(rhs, x0, autodiff))
end

# `x` in the array type and shape of `template`, with the element type of `x`.
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

# The Jacobian of the underlying ODE on the differential unknowns: `M⁻¹ J` for an
# identity/nonsingular mass matrix, or the Schur complement eliminating the algebraic
# unknowns for a singular diagonal mass matrix.
function _reduced_jacobian(J, M, diffidxs, t)
    n = size(J, 1)
    if M isa UniformScaling
        return J ./ M.λ
    end
    Md = Matrix{Float64}(M)
    if !LinearAlgebra.isdiag(Md)
        F = lu(Md; check = false)
        issuccess(F) || throw(
            ArgumentError(
                "`analyze_solution_jacobian` supports singular mass matrices only when they are diagonal."
            )
        )
        return F \ J
    end
    m = LinearAlgebra.diag(Md)
    d = diffidxs
    length(d) == n && return J ./ m
    a = setdiff(1:n, d)
    F = lu(J[a, a]; check = false)
    issuccess(F) || throw(
        ArgumentError(
            "The Jacobian of the algebraic equations with respect to the algebraic unknowns is singular at t = $t, so the algebraic unknowns cannot be eliminated (the DAE is of index > 1 there)."
        )
    )
    return (J[d, d] .- J[d, a] * (F \ J[a, d])) ./ m[d]
end

# Eigenvalues of `J`, sorted from fastest to slowest (decreasing `|real(λ)|`, ties broken
# by decreasing `|λ|`), and the matrix of normalized participation factors (column `k` for
# mode `k`). The two eigenvalues of a complex-conjugate pair are kept adjacent, the one with
# positive imaginary part first.
function _modes(J)
    n = size(J, 1)
    n == 0 && return ComplexF64[], zeros(0, 0)
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
        # Defective Jacobian: the eigenvectors do not form a basis. Fall back to the
        # pseudo-inverse, which still attributes each mode to the unknowns spanning it.
        LinearAlgebra.pinv(V)
    end
    P = abs.(V .* transpose(W))
    for k in 1:n
        s = sum(@view P[:, k])
        s > 0 && (P[:, k] ./= s)
    end
    groups = _conjugate_groups(λ)
    order = sortperm(
        [(abs(real(λ[first(g)])), abs(λ[first(g)])) for g in groups]; rev = true
    )
    perm = reduce(vcat, groups[order])
    return λ[perm], P[:, perm]
end

# Real parts below this are treated as zero (constant/purely oscillatory modes).
_zero_tol(λ) = isempty(λ) ? 0.0 : 1.0e3 * eps(Float64) * max(1.0, maximum(abs, λ))

function _stiffness_ratio(λ)
    tol = _zero_tol(λ)
    rates = [-real(l) for l in λ if real(l) < -tol]
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
