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

Determine which unknowns (state variables) contribute to the fast, slow and poorly damped
modes along an ODE solution. At each analyzed time point, the Jacobian `∂f/∂u` of the
right-hand side is evaluated at `(sol(t), p, t)`; its eigenvalues are the local *modes*,
with time constants `1 / |real(λ)|`, natural frequencies `|λ|` and relative damping
`ζ = -real(λ) / |λ|`, and **participation factors** attribute each mode to the unknowns.

For a mass-matrix DAE `M u' = f(u, p, t)` with a diagonal mass matrix, the algebraic
unknowns are eliminated by a Schur complement, so the modes are those of the underlying
ODE on the differential unknowns, `J_dd - J_da J_aa⁻¹ J_ad` (each row scaled by the inverse
of its mass entry). A nonsingular non-diagonal mass matrix is handled as `M⁻¹ J`. The states
must be real-valued; solutions with complex states are rejected with an `ArgumentError`.

# Participation factors

With right and left eigenvectors `v_k` and `w_k` of the mode `λ_k`, normalized such that
`w_k' v_k = 1`, the participation of unknown `i` is `|v_k[i] w_k[i]|`, normalized so that the
participations of each mode sum to one. Participation factors do not depend on the scaling
of the unknowns (their units) or of the eigenvectors. A participation close to one means
that the mode is (locally) essentially the dynamics of that single unknown. The two
eigenvalues of a complex-conjugate pair describe one oscillatory mode in a two-dimensional
real subspace; their participations are identical and typically distributed over at least
two unknowns (e.g. a position and a velocity).

A defective or nearly defective eigenvalue (e.g. the double zero eigenvalue of a double
integrator, or the double eigenvalue of two identical first-order lags in series) has no
well-defined individual eigenvector. Such eigenvalues are grouped with their nearest
eigenvalues into the smallest cluster whose spectral projector `Π` is well-conditioned, and
the participation of unknown `i` in each mode of the cluster is `|Π[i, i]|` (normalized),
the participation in the invariant subspace of the cluster. The same applies to coincident
eigenvalues (closer than `10 ϵ ‖J‖`, e.g. a repeated eigenvalue of two identical subsystems),
whose individual eigenvectors are not determined; their eigenvalues and damping remain accurate. `cluster_sizes`
reports the size of the cluster of each mode.

# Numerical resolution

The eigenvalues are computed to an accuracy of about `κ ϵ ‖J‖` in suitably scaled
coordinates, where `κ` is the condition number in `eigenvalue_conditions`. Integrators and conservation
laws appear as eigenvalues close to zero, i.e. as the slowest modes; the sign of a real part
below this accuracy is not meaningful. The computed eigenvalues of a cluster of size greater
than one are split by rounding errors, so their individual time constants and damping are
not meaningful either. Modes are ranked by their natural frequency `|λ|`, so that a
lightly damped high-frequency oscillation counts as a fast mode; its time constant may
nevertheless be long.

# Arguments

  - `sol`: an `ODESolution` of an `ODEProblem` (including mass-matrix DAEs, but not a
    `DiscreteProblem` or `ImplicitDiscreteProblem`), typically
    built from a `System`.

# Keywords

  - `ts`: the time points to analyze. Defaults to the saved time points `sol.t` (using
    the saved states `sol.u`); otherwise `sol` is interpolated at `ts`.
  - `threshold`: only unknowns whose participation in a mode exceeds this value are
    reported in `fast_unknowns`, `slow_unknowns`, `least_damped` and the printed report.
  - `nmodes`: the number of fastest and slowest modes listed in the printed report.
  - `verbose`: print a human-readable report. Accepts a `Bool`, a `SciMLLogging` preset,
    or an [`MTKVerbosity`](@ref); the report is controlled by its
    `solution_jacobian_analysis` toggle.
  - `autodiff`: `nothing` (default) uses the analytic Jacobian of the problem if one was
    generated (`ODEProblem(sys, ...; jac = true)`) and `AutoForwardDiff()` otherwise.
    Alternatively an ADTypes backend evaluated through `DifferentiationInterface`. With a
    finite-difference backend, the Jacobian is accurate only to about `√ϵ`; in particular, a
    singular algebraic block of a DAE may then not be detected.

# Returns

A `NamedTuple` with fields

  - `t`: the analyzed time points.
  - `unknowns`: the unknowns the modes are expressed in (the differential unknowns for a
    mass-matrix DAE), in the order of the rows of `participation`.
  - `jacobians`: the (reduced) Jacobian at each time point.
  - `eigenvalues`: the eigenvalues at each time point, sorted from fastest to slowest
    (by decreasing `|λ|`, conjugate pairs adjacent).
  - `eigenvalue_conditions`: the condition number `Σ_i |v[i] w[i]|` of each eigenvalue
    (`Σ_i |Π[i, i]| / m` for a cluster of size `m`), which does not depend on the units of
    the unknowns.
  - `cluster_sizes`: the size of the cluster of each mode (one for a simple, well-conditioned
    eigenvalue).
  - `time_constants`: `1 / |real(λ)|` for each eigenvalue (`Inf` for a zero real part).
  - `natural_frequencies`: `|λ|` for each eigenvalue.
  - `damping`: the relative damping `ζ = -real(λ) / |λ|` of each eigenvalue (`NaN` for
    `λ = 0`).
  - `participation`: at each time point, a matrix whose column `k` holds the participation
    of every unknown in mode `k` (columns ordered as `eigenvalues`).
  - `fast_unknowns` / `slow_unknowns`: at each time point, the `unknown => participation`
    pairs of the fastest / slowest mode, sorted by decreasing participation.
  - `least_damped`: the oscillatory mode (`imag(λ) ≠ 0`, not a member of a cluster of
    split eigenvalues) with the smallest relative damping over all analyzed time points, as a `NamedTuple`
    `(; index, t, eigenvalue, damping, natural_frequency, unknowns)`, where `index` is the
    time-point index and `unknowns` holds the `unknown => participation` pairs, or `nothing`
    if no analyzed time point has an oscillatory mode.

The printed report lists the fastest and slowest modes at the time point with the fastest
mode, and the least-damped oscillatory mode.

# Examples

```julia
sol = solve(prob, Rodas5P())
res = analyze_solution_jacobian(sol)
res.fast_unknowns[end]     # unknowns of the fastest mode at the final time
res.slow_unknowns[end]     # unknowns of the slowest mode at the final time
res.least_damped           # the least-damped oscillatory mode and its unknowns
```

!!! note
    The analysis is local: eigenvalues of the Jacobian along a nonlinear (or
    time-varying) trajectory describe the linearized dynamics at each point. Parameters
    are taken from `sol.prob`, so parameter changes made by callbacks during the solve are
    not reflected. The solution must contain the complete state vector in problem order,
    i.e. it must not be solved with `save_idxs`. For a problem with a symbolic system
    (e.g. built from a `System`), the solution records `save_idxs`, and such solutions are
    rejected, as are solutions that save a subset of the unknowns. A `save_idxs` that
    permutes all unknowns of a problem without a symbolic system is not recorded by the
    solution and cannot be detected.
"""
function analyze_solution_jacobian(
        sol::SciMLBase.AbstractODESolution; ts = nothing, threshold = 0.05,
        nmodes = 3, verbose = true, autodiff = nothing
    )
    verbose = MTKVerbosity(verbose)
    prob = sol.prob
    (prob isa SciMLBase.AbstractODEProblem && !(prob isa SciMLBase.AbstractDiscreteProblem)) || throw(
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
    if !isempty(sol.u) && length(first(sol.u)) != n
        throw(
            ArgumentError(
                "The saved state vector has length $(length(first(sol.u))) but the problem has state dimension $n; solve without `save_idxs`."
            )
        )
    end
    us = ts === nothing ? sol.u : [sol(t) for t in ts_]
    M = f.mass_matrix
    syms = _state_names(sol, n)
    diffidxs = _differential_indices(M, n)
    names = syms[diffidxs]

    jacobians = Matrix{Float64}[]
    eigenvalues = Vector{ComplexF64}[]
    eigenvalue_conditions = Vector{Float64}[]
    cluster_sizes = Vector{Int}[]
    time_constants = Vector{Float64}[]
    natural_frequencies = Vector{Float64}[]
    damping = Vector{Float64}[]
    participation = Matrix{Float64}[]
    resolved_modes = BitVector[]
    fast_unknowns = Vector{Pair}[]
    slow_unknowns = Vector{Pair}[]
    for (u, t) in zip(us, ts_)
        J, δ = _rhs_jacobian(f, u0, u, p, t, autodiff)
        size(J) == (n, n) || throw(
            DimensionMismatch(
                "The Jacobian at t = $t has size $(size(J)), expected ($n, $n)."
            )
        )
        Jr = _reduced_jacobian(J, M, diffidxs, t, δ)
        λ, P, κ, csize, resolved = _modes(Jr)
        push!(jacobians, Jr)
        push!(eigenvalues, λ)
        push!(eigenvalue_conditions, κ)
        push!(cluster_sizes, csize)
        push!(time_constants, [iszero(real(l)) ? Inf : 1 / abs(real(l)) for l in λ])
        push!(natural_frequencies, abs.(λ))
        push!(damping, [iszero(l) ? NaN : -real(l) / abs(l) for l in λ])
        push!(participation, P)
        push!(resolved_modes, resolved)
        push!(fast_unknowns, isempty(λ) ? Pair[] : _weighted_pairs(names, P[:, 1], threshold))
        push!(slow_unknowns, isempty(λ) ? Pair[] : _weighted_pairs(names, P[:, end], threshold))
    end
    least_damped = _least_damped_mode(
        ts_, names, eigenvalues, resolved_modes, damping, participation, threshold
    )

    @SciMLMessage(verbose, :solution_jacobian_analysis) do
        _solution_jacobian_report(
            ts_, names, eigenvalues, cluster_sizes, time_constants, damping, participation,
            least_damped, threshold, nmodes
        )
    end
    return (;
        t = ts_, unknowns = names, jacobians, eigenvalues, eigenvalue_conditions,
        cluster_sizes, time_constants, natural_frequencies, damping, participation,
        fast_unknowns, slow_unknowns, least_damped,
    )
end

# The oscillatory mode with a resolved eigenvalue (not a member of a cluster of split
# eigenvalues) with the smallest relative damping over all time points, or `nothing` if there
# is no such mode.
function _least_damped_mode(
        ts, names, eigenvalues, resolved_modes, damping, participation, threshold
    )
    best = nothing
    for i in eachindex(eigenvalues)
        λ = eigenvalues[i]
        for k in eachindex(λ)
            (iszero(imag(λ[k])) || !resolved_modes[i][k]) && continue
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
# its relative accuracy `δ`. Without an explicit `autodiff` backend, the analytic Jacobian of
# the problem is used if it has one (except for a `SplitFunction`, whose analytic Jacobian
# covers only `f1`), and ForwardDiff otherwise; both are accurate to the working precision.
# A finite-difference backend is accurate only to about `√ϵ`.
function _rhs_jacobian(f, u0, u, p, t, autodiff)
    T = float(real(eltype(u0)))
    x0 = u isa Number ? [float(u)] : collect(float.(vec(u)))
    precision = max(eps(T), eps(Float64))
    if autodiff === nothing && SciMLBase.has_jac(f) && !(f isa SciMLBase.SplitFunction) &&
            !(u0 isa Number)
        J = if SciMLBase.isinplace(f)
            Jbuf = f.jac_prototype === nothing ? zeros(eltype(u), length(u), length(u)) :
                fill!(similar(f.jac_prototype, eltype(u)), 0)
            f.jac(Jbuf, u, p, t)
            Jbuf
        else
            f.jac(u, p, t)
        end
        return Matrix{Float64}(J), precision
    end
    backend = autodiff === nothing ? ADTypes.AutoForwardDiff() : autodiff
    J = DI.jacobian(_rhs_function(f, u0, p, t), backend, x0)
    dense = backend isa ADTypes.AutoSparse ? ADTypes.dense_ad(backend) : backend
    fd = dense isa Union{ADTypes.AutoFiniteDiff, ADTypes.AutoFiniteDifferences}
    return Matrix{Float64}(J), fd ? sqrt(precision) : precision
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

# The Jacobian of the underlying ODE on the differential unknowns: `M⁻¹ J` for an
# identity/nonsingular mass matrix, or the Schur complement `J_dd - J_da J_aa⁻¹ J_ad` (each
# differential row scaled by the inverse of its mass entry) for a singular diagonal mass
# matrix. The algebraic block `J_aa` is rejected when, after scaling its rows and then the
# columns of the row-scaled block to unit norm (so that the test does not depend on the
# scaling of the algebraic equations), its smallest singular value is not resolved by the
# relative accuracy `δ` of the Jacobian.
function _reduced_jacobian(J, M, diffidxs, t, δ)
    n = size(J, 1)
    if M isa UniformScaling
        return J ./ M.λ
    end
    Md = Matrix{Float64}(M)
    if !LinearAlgebra.isdiag(Md)
        LinearAlgebra.cond(Md) * eps(Float64) < 1.0e-2 || throw(
            ArgumentError(
                "`analyze_solution_jacobian` supports singular mass matrices only when they are diagonal."
            )
        )
        return Md \ J
    end
    m = LinearAlgebra.diag(Md)
    d = diffidxs
    length(d) == n && return J ./ m
    a = setdiff(1:n, d)
    Jaa = J[a, a]
    # Equilibrate the rows, then the columns of the row-scaled block.
    rownorms = [LinearAlgebra.norm(@view Jaa[i, :]) for i in axes(Jaa, 1)]
    σmin = 0.0
    if !any(iszero, rownorms)
        Jr = Jaa ./ rownorms
        colnorms = [LinearAlgebra.norm(@view Jr[:, j]) for j in axes(Jr, 2)]
        any(iszero, colnorms) ||
            (σmin = minimum(LinearAlgebra.svdvals(Jr ./ transpose(colnorms))))
    end
    σmin > 100 * δ * sqrt(length(a)) || throw(
        ArgumentError(
            "The Jacobian of the algebraic equations with respect to the algebraic unknowns is singular to the accuracy of the Jacobian at t = $t (smallest singular value $(round(σmin, sigdigits = 3)) after scaling the algebraic equations to unit norm), so the algebraic unknowns cannot be eliminated. The DAE may be of index greater than one there."
        )
    )
    isempty(d) && return zeros(0, 0)
    return (J[d, d] .- J[d, a] * (Jaa \ J[a, d])) ./ m[d]
end

# Modes of `J`: the eigenvalues sorted from fastest to slowest (decreasing `|λ|`, ties broken
# by decreasing `|real(λ)|`, conjugate pairs adjacent with positive imaginary part first),
# the matrix of normalized participation factors (column `k` for mode `k`), the condition
# numbers, the cluster sizes, and whether the eigenvalue of each mode is resolved (not a
# member of a cluster of split eigenvalues). The eigendecomposition is computed for the
# balanced matrix `D⁻¹ J D` (diagonal `D`), which has the same eigenvalues and participation
# factors. The condition numbers `Σ_i |v[i] w[i]|` do not depend on the units of the
# unknowns.
#
# Modes whose condition number is at most `0.01 / √ϵ` and whose eigenvalue is not coincident
# with another one get the participation
# `|v[i] w[i]|` from their right and left eigenvectors. The left eigenvectors are the rows
# of `V⁻¹` when the eigenbasis is well-conditioned, and are otherwise computed from `J'`
# (a defective eigenvalue makes `V` singular, which would corrupt `V⁻¹` for all modes). The
# remaining, ill-conditioned modes are grouped with their nearest eigenvalues (including all
# eigenvalues at the same distance) into the smallest cluster whose spectral projector `Π`
# satisfies `Σ_i |Π[i, i]| / m ≤ 0.01 / √ϵ` (`m` the cluster size); the participation of each
# mode of the cluster is then `|Π[i, i]|` and its condition number `Σ_i |Π[i, i]| / m`, which
# for a single eigenvalue equals its condition number.
function _modes(J)
    n = size(J, 1)
    n == 0 && return ComplexF64[], zeros(0, 0), Float64[], Int[], BitVector()
    B = Matrix{Float64}(J)
    LinearAlgebra.LAPACK.gebal!('S', B)
    # Unsorted, LAPACK returns each conjugate pair of a real matrix consecutively, as exact
    # conjugates with conjugate eigenvectors. Julia's default sorting by (real, imag) would
    # separate repeated pairs, e.g. [a - bi, a - bi, a + bi, a + bi].
    E = LinearAlgebra.eigen(B; sortby = nothing)
    λ = ComplexF64.(E.values)
    V = Matrix{ComplexF64}(E.vectors)
    condV = LinearAlgebra.cond(V)
    wellconditioned = isfinite(condV) && condV * eps(Float64) < 1.0e-2
    W = Matrix{ComplexF64}(wellconditioned ? inv(V) : _left_eigenvectors(B, λ, V))
    κmax = 1.0e-2 / sqrt(eps(Float64))
    P = abs.(V .* transpose(W))
    # `Σ_i |v[i] w[i]|` (with `w' v = 1`) is the minimum of `‖D⁻¹ v‖ ‖D w‖` over diagonal `D`,
    # i.e. the eigenvalue condition number in the most favourable units of the unknowns.
    κ = [(c = sum(@view P[:, k]); isfinite(c) ? c : Inf) for k in 1:n]
    if !wellconditioned
        # A left eigenvector from `J'` matched to the wrong right eigenvector (e.g. of the
        # other one of two identical subsystems with split eigenvalues) is not orthogonal to
        # the right eigenvectors of the other eigenvalues. `W V` does not depend on the
        # scaling of the unknowns. Only row `k` is tested for mode `k`: a large entry in
        # column `k` may stem from the inaccurate left eigenvector of another, defective mode.
        S = W * V
        for k in 1:n
            any(j -> j != k && abs(S[k, j]) > 1.0e-3, 1:n) && (κ[k] = Inf)
        end
    end
    # Eigenvalues closer than the eigenvector perturbation bound `ϵ ‖J‖ / gap` permits do not
    # determine their individual eigenvectors (e.g. a repeated eigenvalue of two identical
    # subsystems, split by rounding errors), so they are always treated as a cluster.
    tolc = 10 * eps(Float64) * LinearAlgebra.norm(B)
    for k in 1:n
        count(l -> abs(l - λ[k]) <= tolc, λ) > 1 && (κ[k] = Inf)
    end
    # The conjugate of each eigenvalue (itself if real); clusters are formed consistently
    # for both eigenvalues of a pair.
    partner = collect(1:n)
    for g in _conjugate_groups(λ)
        length(g) == 2 && (partner[g[1]] = g[2]; partner[g[2]] = g[1])
    end
    for k in 1:n
        isinf(κ[k]) && (κ[partner[k]] = Inf)
    end
    csize = ones(Int, n)
    cluster = zeros(Int, n) # cluster id of each mode, 0 if not clustered
    for k in 1:n
        (κ[k] <= κmax || cluster[k] > 0) && continue
        # Candidate clusters grow by the distance to `λ_k`; all eigenvalues at the same
        # distance are added together, since coincident eigenvalues cannot be assigned to
        # different invariant subspaces by their values.
        distances = abs.(λ .- λ[k])
        for r in vcat(tolc, sort(unique(filter(>(tolc), distances))))
            c = findall(<=(r), distances)
            length(c) == 1 && continue
            # A candidate that touches an existing cluster absorbs that cluster entirely.
            touched = setdiff(unique(cluster[c]), 0)
            c = union(c, findall(in(touched), cluster))
            # A cluster that contains one eigenvalue of a conjugate pair but not the other
            # is merged with its conjugate cluster; otherwise the conjugate cluster, whose
            # projector is `conj(Π)`, is formed together with it.
            cb = partner[c]
            isempty(intersect(c, cb)) || (c = union(c, cb); cb = c)
            m = length(c)
            Π = _spectral_projector(B, λ[c])
            Π === nothing && continue
            # The analogue of `Σ_i |v[i] w[i]|` for the projector of the cluster.
            κΠ = sum(abs, LinearAlgebra.diag(Π)) / m
            (κΠ <= κmax || m == n) || continue
            for (id, cc) in ((k, c), (partner[k], cb)), j in cc
                P[:, j] .= abs.(LinearAlgebra.diag(Π))
                κ[j] = κΠ
                csize[j] = m
                cluster[j] = id
            end
            break
        end
    end
    for k in 1:n
        s = sum(@view P[:, k])
        s > 0 && (P[:, k] ./= s)
    end
    # The eigenvalues of a cluster of coincident eigenvalues (a repeated eigenvalue with
    # independent eigenvectors) are accurate; those of a cluster of split eigenvalues (a
    # defective or nearly defective eigenvalue) are not.
    resolved = BitVector(
        [
            csize[k] == 1 || all(
                j -> abs(λ[j] - λ[k]) <= tolc, findall(==(cluster[k]), cluster)
            ) for k in 1:n
        ]
    )
    groups = _conjugate_groups(λ)
    order = sortperm(
        [(abs(λ[first(g)]), abs(real(λ[first(g)]))) for g in groups]; rev = true
    )
    perm = reduce(vcat, groups[order])
    return λ[perm], P[:, perm], κ[perm], csize[perm], resolved[perm]
end

# Left eigenvectors of `J` for the eigenvalues `λ` with right eigenvectors `V`, as the rows
# of `W` normalized such that `W[k, :]` times `V[:, k]` is one. Each eigenvalue is matched
# with the closest unmatched eigenvalue of `J'`. For a defective or repeated eigenvalue the
# match may be arbitrary and the normalization ill-conditioned; such modes have a large
# condition number and are treated as clusters by `_modes`.
function _left_eigenvectors(J, λ, V)
    n = length(λ)
    E = LinearAlgebra.eigen(Matrix(J'); sortby = nothing)
    μ = ComplexF64.(E.values)
    Y = ComplexF64.(E.vectors)
    W = similar(V)
    used = falses(n)
    for k in 1:n
        j = argmin(i -> used[i] ? Inf : abs(μ[i] - conj(λ[k])), 1:n)
        used[j] = true
        y = @view Y[:, j]
        W[k, :] .= conj.(y) ./ LinearAlgebra.dot(y, @view V[:, k])
    end
    return W
end

# Spectral projector of `J` onto the invariant subspace of the eigenvalues `λc`,
# `Π = X (Y' X)⁻¹ Y'`, where the columns of `X` and `Y` are orthonormal bases of the right
# and left invariant subspaces, obtained from reordered complex Schur decompositions of `J`
# and `J'`. Unlike eigenvectors, invariant subspaces are well defined for defective and
# repeated eigenvalues. If the selected eigenvalues split a defective cluster, the right
# and left invariant subspaces are not complementary, `Y' X` is singular, and `nothing` is
# returned.
function _spectral_projector(J, λc)
    X = _invariant_subspace(complex(J), λc)
    Y = _invariant_subspace(complex(Matrix(J')), conj.(λc))
    G = Y' * X
    LinearAlgebra.cond(G) * eps(Float64) < 1 || return nothing
    return X * (G \ Y')
end

# Orthonormal basis of the invariant subspace of `A` for the eigenvalues closest to `λc`.
function _invariant_subspace(A, λc)
    F = LinearAlgebra.schur(A)
    select = falses(length(F.values))
    for l in λc
        j = argmin(i -> select[i] ? Inf : abs(F.values[i] - l), eachindex(F.values))
        select[j] = true
    end
    F = LinearAlgebra.ordschur(F, select)
    return F.Z[:, 1:length(λc)]
end

# The human-readable report of `analyze_solution_jacobian`: the fastest and slowest modes at
# the time point with the fastest mode, and the least-damped oscillatory mode.
function _solution_jacobian_report(
        ts, names, eigenvalues, cluster_sizes, time_constants, damping, participation,
        least_damped, threshold, nmodes
    )
    io = IOBuffer()
    println(io, "Solution Jacobian mode analysis")
    print(io, "  state dimension $(length(names)), $(length(ts)) time point(s)")
    isempty(ts) || print(io, " in [$(first(ts)), $(last(ts))]")
    nonempty = findall(!isempty, eigenvalues)
    if !isempty(nonempty)
        i = nonempty[argmax([abs(first(eigenvalues[j])) for j in nonempty])]
        λ = eigenvalues[i]
        groups = _conjugate_groups(λ)
        ng = length(groups)
        args = (λ, cluster_sizes[i], time_constants[i], damping[i], participation[i], names, threshold)
        if ng <= 2nmodes
            print(io, "\n  Modes at t = $(ts[i]) (the time point with the fastest mode), fastest first:")
            _print_modes(io, groups, args...)
        else
            print(io, "\n  Fastest $nmodes mode(s) at t = $(ts[i]) (the time point with the fastest mode):")
            _print_modes(io, groups[1:nmodes], args...)
            print(io, "\n  Slowest $nmodes mode(s) at t = $(ts[i]):")
            _print_modes(io, groups[(ng - nmodes + 1):ng], args...)
        end
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

function _print_modes(io, groups, λ, csize, τ, ζ, P, names, threshold)
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
        csize[j] > 1 && print(io, " (cluster of $(csize[j]) eigenvalues)")
        for (s, w) in _weighted_pairs(names, view(P, :, j), threshold)
            print(io, "\n      ", lpad(string(round(w, digits = 3)), 6), "  ", s)
        end
    end
    return nothing
end
