using ModelingToolkitBase, Test
using ModelingToolkitBase: unwrap, complete, unknowns, default_toterm
using ModelingToolkitBase: has_array_equations, accepts_array_equations
using Symbolics
using SciMLBase
using OrdinaryDiffEqBDF: DFBDF
using Sundials: IDA
using DiffEqBase: BrownFullBasicInit

# A system whose interior is written as one array equation over slices, as produced by a
# finite-difference PDE discretization that does not scalarize.
function heat_array_system(n)
    @independent_variables t
    @variables u(t)[1:n]
    D = Differential(t)
    dx = 1 / (n - 1)
    # Residual (cardinalized) form, as a finite-difference discretization emits it: the
    # derivative sits inside the expression rather than being the equation's whole LHS.
    lap = (u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]) ./ dx^2
    interior = broadcast(-, D(u[2:(n - 1)]), lap) ~ zeros(n - 2)
    eqs = [interior, u[1] ~ 0.0, u[n] ~ 0.0]
    @named sys = System(eqs, t, collect(u), [])
    return complete(sys), u, t, D
end

@testset "array equations reach DAEProblem" begin
    n = 11
    sys, u, t, D = heat_array_system(n)
    xs = range(0.0, 1.0, length = n)
    op = vcat(
        [u[i] => sinpi(xs[i]) for i in 1:n],
        [D(u[i]) => 0.0 for i in 1:n]
    )

    prob = DAEProblem(sys, op, (0.0, 0.1); build_initializeprob = false)

    # one output row per element of the array equation, not one per equation
    @test length(prob.u0) == n
    @test prob.u0 isa Vector{Float64}

    # the interior points are differential, the two boundary points algebraic
    @test prob.differential_vars !== nothing
    @test count(prob.differential_vars) == n - 2

    # the residual evaluates: no `Differential` survives into the generated code
    out = zeros(n)
    du = zeros(n)
    prob.f(out, du, prob.u0, prob.p, 0.0)
    @test all(isfinite, out)
    # with du = 0 the interior residual is minus the Laplacian, which is nonzero here
    @test any(!iszero, out)
end

@testset "DAE initialization accepts derivative guesses" begin
    @independent_variables t
    @variables y(t)[1:3] z(t)
    D = Differential(t)
    @named sys = System(
        [zeros(3) ~ D(y[1:3]) .+ z .* y[1:3], 0 ~ z - sum(y)],
        t, [collect(y); z], []
    )
    sys = complete(sys)
    y0 = [1.0, 2.0, 3.0]
    op = [y[i] => y0[i] for i in 1:3]

    expected_du = -sum(y0) .* y0
    for derivative_guesses in (
            [D(y[i]) => 0.5 for i in 1:3],
            [D(y) => fill(0.5, 3)],
            [D(y) => 0.5],
            [],
        )
        guesses = [z => 0.0; derivative_guesses]
        prob = DAEProblem(sys, op, (0.0, 0.1); guesses)
        @test prob.u0[4] == 0.0
        @test prob.du0[1:3] == fill(isempty(derivative_guesses) ? 0.0 : 0.5, 3)
        initprob = prob.f.initialization_data.initializeprob
        init_sol = solve(initprob)

        @test SciMLBase.successful_retcode(init_sol)
        init_unknowns = unknowns(initprob.f.sys)
        z_index = findfirst(isequal(z), init_unknowns)
        for i in 1:3
            idx = findfirst(isequal(default_toterm(unwrap(D(y[i])))), init_unknowns)
            @test idx !== nothing
            idx === nothing || @test init_sol.u[idx] ≈ expected_du[i]
        end
        @test init_sol.u[z_index] ≈ sum(y0)

        residual = zeros(4)
        for initializealg in (nothing, SciMLBase.OverrideInit(), BrownFullBasicInit())
            integ = initializealg === nothing ?
                init(prob, DFBDF()) :
                init(prob, DFBDF(); initializealg)
            @test integ.du[1:3] ≈ expected_du
            prob.f(residual, integ.du, integ.u, integ.p, 0.0)
            @test maximum(abs, residual) < 1.0e-8
        end
    end

    # a derivative fixed in `op` for only part of an array keeps its value in
    # `du0` and in the integrator's `du`; the remaining elements are solved for
    op_partial = [op; D(y[1]) => expected_du[1]]
    for derivative_guesses in (
            [D(y) => 0.5],
            [D(y) => fill(0.5, 3)],
            [D(y[2]) => 0.5, D(y[3]) => 0.5],
            [],
        )
        guesses = [z => 0.0; derivative_guesses]
        prob = DAEProblem(sys, op_partial, (0.0, 0.1); guesses)
        @test prob.du0[1] == expected_du[1]
        residual = zeros(4)
        for initializealg in (nothing, SciMLBase.OverrideInit(), BrownFullBasicInit())
            integ = initializealg === nothing ?
                init(prob, DFBDF()) :
                init(prob, DFBDF(); initializealg)
            @test integ.du ≈ [expected_du; 0.0]
            prob.f(residual, integ.du, integ.u, integ.p, 0.0)
            @test maximum(abs, residual) < 1.0e-8
        end
    end

    @variables yk(t)[1:3] zk(t)
    @parameters k
    sysk = complete(
        System(
            [zeros(3) ~ D(yk[1:3]) .+ k * zk .* yk[1:3], 0 ~ zk - sum(yk)],
            t, [collect(yk); zk], [k]; name = :partial_fixed_derivative
        )
    )
    opk = [[yk[i] => Float64(i) for i in 1:3]; k => 1.0; D(yk[1]) => -6.0]
    probk = DAEProblem(
        sysk, opk, (0.0, 0.1); guesses = [zk => 0.0, D(yk[2]) => 0.5, D(yk[3]) => 0.5]
    )
    @test probk.du0[1] == -6.0
    integk = init(probk, DFBDF())
    @test integk.du ≈ [-6.0, -12.0, -18.0, 0.0]
    residualk = zeros(4)
    probk.f(residualk, integk.du, integk.u, integk.p, 0.0)
    @test maximum(abs, residualk) < 1.0e-8

    # a scalar guess for an array derivative broadcasts to the array shape; storing a
    # scalar under the array-shaped key breaks `get_possibly_indexed` readback
    prob = DAEProblem(
        sys, op, (0.0, 0.1); guesses = [z => 0.0, D(y) => 0.5, D(z) => 0.0],
        build_initializeprob = false
    )
    @test prob.du0[1:3] == fill(0.5, 3)

    # the same broadcast fills only the elements an operating-point entry did not
    prob = DAEProblem(
        sys, op_partial, (0.0, 0.1); guesses = [z => 0.0, D(y) => 0.5, D(z) => 0.0],
        build_initializeprob = false
    )
    @test prob.du0[1:3] == [expected_du[1], 0.5, 0.5]

    # under `FullSpecialize` the solved derivatives are written back through the
    # generated `initializeprobpmap` rather than a wrapper closure
    prob = DAEProblem{true, SciMLBase.FullSpecialize}(
        sys, op, (0.0, 0.1); guesses = [z => 0.0]
    )
    integ = init(prob, DFBDF())
    @test integ.du ≈ [expected_du; 0.0]

    # omitted derivative values with no initialization problem still error, matching
    # the pre-change contract
    @test_throws ModelingToolkitBase.MissingVariablesError DAEProblem(
        sys, op, (0.0, 0.1); guesses = [z => 0.0], build_initializeprob = false
    )
    # an explicit `missing_guess_value = Error()` opts out of the zero starting
    # guess even when an initialization problem would solve for them
    @test_throws ModelingToolkitBase.MissingVariablesError DAEProblem(
        sys, op, (0.0, 0.1); guesses = [z => 0.0],
        missing_guess_value = MissingGuessValue.Error()
    )

    fixed_du0 = [D(y[i]) => 0.25 for i in 1:3]
    prob = DAEProblem(
        sys, [op; [z => 0.0, D(z) => 0.0]; fixed_du0], (0.0, 0.1);
        build_initializeprob = false
    )
    @test prob.du0[1:3] == fill(0.25, 3)
end

@testset "operating-point derivatives win over derivative guesses" begin
    @independent_variables t
    D = Differential(t)
    @variables x1(t) x2(t) w(t) y(t)[1:3] z(t) v(t) a(t) b(t) xc(t)
    @parameters k
    ss = complete(
        System(
            [0 ~ D(x1) + k * w * x1, 0 ~ D(x2) + k * w * x2, 0 ~ w - x1 - x2],
            t, [x1, x2, w], [k]; name = :scalar_dae
        )
    )
    sa = complete(
        System(
            [zeros(3) ~ D(y[1:3]) .+ k * z .* y[1:3], 0 ~ z - sum(y)],
            t, [collect(y); z], [k]; name = :array_dae
        )
    )
    s1 = complete(System([0 ~ D(v) + v], t, [v], []; name = :linear_dae))
    s2 = complete(System([0 ~ D(a) - b, 0 ~ D(b) + a], t, [a, b], []; name = :oscillator))
    sc = complete(System([0 ~ D(xc)^3 + D(xc) - xc], t, [xc], []; name = :implicit_cubic))
    heat, u, _, Dh = heat_array_system(11)

    ops = [x1 => 1.0, x2 => 2.0, k => 1.0]
    opa = [[y[i] => Float64(i) for i in 1:3]; k => 1.0]
    # consistent `integ.du`: w = x1 + x2 = 3 and z = sum(y) = 6; the algebraic slot keeps
    # its `du0` value, which is 0 in every case below
    dus = [-3.0, -6.0, 0.0]
    dua = [-6.0, -12.0, -18.0, 0.0]
    fixa = [D(y[i]) => dua[i] for i in 1:3]
    xs = range(0.0, 1.0, length = 11)
    u0h = sinpi.(xs)
    duh = [0.0; (u0h[1:9] .- 2 .* u0h[2:10] .+ u0h[3:11]) ./ step(xs)^2; 0.0]
    x1ˍt = default_toterm(unwrap(D(x1)))
    yˍt = default_toterm(unwrap(D(y)))

    # (label, system, op, guesses, expected `prob.du0`, `du` slots fixed by `op`, expected
    # `integ.du` or `nothing` when the fixed values are inconsistent)
    cases = [
        ("scalar both fixed, both guessed", ss, [ops; D(x1) => -3.0; D(x2) => -6.0], [w => 0.0, D(x1) => 9.0, D(x2) => 9.0], [-3.0, -6.0, 0.0], [1, 2], dus),
        ("scalar x1 fixed and guessed", ss, [ops; D(x1) => -3.0], [w => 0.0, D(x1) => 9.0], [-3.0, 0.0, 0.0], [1], dus),
        ("scalar all fixed incl D(w), guessed", ss, [ops; D(x1) => -3.0; D(x2) => -6.0; D(w) => 0.0], [w => 0.0, D(x1) => 9.0, D(x2) => 9.0], [-3.0, -6.0, 0.0], [1, 2, 3], dus),
        ("scalar x1 fixed, x2 guessed", ss, [ops; D(x1) => -3.0], [w => 0.0, D(x2) => 0.5], [-3.0, 0.5, 0.0], [1], dus),
        ("scalar x1 fixed, x2 omitted", ss, [ops; D(x1) => -3.0], [w => 0.0], [-3.0, 0.0, 0.0], [1], dus),
        ("scalar both guessed", ss, ops, [w => 0.0, D(x1) => 0.5, D(x2) => 0.5], [0.5, 0.5, 0.0], Int[], dus),
        ("scalar both omitted", ss, ops, [w => 0.0], [0.0, 0.0, 0.0], Int[], dus),
        ("scalar x1 guessed, x2 omitted", ss, ops, [w => 0.0, D(x1) => 0.5], [0.5, 0.0, 0.0], Int[], dus),
        ("scalar both fixed", ss, [ops; D(x1) => -3.0; D(x2) => -6.0], [w => 0.0], [-3.0, -6.0, 0.0], [1, 2], dus),
        ("scalar x2 fixed, x1 and w guessed", ss, [ops; D(x2) => -6.0], [w => 0.0, D(x1) => 0.5, D(w) => 0.0], [0.5, -6.0, 0.0], [2], dus),
        ("scalar toterm key fixed and guessed", ss, [ops; x1ˍt => -3.0], [w => 0.0, D(x1) => 9.0], [-3.0, 0.0, 0.0], [1], dus),
        ("scalar x1 inconsistent", ss, [ops; D(x1) => 5.0], [w => 0.0], [5.0, 0.0, 0.0], [1], nothing),
        ("scalar x1 inconsistent and guessed", ss, [ops; D(x1) => 5.0], [w => 0.0, D(x1) => 9.0], [5.0, 0.0, 0.0], [1], nothing),
        ("array per-element fixed, scalar guess", sa, [opa; fixa], [z => 0.0, D(y) => 0.5], dua, [1, 2, 3], dua),
        ("array per-element fixed, per-element guesses", sa, [opa; fixa], [z => 0.0, [D(y[i]) => 9.0 for i in 1:3]...], dua, [1, 2, 3], dua),
        ("array whole fixed, vector guess", sa, [opa; D(y) => dua[1:3]], [z => 0.0, D(y) => fill(0.5, 3)], dua, [1, 2, 3], dua),
        ("array whole fixed incl D(z), scalar guess", sa, [opa; D(y) => dua[1:3]; D(z) => 0.0], [z => 0.0, D(y) => 0.5], dua, [1, 2, 3, 4], dua),
        ("array whole fixed, no derivative guess", sa, [opa; D(y) => dua[1:3]], [z => 0.0], dua, [1, 2, 3], dua),
        ("array whole fixed, per-element guess", sa, [opa; D(y) => dua[1:3]], [z => 0.0, D(y[2]) => 9.0], dua, [1, 2, 3], dua),
        ("array toterm key fixed, scalar guess", sa, [opa; yˍt => dua[1:3]], [z => 0.0, D(y) => 0.5], dua, [1, 2, 3], dua),
        ("array y[1] fixed and guessed", sa, [opa; D(y[1]) => -6.0], [z => 0.0, D(y[1]) => 9.0], [-6.0, 0.0, 0.0, 0.0], [1], dua),
        ("array per-element fixed incl D(z), per-element guesses", sa, [opa; fixa; D(z) => 0.0], [z => 0.0, [D(y[i]) => 9.0 for i in 1:3]...], dua, [1, 2, 3, 4], dua),
        ("array per-element fixed incl D(z)", sa, [opa; fixa; D(z) => 0.0], [z => 0.0], dua, [1, 2, 3, 4], dua),
        ("array per-element fixed, vector guess", sa, [opa; fixa], [z => 0.0, D(y) => fill(0.5, 3)], dua, [1, 2, 3], dua),
        ("array y[2] fixed, rest omitted", sa, [opa; D(y[2]) => -12.0], [z => 0.0], [0.0, -12.0, 0.0, 0.0], [2], dua),
        ("array y[2] fixed, scalar guess", sa, [opa; D(y[2]) => -12.0], [z => 0.0, D(y) => 0.5], [0.5, -12.0, 0.5, 0.0], [2], dua),
        ("array y[1], y[3] fixed, y[2] guessed", sa, [opa; D(y[1]) => -6.0; D(y[3]) => -18.0], [z => 0.0, D(y[2]) => 0.5], [-6.0, 0.5, -18.0, 0.0], [1, 3], dua),
        ("array y[1], y[3] fixed, y[2] omitted", sa, [opa; D(y[1]) => -6.0; D(y[3]) => -18.0], [z => 0.0], [-6.0, 0.0, -18.0, 0.0], [1, 3], dua),
        ("array y[3] fixed, y[1] guessed", sa, [opa; D(y[3]) => -18.0], [z => 0.0, D(y[1]) => 0.5], [0.5, 0.0, -18.0, 0.0], [3], dua),
        ("array y[1] fixed incl D(z), scalar guess", sa, [opa; D(y[1]) => -6.0; D(z) => 0.0], [z => 0.0, D(y) => 0.5], [-6.0, 0.5, 0.5, 0.0], [1, 4], dua),
        ("array y[1] fixed, conflicting vector guess", sa, [opa; D(y[1]) => -6.0], [z => 0.0, D(y) => [100.0, 0.5, 0.5]], [-6.0, 0.5, 0.5, 0.0], [1], dua),
        ("array y[1] fixed, scalar guess", sa, [opa; D(y[1]) => -6.0], [z => 0.0, D(y) => 0.5], [-6.0, 0.5, 0.5, 0.0], [1], dua),
        ("array y[1] fixed, rest omitted", sa, [opa; D(y[1]) => -6.0], [z => 0.0], [-6.0, 0.0, 0.0, 0.0], [1], dua),
        ("array y[1] fixed, per-element guesses", sa, [opa; D(y[1]) => -6.0], [z => 0.0, D(y[2]) => 0.5, D(y[3]) => 0.5], [-6.0, 0.5, 0.5, 0.0], [1], dua),
        ("array nothing fixed, y[2] guessed", sa, opa, [z => 0.0, D(y[2]) => 0.5], [0.0, 0.5, 0.0, 0.0], Int[], dua),
        ("array nothing fixed, all omitted", sa, opa, [z => 0.0], [0.0, 0.0, 0.0, 0.0], Int[], dua),
        ("array per-element guesses", sa, opa, [z => 0.0, [D(y[i]) => 0.5 for i in 1:3]...], [0.5, 0.5, 0.5, 0.0], Int[], dua),
        ("array vector guess", sa, opa, [z => 0.0, D(y) => fill(0.5, 3)], [0.5, 0.5, 0.5, 0.0], Int[], dua),
        ("array scalar guess", sa, opa, [z => 0.0, D(y) => 0.5], [0.5, 0.5, 0.5, 0.0], Int[], dua),
        ("array zero vector guess", sa, opa, [z => 0.0, D(y) => zeros(3)], zeros(4), Int[], dua),
        ("array zero scalar guess", sa, opa, [z => 0.0, D(y) => 0.0], zeros(4), Int[], dua),
        ("array nonuniform vector guess", sa, opa, [z => 0.0, D(y) => [0.2, 0.4, 0.6]], [0.2, 0.4, 0.6, 0.0], Int[], dua),
        ("array per-element inconsistent incl D(z)", sa, [opa; D(y[1]) => 5.0; D(y[2]) => -12.0; D(y[3]) => -18.0; D(z) => 0.0], [z => 0.0], [5.0, -12.0, -18.0, 0.0], [1, 2, 3, 4], nothing),
        ("array y[1] inconsistent, per-element guesses", sa, [opa; D(y[1]) => 5.0], [z => 0.0, D(y[2]) => 0.5, D(y[3]) => 0.5], [5.0, 0.5, 0.5, 0.0], [1], nothing),
        ("array y[1] inconsistent, rest omitted", sa, [opa; D(y[1]) => 5.0], [z => 0.0], [5.0, 0.0, 0.0, 0.0], [1], nothing),
        ("array whole inconsistent", sa, [opa; D(y) => [5.0, -12.0, -18.0]], [z => 0.0], [5.0, -12.0, -18.0, 0.0], [1, 2, 3], nothing),
        ("array fixed with fixed z inconsistent, guessed", sa, [opa; z => 0.0; [D(y[i]) => 0.25 for i in 1:3]], [D(y[i]) => 9.0 for i in 1:3], [0.25, 0.25, 0.25, 0.0], [1, 2, 3], nothing),
        ("linear scalar guessed", s1, [v => 1.0], [D(v) => 0.5], [0.5], Int[], [-1.0]),
        ("linear scalar omitted", s1, [v => 1.0], Pair[], [0.0], Int[], [-1.0]),
        ("linear scalar fixed and guessed", s1, [v => 1.0, D(v) => -1.0], [D(v) => 0.5], [-1.0], [1], [-1.0]),
        ("oscillator guessed", s2, [a => 1.0, b => 2.0], [D(a) => 0.0, D(b) => 0.0], [0.0, 0.0], Int[], [2.0, -1.0]),
        ("implicit cubic guessed", sc, [xc => 2.0], [D(xc) => 0.5], [0.5], Int[], [1.0]),
        ("implicit cubic fixed and guessed", sc, [xc => 2.0, D(xc) => 1.0], [D(xc) => 0.5], [1.0], [1], [1.0]),
        ("heat scalar zero guess", heat, [u[i] => u0h[i] for i in 1:11], [Dh(u) => 0.0], zeros(11), Int[], duh),
        ("heat vector zero guess", heat, [u[i] => u0h[i] for i in 1:11], [Dh(u) => zeros(11)], zeros(11), Int[], duh),
        ("heat per-element zero guesses", heat, [u[i] => u0h[i] for i in 1:11], [Dh(u[i]) => 0.0 for i in 1:11], zeros(11), Int[], duh),
        ("heat omitted", heat, [u[i] => u0h[i] for i in 1:11], Pair[], zeros(11), Int[], duh),
    ]
    initializealgs = (
        ("DFBDF", DFBDF, nothing),
        ("DFBDF + OverrideInit", DFBDF, SciMLBase.OverrideInit()),
        ("IDA", IDA, nothing),
        ("IDA + OverrideInit", IDA, SciMLBase.OverrideInit()),
    )
    # tight integrator tolerances so the initialization solve converges well below the
    # residual threshold
    function init_du(prob, solver, initializealg)
        tols = (; abstol = 1.0e-10, reltol = 1.0e-10)
        integ = initializealg === nothing ? init(prob, solver(); tols...) :
            init(prob, solver(); initializealg, tols...)
        return collect(integ.du), collect(integ.u), integ.p
    end
    for (label, sys, op, guesses, du0, fixed, du) in cases
        @testset "$label" begin
            prob = DAEProblem(sys, op, (0.0, 0.2); guesses, warn_initialize_determined = false)
            @test prob.du0 == du0
            for (algname, solver, initializealg) in initializealgs
                @testset "$algname" begin
                    integ_du, integ_u, integ_p = init_du(prob, solver, initializealg)
                    @test integ_du[fixed] ≈ du0[fixed]
                    if du !== nothing
                        @test integ_du ≈ du
                        residual = zeros(length(du))
                        prob.f(residual, integ_du, integ_u, integ_p, 0.0)
                        @test maximum(abs, residual) < 1.0e-8
                    end
                end
            end
        end
    end

    # a changed state or parameter re-solves the derivatives the operating point left free
    prob = DAEProblem(
        sa, opa, (0.0, 0.2); guesses = [z => 0.0, D(y) => 0.5], warn_initialize_determined = false
    )
    for (label, newprob, expected) in (
            ("numeric u0", remake(prob; u0 = [4.0, 5.0, 6.0, 0.0]), -15.0 .* [4.0, 5.0, 6.0]),
            ("symbolic u0", remake(prob; u0 = [y[i] => 3.0 + i for i in 1:3]), -15.0 .* [4.0, 5.0, 6.0]),
            ("parameter", remake(prob; p = [k => 2.0]), -12.0 .* [1.0, 2.0, 3.0]),
        )
        for (algname, solver, initializealg) in initializealgs
            @testset "remake $label: $algname" begin
                integ_du, integ_u, integ_p = init_du(newprob, solver, initializealg)
                @test integ_du[1:3] ≈ expected
                residual = zeros(4)
                newprob.f(residual, integ_du, integ_u, integ_p, 0.0)
                @test maximum(abs, residual) < 1.0e-8
            end
        end
    end

    # `CheckInit` consumes `du0` directly, so it accepts consistent fixed derivatives even
    # when conflicting guesses are supplied, and rejects inconsistent ones
    opz = [opa; z => 6.0]
    for (label, sys, op, guesses) in (
            ("array per-element fixed", sa, [opz; fixa; D(z) => 0.0], Pair[]),
            ("array whole fixed", sa, [opz; D(y) => dua[1:3]; D(z) => 0.0], Pair[]),
            ("array per-element fixed, D(z) omitted", sa, [opz; fixa], Pair[]),
            ("array per-element fixed, conflicting guess", sa, [opz; fixa; D(z) => 0.0], [D(y) => 9.0]),
            ("scalar fixed", ss, [ops; w => 3.0; D(x1) => -3.0; D(x2) => -6.0; D(w) => 0.0], Pair[]),
            ("scalar fixed, conflicting guesses", ss, [ops; w => 3.0; D(x1) => -3.0; D(x2) => -6.0; D(w) => 0.0], [D(x1) => 9.0, D(x2) => 9.0]),
        )
        @testset "CheckInit $label" begin
            prob = DAEProblem(sys, op, (0.0, 0.2); guesses, warn_initialize_determined = false)
            integ = init(prob, DFBDF(); initializealg = SciMLBase.CheckInit())
            @test collect(integ.du) ≈ (sys === sa ? dua : dus)
        end
    end
    prob = DAEProblem(
        sa, [opz; D(y[1]) => 5.0; D(y[2]) => -12.0; D(y[3]) => -18.0; D(z) => 0.0],
        (0.0, 0.2); warn_initialize_determined = false
    )
    @test_throws SciMLBase.CheckInitFailureError init(
        prob, DFBDF(); initializealg = SciMLBase.CheckInit()
    )
end

@testset "array unknowns flatten under an array equation" begin
    n = 11
    @independent_variables t
    @variables u(t)[1:n]
    D = Differential(t)
    @named sys = System([D(u) ~ -u], t, [u], [])
    sys = complete(sys)
    op = [u => zeros(n), D(u) => zeros(n)]
    # `flat_unknowns` flattens `u`, so the array equation over the array unknown
    # constructs directly
    prob = ODEProblem(sys, op, (0.0, 0.1); build_initializeprob = false)
    @test length(prob.u0) == n
end

@testset "array-equation DAE solves to the analytic solution" begin
    n = 21
    sys, u, t, D = heat_array_system(n)
    xs = range(0.0, 1.0, length = n)
    op = vcat(
        [u[i] => sinpi(xs[i]) for i in 1:n],
        [D(u[i]) => 0.0 for i in 1:n]
    )
    tend = 0.1
    prob = DAEProblem(sys, op, (0.0, tend); build_initializeprob = false)
    # `du0` above is not consistent; the solver's own DAE initialization supplies it.
    sol = solve(
        prob, DFBDF(); initializealg = BrownFullBasicInit(),
        reltol = 1.0e-8, abstol = 1.0e-8, saveat = [tend]
    )
    @test SciMLBase.successful_retcode(sol)
    exact = [exp(-pi^2 * tend) * sinpi(xi) for xi in xs]
    # second-order spatial discretization on 21 points
    @test maximum(abs.(sol.u[end] .- exact)) < 5.0e-3
end

@testset "array equations over a 2D slice keep their shape" begin
    # A derivative of a 2D slice must substitute a 2D array of scalar derivatives; a
    # flattened one does not broadcast against the surrounding slices and codegen fails
    # with a DimensionMismatch.
    n = 6
    @independent_variables t
    @variables w(t)[1:n, 1:n]
    D = Differential(t)
    dx = 1 / (n - 1)
    inner = 2:(n - 1)
    lap = (
        w[1:(n - 2), inner] .+ w[3:n, inner] .+ w[inner, 1:(n - 2)] .+
            w[inner, 3:n] .- 4 .* w[inner, inner]
    ) ./ dx^2
    eqs = Equation[broadcast(-, D(w[inner, inner]), lap) ~ zeros(n - 2, n - 2)]
    for i in 1:n
        push!(eqs, w[i, 1] ~ 0.0)
        push!(eqs, w[i, n] ~ 0.0)
    end
    for j in inner
        push!(eqs, w[1, j] ~ 0.0)
        push!(eqs, w[n, j] ~ 0.0)
    end
    @named sys2d = System(eqs, t, vec(collect(w)), [])
    sys2d = complete(sys2d)

    op = vcat(
        [w[i, j] => 0.25 for i in 1:n, j in 1:n] |> vec,
        [D(w[i, j]) => 0.0 for i in 1:n, j in 1:n] |> vec
    )
    prob = DAEProblem(sys2d, op, (0.0, 0.01); build_initializeprob = false)
    @test length(prob.u0) == n * n
    out = zeros(n * n)
    prob.f(out, zeros(n * n), prob.u0, prob.p, 0.0)
    @test all(isfinite, out)
end

@testset "array equations written as `D(slice) ~ rhs`" begin
    # The residual form above puts the derivative inside the expression. The equivalent
    # `D(u[2:n-1]) ~ rhs` form has no scalar `toterm` name for its LHS, which the
    # derivative-substitution machinery has to skip rather than trip over.
    n = 11
    @independent_variables t
    @variables u(t)[1:n]
    D = Differential(t)
    dx = 1 / (n - 1)
    lap = (u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]) ./ dx^2
    eqs = [D(u[2:(n - 1)]) ~ lap, u[1] ~ 0.0, u[n] ~ 0.0]
    @named sys = System(eqs, t, collect(u), [])
    sys = complete(sys)

    xs = range(0.0, 1.0, length = n)
    op = vcat([u[i] => sinpi(xs[i]) for i in 1:n], [D(u[i]) => 0.0 for i in 1:n])
    prob = DAEProblem(sys, op, (0.0, 0.1); build_initializeprob = false)
    @test length(prob.u0) == n

    # the residual matches the analytic derivative of the initial condition
    out = zeros(n)
    du = zeros(n)
    du[2:(n - 1)] .= [-pi^2 * sinpi(x) for x in xs[2:(n - 1)]]
    prob.f(out, du, prob.u0, prob.p, 0.0)
    @test maximum(abs, out) < 1.0e-1

    sol = solve(
        prob, DFBDF(); initializealg = BrownFullBasicInit(), reltol = 1.0e-8,
        abstol = 1.0e-8, saveat = [0.1]
    )
    @test SciMLBase.successful_retcode(sol)
    @test maximum(abs, sol.u[end] .- [exp(-pi^2 * 0.1) * sinpi(x) for x in xs]) < 1.0e-2
end

@testset "has_array_equations detects every array-equation form" begin
    n = 5
    @independent_variables t
    @variables u(t)[1:n]
    D = Differential(t)
    lap = u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]
    @test has_array_equations([zeros(n - 2) ~ broadcast(-, lap)])
    @test has_array_equations([broadcast(-, lap) ~ zeros(n - 2)])
    @test has_array_equations([D(u[2:(n - 1)]) ~ lap])
    @test !has_array_equations([u[1] ~ 0.0, u[n] ~ 0.0])
    @test !has_array_equations(Equation[])
end

@testset "accepting array equations is a per-constructor capability" begin
    @test accepts_array_equations(DAEFunction)
    @test accepts_array_equations(NonlinearFunction)
    @test accepts_array_equations(ODEFunction)
    @test !accepts_array_equations(SDEFunction)
    @test !accepts_array_equations(ImplicitDiscreteFunction)
    # Optimization vectorizes `costs` and `constraints`, not `equations`:
    # `check_no_equations` rejects them before this gate is reached.
    @test !accepts_array_equations(OptimizationFunction)
    @test !accepts_array_equations(MultiObjectiveOptimizationFunction)
end

@testset "symbolic jacobian from array residuals" begin
    n = 11
    sys, u, t, D = heat_array_system(n)
    xs = range(0.0, 1.0, length = n)
    op = vcat(
        [u[i] => sinpi(xs[i]) for i in 1:n],
        [D(u[i]) => 0.0 for i in 1:n]
    )
    # the jacobian is built from the scalarized `full_equations`, one row per residual;
    # `sparse = true` is not tested because `W_sparsity` requires a semi-explicit mass
    # matrix, which no residual-form DAE has
    prob = DAEProblem(sys, op, (0.0, 0.1); build_initializeprob = false, jac = true)
    J = zeros(n, n)
    γ = 2.0
    prob.f.jac(J, prob.du0, prob.u0, prob.p, γ, 0.0)
    dx = 1 / (n - 1)
    # residual row 1 is `lap[1] - D(u[2])`
    @test J[1, 1:3] ≈ [1, -2 - γ * dx^2, 1] ./ dx^2
    # residual row `n - 1` is `0 - u[1]`
    @test J[n - 1, 1] ≈ -1
    @test count(!iszero, J) == 3 * (n - 2) + 2
    sol = solve(
        prob, DFBDF(); initializealg = BrownFullBasicInit(), reltol = 1.0e-8,
        abstol = 1.0e-8, saveat = [0.1]
    )
    @test SciMLBase.successful_retcode(sol)
    @test maximum(abs, sol.u[end] .- [exp(-pi^2 * 0.1) * sinpi(x) for x in xs]) < 1.0e-2
end
