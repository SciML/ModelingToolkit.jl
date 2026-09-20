using ModelingToolkit, Test, ADTypes, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: ContinuousClock, unwrap
using SciMLBase: PeriodicClock, SolverStepClock
using SymbolicIndexingInterface: is_parameter
import ControlSystemsBase as CS
using OrdinaryDiffEqDefault

# Convert a partition to a ControlSystemsBase state-space model.
to_ss(p) = p.Ts === nothing ? CS.ss(p.A, p.B, p.C, p.D) : CS.ss(p.A, p.B, p.C, p.D, p.Ts)

# Compare the frequency responses of two systems on a common time domain.
function freqresp_isapprox(G1, G2; w = exp10.(range(-2, 0.5, length = 20)), atol = 1.0e-8)
    return isapprox(CS.freqresp(G1, w), CS.freqresp(G2, w); atol)
end

const dt = 0.1

@testset "Continuous-only model" begin
    @variables x(t) u(t) y(t)
    @named sys = System([D(x) ~ -x + u, y ~ 2x], t)
    op = Dict(x => 0.0, u => 0.0)
    hl = linearize_hybrid(sys, [u], [y]; op)
    @test length(hl.partitions) == 1
    @test hl.continuous_index == 1
    @test isempty(hl.connections)
    p = hl.partitions[1]
    @test p.clock isa ContinuousClock
    @test p.Ts === nothing
    @test isequal(p.inputs, unwrap.([u]))
    @test isequal(p.outputs, unwrap.([y]))
    @test isequal(p.unknowns, unwrap.([x]))
    @test p.nu_user == 1
    @test p.ny_user == 1
    mats, _, _ = linearize(sys, [u], [y]; op)
    @test p.A == mats.A
    @test p.B == mats.B
    @test p.C == mats.C
    @test p.D == mats.D
    @test p.x0 == [0.0]
end

@testset "Single clock hybrid loop" begin
    @variables x(t) y(t) u(t) d(t) yd(t) ud(t) r(t)
    @parameters kp = 2.0
    eqs = [
        yd ~ Sample(dt)(y)
        ud ~ kp * (r - yd)
        r ~ 1.0
        u ~ Hold(ud)
        D(x) ~ -x + u + d
        y ~ x
    ]
    @named sys = System(eqs, t)
    hl = @test_nowarn linearize_hybrid(sys, [d], [x]; op = Dict(x => 0.0, d => 0.0))
    @test length(hl.partitions) == 2
    @test hl.continuous_index == 1
    cont, disc = hl.partitions
    @test cont.clock isa ContinuousClock
    @test disc.clock == PeriodicClock(dt, 0.0)
    @test disc.Ts == dt

    # Continuous partition: the held controller output enters as a boundary input, the
    # sampled plant output leaves as a boundary output.
    @test isequal(cont.inputs, unwrap.([d, Hold(ud)]))
    @test isequal(cont.outputs, unwrap.([x, y]))
    @test cont.nu_user == 1
    @test cont.ny_user == 1
    @test cont.A == [-1.0;;]
    @test cont.B == [1.0 1.0]
    @test cont.C == [1.0; 1.0;;]
    @test cont.D == zeros(2, 2)

    # Discrete partition: a static gain from the sampled output to the controller output.
    @test isequal(disc.inputs, unwrap.([Sample(dt)(y)]))
    @test isequal(disc.outputs, unwrap.([ud]))
    @test disc.nu_user == 0
    @test disc.ny_user == 0
    @test isempty(disc.unknowns)
    @test size(disc.A) == (0, 0)
    @test size(disc.B) == (0, 1)
    @test size(disc.C) == (1, 0)
    @test disc.D == [-2.0;;]

    @test length(hl.connections) == 2
    @test (; from = 2, output = 1, to = 1, input = 2) in hl.connections
    @test (; from = 1, output = 2, to = 2, input = 1) in hl.connections

    # Assemble the sampled-data loop with the advanced interface of `feedback` and compare
    # with the closed-loop recursion x(k+1) = e^{-dt} x(k) + (1 - e^{-dt}) (u(k) + d(k)),
    # u(k) = -kp x(k).
    Pd = CS.c2d(to_ss(cont), dt)
    Cd = to_ss(disc)
    G = CS.feedback(
        Pd, Cd; U1 = [2], Y1 = [2], U2 = [1], Y2 = [1], W1 = [1], Z1 = [1], pos_feedback = true
    )
    a = exp(-dt)
    Gref = CS.ss(a - 2 * (1 - a), 1 - a, 1, 0, dt)
    @test freqresp_isapprox(G, Gref)
end

@testset "Discrete-only model" begin
    k = ShiftIndex(Clock(dt))
    @variables xd(t) ud(t) yd(t)
    eqs = [xd(k) ~ 0.5xd(k - 1) + ud(k), yd ~ xd]
    @named sys = System(eqs, t)
    hl = linearize_hybrid(sys, [ud], [yd]; op = Dict(xd => 0.0, ud => 0.0))
    @test length(hl.partitions) == 1
    @test hl.continuous_index === nothing
    @test isempty(hl.connections)
    p = hl.partitions[1]
    @test p.Ts == dt
    @test length(p.unknowns) == 1
    @test isequal(ModelingToolkit.getunshifted(p.unknowns[1]), unwrap(xd))
    @test p.A == [0.5;;]
    @test p.B == [1.0;;]
    @test p.C == [0.5;;]
    @test p.D == [1.0;;]
    # x(k) = 0.5 x(k-1) + u(k), y(k) = x(k) has the transfer function z / (z - 0.5)
    @test freqresp_isapprox(to_ss(p), CS.tf([1.0, 0.0], [1.0, -0.5], dt))
end

@testset "Deeper history and array discrete state" begin
    k = ShiftIndex(Clock(dt))
    @testset "Second-order difference equation" begin
        @variables xd(t) ud(t)
        eqs = [xd(k) ~ 0.5xd(k - 1) + 0.2xd(k - 2) + ud(k)]
        @named sys = System(eqs, t)
        hl = linearize_hybrid(sys, [ud], [xd]; op = Dict(xd => 0.0, ud => 0.0))
        p = hl.partitions[1]
        @test length(p.unknowns) == 2
        @test freqresp_isapprox(to_ss(p), CS.tf([1.0, 0.0, 0.0], [1.0, -0.5, -0.2], dt))
    end
    @testset "Array-valued discrete state space" begin
        Adv = [0.5 0.1; 0.0 0.8]
        Bdv = [1.0, 0.5]
        Cdv = [1.0, 1.0]
        @variables (xa(t))[1:2] ua(t) ya(t)
        @parameters Ad[1:2, 1:2] = Adv Bd[1:2] = Bdv Cd[1:2] = Cdv
        eqs = [
            xa(k) ~ Ad * xa(k - 1) + Bd * ua(k)
            ya ~ Cd' * xa
        ]
        @named sys = System(eqs, t)
        hl = @test_nowarn linearize_hybrid(sys, [ua], [ya]; op = Dict(xa => zeros(2), ua => 0.0))
        p = hl.partitions[1]
        @test length(p.unknowns) == 2
        @test isequal(p.inputs, unwrap.([ua]))
        @test isequal(p.outputs, unwrap.([ya]))
        # y(k) = Cd' x(k) = Cd' Ad x(k-1) + Cd' Bd u(k)
        @test freqresp_isapprox(to_ss(p), CS.ss(Adv, Bdv, Cdv' * Adv, Cdv' * Bdv, dt))
        @test sort(eigvals(p.A)) ≈ sort(eigvals(Adv))
    end
end

@testset "Multi-rate model with linear partitions" begin
    dt1 = 0.1
    dt2 = 0.3
    k1 = ShiftIndex(Clock(dt1))
    k2 = ShiftIndex(Clock(dt2))
    @variables xc(t) d(t) w(t) z(t)
    eqs = [
        D(xc) ~ -xc + Hold(z) + d
        w(k1) ~ 0.5w(k1 - 1) + 0.5Sample(dt1)(xc)
        z(k2) ~ 0.9z(k2 - 1) - 0.1Sample(dt2)(Hold(w))
    ]
    @named sys = System(eqs, t)
    hl = @test_nowarn linearize_hybrid(sys, [d], [xc]; op = Dict(xc => 0.0, d => 0.0, w => 0.0, z => 0.0))
    @test length(hl.partitions) == 3
    @test hl.continuous_index == 1
    cont = hl.partitions[1]
    fast = hl.partitions[findfirst(p -> p.Ts == dt1, hl.partitions)]
    slow = hl.partitions[findfirst(p -> p.Ts == dt2, hl.partitions)]

    # The held fast signal `Hold(w)` only feeds the slow partition and is not an input of
    # the continuous partition.
    @test isequal(cont.inputs, unwrap.([d, Hold(z)]))
    @test isequal(cont.outputs, unwrap.([xc]))
    @test cont.A == [-1.0;;]
    @test cont.B == [1.0 1.0]
    @test cont.C == [1.0;;]

    @test isequal(fast.inputs, unwrap.([Sample(dt1)(xc)]))
    @test isequal(fast.outputs, unwrap.([w]))
    @test fast.A == [0.5;;]
    @test fast.B == [0.5;;]
    @test fast.C == [0.5;;]
    @test fast.D == [0.5;;]

    @test isequal(slow.inputs, unwrap.([Sample(dt2)(Hold(w))]))
    @test isequal(slow.outputs, unwrap.([z]))
    @test slow.A == [0.9;;]
    @test slow.B == [-0.1;;]
    @test slow.C == [0.9;;]
    @test slow.D == [-0.1;;]

    @test length(hl.connections) == 3
    ifast = findfirst(p -> p.Ts == dt1, hl.partitions)
    islow = findfirst(p -> p.Ts == dt2, hl.partitions)
    @test (; from = 1, output = 1, to = ifast, input = 1) in hl.connections
    @test (; from = ifast, output = 1, to = islow, input = 1) in hl.connections
    @test (; from = islow, output = 1, to = 1, input = 2) in hl.connections
end

@testset "Operating point" begin
    @variables x(t) y(t) u(t) yd(t) ud(t)
    @parameters kp = 2.0
    eqs = [
        yd ~ Sample(dt)(y)
        ud ~ -kp * yd
        u ~ Hold(ud)
        D(x) ~ -x^3 + u
        y ~ x
    ]
    @named sys = System(eqs, t)
    @testset "Dictionary" begin
        hl = @test_nowarn linearize_hybrid(sys, [], [y]; op = Dict(x => 2.0))
        cont = hl.partitions[hl.continuous_index]
        @test cont.A ≈ [-12.0;;]
        @test cont.x0 == [2.0]
    end
    @testset "Boundary values are resolved from the source partition" begin
        # The held input enters the plant through a nonlinearity, so the value of the
        # controller output at the operating point determines the input matrix.
        @variables xn(t) yn(t) un(t) ydn(t) udn(t)
        eqs_n = [
            ydn ~ Sample(dt)(yn)
            udn ~ -kp * ydn^2
            un ~ Hold(udn)
            D(xn) ~ -xn + un^2
            yn ~ xn
        ]
        @named sysn = System(eqs_n, t)
        hl = @test_nowarn linearize_hybrid(sysn, [], [yn]; op = Dict(xn => 2.0))
        cont = hl.partitions[hl.continuous_index]
        # u = -kp yn^2 = -8 and B = 2u
        @test cont.B ≈ [-16.0;;]
        disc = hl.partitions[2]
        # D = -2 kp yn = -8
        @test disc.D ≈ [-8.0;;]
    end
    @testset "Missing values default to zero with a warning" begin
        k = ShiftIndex(Clock(dt))
        @variables xm(t) ym(t) ydm(t) udm(t)
        eqs_m = [
            ydm ~ Sample(dt)(ym)
            udm(k) ~ 0.5udm(k - 1) - kp * ydm
            D(xm) ~ -xm^3 + Hold(udm)
            ym ~ xm
        ]
        @named sysm = System(eqs_m, t)
        # The history of the discrete state has no value.
        hl = @test_logs (:warn, r"No operating-point value") match_mode = :any linearize_hybrid(sysm, [], [ym]; op = Dict(xm => 1.0))
        cont = hl.partitions[hl.continuous_index]
        @test cont.A ≈ [-3.0;;]
        @test_nowarn linearize_hybrid(sysm, [], [ym]; op = Dict(xm => 1.0), warn_missing_op = false)
        @test_nowarn linearize_hybrid(sysm, [], [ym]; op = Dict(xm => 1.0, udm => 0.0))
    end
    @testset "Solution" begin
        @variables xs(t) us(t) ys(t)
        @named sim = System([D(xs) ~ -xs^3 + 0.5], t)
        prob = ODEProblem(mtkcompile(sim), [xs => 2.0], (0.0, 1.0))
        sol = solve(prob)
        @named lin = System([D(xs) ~ -xs^3 + us, ys ~ xs], t)
        hl = linearize_hybrid(lin, [us], [ys]; op = LinearizationOpPoint(sol, 0.7; op = Dict(us => 0.5)))
        p = hl.partitions[1]
        @test p.x0 ≈ [sol(0.7; idxs = xs)]
        @test p.A ≈ [-3 * sol(0.7; idxs = xs)^2;;]
    end
end

nodual(x::Float64) = 2x
@register_symbolic nodual(x)

@testset "Automatic differentiation backends" begin
    @variables x(t) y(t) d(t) yd(t) ud(t)
    @parameters kp = 2.0
    eqs = [
        yd ~ Sample(dt)(y)
        ud ~ -kp * yd
        D(x) ~ -x^2 + Hold(ud) + d
        y ~ x
    ]
    @named sys = System(eqs, t)
    op = Dict(x => 1.5, d => 0.0)
    hl_fd = linearize_hybrid(sys, [d], [y]; op)
    hl_fdiff = linearize_hybrid(sys, [d], [y]; op, autodiff = AutoFiniteDiff())
    for (p1, p2) in zip(hl_fd.partitions, hl_fdiff.partitions)
        @test p1.A ≈ p2.A atol = 1.0e-6
        @test p1.B ≈ p2.B atol = 1.0e-6
        @test p1.C ≈ p2.C atol = 1.0e-6
        @test p1.D ≈ p2.D atol = 1.0e-6
    end

    @testset "Fallback to finite differences" begin
        @variables yn(t)
        eqs_nd = [
            yd ~ Sample(dt)(y)
            ud ~ -kp * nodual(yd)
            D(x) ~ -x^2 + Hold(ud) + d
            y ~ x
        ]
        @named sysn = System(eqs_nd, t)
        hl = @test_logs (:warn, r"retrying with") match_mode = :any linearize_hybrid(sysn, [d], [y]; op)
        disc = hl.partitions[2]
        @test disc.D ≈ [-4.0;;] atol = 1.0e-6
        @test_throws Exception linearize_hybrid(sysn, [d], [y]; op, fallback_autodiff = nothing)
        # Errors unrelated to differentiation are not retried.
        @variables zq(t)
        @named sysq = System([eqs_nd; D(zq) ~ zq], t)
        @test_throws ModelingToolkit.IONotFoundError linearize_hybrid(sysq, [d], [zq^2]; op)
    end
end

@testset "Warnings for unsupported constructs" begin
    @variables x(t) y(t) yd(t) ud(t)
    @parameters kp = 2.0
    eqs = [
        yd ~ Sample(dt)(y)
        ud ~ -kp * yd
        D(x) ~ -x + Hold(ud)
        y ~ x
    ]
    @named sys = System(eqs, t; continuous_events = [[x ~ 0.5] => [x ~ 0.0]])
    op = Dict(x => 0.0)
    @test_logs (:warn, r"Events are not accounted for") match_mode = :any linearize_hybrid(sys, [], [y]; op)
    @test_nowarn linearize_hybrid(sys, [], [y]; op, warn_unsupported = false)
end

@testset "Model constructs that partitions must not inherit" begin
    @variables x(t) y(t) u(t) yd(t) ud(t)
    @parameters kp = 2.0
    k = ShiftIndex(Clock(dt))
    @testset "Assertions" begin
        # An assertion on a clocked variable inherited by the continuous partition
        eqs = [yd ~ Sample(dt)(y), ud ~ -kp * yd, u ~ Hold(ud), D(x) ~ -x + u, y ~ x]
        @named sys = System(eqs, t; assertions = Dict(ud < 10.0 => "ud out of range"))
        hl = @test_logs (:warn, r"Assertions are not accounted for") match_mode = :any linearize_hybrid(sys, [], [x]; op = Dict(x => 0.0))
        @test hl.partitions[1].A == [-1.0;;]
        # An assertion on a continuous variable inherited by a clocked partition with a state
        eqs = [yd ~ Sample(dt)(y), ud(k) ~ 0.5ud(k - 1) - kp * yd, u ~ Hold(ud), D(x) ~ -x + u, y ~ x]
        @named sys = System(eqs, t; assertions = Dict(x < 10.0 => "x out of range"))
        hl = linearize_hybrid(sys, [], [x]; op = Dict(x => 0.0, ud => 0.0), warn_unsupported = false)
        @test hl.partitions[2].A == [0.5;;]
        @test isempty(ModelingToolkit.assertions(hl.partitions[2].sys))
    end
    @testset "Parameters of other partitions" begin
        # `pm` is bound to `missing` and used by the plant only. The clocked partition must
        # neither require a value for it nor keep it as a parameter.
        @parameters pm
        eqs = [yd ~ Sample(dt)(y), ud(k) ~ 0.5ud(k - 1) - kp * yd, u ~ Hold(ud), D(x) ~ -pm * x + u, y ~ x]
        @named sys = System(eqs, t; bindings = [pm => missing])
        hl = @test_nowarn linearize_hybrid(sys, [], [x]; op = Dict(x => 0.0, ud => 0.0, pm => 3.0))
        cont, disc = hl.partitions
        @test cont.A == [-3.0;;]
        @test !is_parameter(disc.sys, pm)
        @test is_parameter(cont.sys, pm)
        @test !is_parameter(cont.sys, kp)
    end
    @testset "Initial value of an observed alias of a boundary signal" begin
        # `u` is an alias of the held controller output and carries an initial value, which
        # must not constrain the initialization of the continuous partition.
        @variables ui(t) = 0.3
        eqs = [yd ~ Sample(dt)(y), ud ~ -kp * yd, ui ~ Hold(ud), D(x) ~ -x + ui, y ~ x]
        @named sys = System(eqs, t)
        hl = @test_nowarn linearize_hybrid(sys, [], [x]; op = Dict(x => 1.0))
        cont = hl.partitions[hl.continuous_index]
        @test cont.A == [-1.0;;]
        @test cont.x0 == [1.0]
        # Likewise for an operating-point value given for such a variable.
        hl = @test_nowarn linearize_hybrid(sys, [], [x]; op = Dict(x => 1.0, ui => 0.3))
        @test hl.partitions[hl.continuous_index].x0 == [1.0]
    end
end

@testset "Analysis points" begin
    @variables x(t) y(t) yd(t) ud(t)
    @variables uh(t) [output = true] u(t) [input = true]
    @parameters kp = 2.0
    eqs = [
        yd ~ Sample(dt)(y)
        ud ~ -kp * yd
        uh ~ Hold(ud)
        connect(uh, :plant_input, u)
        D(x) ~ -x + u
        y ~ x
    ]
    @named sys = System(eqs, t)
    op = Dict(x => 0.0)

    @testset "linearize_hybrid with analysis points" begin
        hl = linearize_hybrid(sys, :plant_input, [y]; op)
        cont = hl.partitions[hl.continuous_index]
        @test cont.nu_user == 1
        @test cont.ny_user == 1
        @test length(cont.inputs) == 2
        @test cont.A == [-1.0;;]
        @test cont.B == [1.0 1.0]
        @test cont.C == [1.0;;]
    end
    @testset "get_sensitivity" begin
        hl = get_sensitivity(sys, :plant_input; op, hybrid = true)
        @test hl isa HybridLinearization
        cont = hl.partitions[hl.continuous_index]
        disc = hl.partitions[2]
        @test cont.nu_user == 1
        @test cont.ny_user == 1
        @test disc.D == [-2.0;;]
        # Sensitivity of the sampled-data loop, assembled in discrete time.
        Pd = CS.c2d(to_ss(cont), dt)
        Cd = to_ss(disc)
        S = CS.feedback(
            Pd, Cd; U1 = [2], Y1 = [2], U2 = [1], Y2 = [1], W1 = [1], Z1 = [1], pos_feedback = true
        )
        a = exp(-dt)
        # u = d + uh, uh = -kp y_sampled, y(k+1) = a y(k) + (1 - a) u(k)
        Sref = CS.ss(a - 2 * (1 - a), 1 - a, -2.0, 1.0, dt)
        @test freqresp_isapprox(S, Sref)
    end
    @testset "get_comp_sensitivity and get_looptransfer" begin
        hl = get_comp_sensitivity(sys, :plant_input; op, hybrid = true)
        @test hl isa HybridLinearization
        @test hl.partitions[hl.continuous_index].nu_user == 1
        # The input of the broken connection needs an operating-point value.
        hl = @test_nowarn get_looptransfer(sys, :plant_input; op = Dict(x => 0.0, u => 0.0), hybrid = true)
        @test hl isa HybridLinearization
        @test hl.partitions[hl.continuous_index].nu_user == 1
        # Without `hybrid`, the standard compiler rejects the hybrid model.
        @test_throws ModelingToolkit.HybridSystemNotSupportedException get_looptransfer(sys, :plant_input; op)
    end
    @testset "Analysis point on a discrete connection" begin
        @variables ydo(t) [output = true] ydi(t) [input = true]
        eqs_d = [
            ydo ~ Sample(dt)(y)
            connect(ydo, :ctrl_input, ydi)
            ud ~ -kp * ydi
            D(x) ~ -x + Hold(ud)
            y ~ x
        ]
        @named sysd = System(eqs_d, t)
        hl = get_sensitivity(sysd, :ctrl_input; op = Dict(x => 0.0), hybrid = true)
        disc = hl.partitions[2]
        @test disc.nu_user == 1
        @test disc.ny_user == 1
        @test disc.D[1, 1] == 1.0
        @test disc.D[1, 2] == 1.0
        @test isequal(disc.outputs[2], unwrap(ud))
        @test disc.D[2, 1] == -2.0
    end
end
