using ModelingToolkit, Test
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkitBase as MTKBase

differential_states(sys) = [
    only(arguments(eq.lhs)) for eq in equations(sys) if ModelingToolkit.isdiffeq(eq)
]

function finite_at_t0(prob)
    du = similar(prob.u0)
    prob.f(du, prob.u0, prob.p, prob.tspan[1])
    return all(isfinite, prob.u0) && all(isfinite, du)
end

# Planar pendulum in Cartesian coordinates: one constraint `x² + y² = L²`, so exactly one
# of `x`, `y` stays a differential state and the other is solved from the constraint.
# Solving for a coordinate whose column of the constraint jacobian, `2x` or `2y`, vanishes
# at the initial point makes the reduced system singular there. The initial configuration
# is passed to `mtkcompile` through `initial_point`.
function pendulum(; x0, y0, vx0 = 0.0, initial_point = true, kwargs...)
    @parameters g = 9.81 L = 1.0
    @variables x(t) y(t) vx(t) vy(t) λ(t) [guess = 0.0]
    eqs = [D(x) ~ vx, D(y) ~ vy, D(vx) ~ λ * x, D(vy) ~ λ * y - g, 0 ~ x^2 + y^2 - L^2]
    point = initial_point ? [x => x0, y => y0, vx => vx0, vy => 0.0] : nothing
    sys = mtkcompile(System(eqs, t; name = :pendulum, kwargs...); initial_point = point)
    return sys, (; x, y, vx, vy, λ)
end

@testset "single vanishing constraint-jacobian column (pendulum)" begin
    # Swinging through the bottom: `∂g/∂x = 2x = 0` there, so `x` must remain a state.
    # The purely structural selection keeps `y`, which cannot represent this motion.
    sys, v = pendulum(; x0 = 0.0, y0 = -1.0, vx0 = 1.0)
    @test issetequal(differential_states(sys), [v.x, v.vx])
    prob = ODEProblem(sys, [v.x => 0.0, v.vx => 1.0, v.y => -1.0, v.vy => 0.0], (0.0, 0.3))
    @test finite_at_t0(prob)
    sol = solve(prob, Rodas5P(); abstol = 1.0e-8, reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
    @test sol[v.x^2 + v.y^2][end] ≈ 1.0 atol = 1.0e-6
    @test sol[v.x][end] > 0.2

    # Horizontal: `∂g/∂y = 2y = 0`, so `y` must remain a state.
    sys, v = pendulum(; x0 = -1.0, y0 = 0.0)
    @test issetequal(differential_states(sys), [v.y, v.vy])
    # (stop before the pendulum reaches the bottom, where `y` as a state is singular)
    prob = ODEProblem(sys, [v.y => 0.0, v.vy => 0.0, v.x => -1.0, v.vx => 0.0], (0.0, 0.3))
    @test finite_at_t0(prob)
    sol = solve(prob, Rodas5P(); abstol = 1.0e-8, reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
end

@testset "the structural selection is kept when it is not singular" begin
    # No initial point: the purely structural choice (`y`) stays.
    sys, v = pendulum(; x0 = 0.0, y0 = -1.0, vx0 = 1.0, initial_point = false)
    @test issetequal(differential_states(sys), [v.y, v.vy])
    # Guesses are not the initial point, even when they describe a singular configuration.
    @parameters g = 9.81 L = 1.0
    @variables x(t) [guess = 0.0] y(t) [guess = -1.0] vx(t) vy(t) λ(t)
    eqs = [D(x) ~ vx, D(y) ~ vy, D(vx) ~ λ * x, D(vy) ~ λ * y - g, 0 ~ x^2 + y^2 - L^2]
    @mtkcompile guessed = System(eqs, t)
    @test issetequal(differential_states(guessed), [y, vy])
    # Generic initial point: both choices are non-singular, the structural choice stays.
    sys, v = pendulum(; x0 = 0.6, y0 = -0.8)
    @test issetequal(differential_states(sys), [v.y, v.vy])
    # `state_priority` still decides between non-singular choices ...
    sys, v = pendulum(; x0 = 0.6, y0 = -0.8, state_priorities = [x => 10, vx => 10])
    @test issetequal(differential_states(sys), [v.x, v.vx])
    # ... but a prioritized choice that is singular at the initial point is not taken.
    sys, v = pendulum(; x0 = 0.0, y0 = -1.0, vx0 = 1.0, state_priorities = [y => 10, vy => 10])
    @test issetequal(differential_states(sys), [v.x, v.vx])
end

# Two constraints on three angles with one degree of freedom. At `a = b = 0` the columns of
# `a` and `b` in the constraint jacobian are parallel, so no single column vanishes, but
# solving the constraints for `{a, b}` is singular there; `{a, c}` or `{b, c}` are fine.
@testset "rank-deficient block without a vanishing column" begin
    @variables a(t) b(t) c(t) va(t) vb(t) vc(t) λ1(t) [guess = 0.0] λ2(t) [guess = 0.0]
    g = [sin(a) - sin(b), cos(a) + cos(b) - c]
    G = Symbolics.jacobian(g, [a, b, c])
    eqs = [
        D(a) ~ va, D(b) ~ vb, D(c) ~ vc,
        D(va) ~ 1 - (G' * [λ1, λ2])[1],
        D(vb) ~ -(G' * [λ1, λ2])[2],
        D(vc) ~ -1 - (G' * [λ1, λ2])[3],
        0 ~ g[1], 0 ~ g[2],
    ]
    point = [a => 0.0, b => 0.0, c => 2.0, va => 0.0, vb => 0.0, vc => 0.0]
    @named sys = System(eqs, t)
    # structurally, `c` is kept as the state
    @test any(isequal(c), differential_states(mtkcompile(sys)))
    sys = mtkcompile(sys; initial_point = point)
    states = differential_states(sys)
    @test length(states) == 2
    @test !any(isequal(c), states)
    prob = ODEProblem(sys, point, (0.0, 0.5))
    @test finite_at_t0(prob)
    sol = solve(prob, Rodas5P(); abstol = 1.0e-8, reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
end

# Car axis problem from the IVP test set (index 3). The coefficient `yb = r sin(ω t)` of the
# first constraint vanishes at `t = 0`. Tearing used to eliminate `λ1` through the `dyl`
# equation, dividing by `yb`, which made `xlˍtt(0)` a `NaN`.
function caraxis(; kwargs...)
    M, ϵ, L, L0, r, ω, g = 10.0, 1.0e-2, 1.0, 0.5, 0.1, 10.0, 1.0
    k = M * ϵ^2 / 2
    @variables xl(t) [guess = 0.0] yl(t) [guess = 0.5] xr(t) [guess = 1.0] yr(t) [guess = 0.5]
    @variables dxl(t) [guess = -0.5] dyl(t) [guess = 0.0] dxr(t) [guess = -0.5] dyr(t) [guess = 0.0]
    @variables λ1(t) [guess = 0.0] λ2(t) [guess = 0.0]
    yb = r * sin(ω * t)
    xb = sqrt(L^2 - yb^2)
    Ll = sqrt(xl^2 + yl^2)
    Lr = sqrt((xr - xb)^2 + (yr - yb)^2)
    eqs = [
        D(xl) ~ dxl, D(yl) ~ dyl, D(xr) ~ dxr, D(yr) ~ dyr,
        k * D(dxl) ~ (L0 - Ll) * xl / Ll + λ1 * xb + 2λ2 * (xl - xr),
        k * D(dyl) ~ (L0 - Ll) * yl / Ll + λ1 * yb + 2λ2 * (yl - yr) - k * g,
        k * D(dxr) ~ (L0 - Lr) * (xr - xb) / Lr - 2λ2 * (xl - xr),
        k * D(dyr) ~ (L0 - Lr) * (yr - yb) / Lr - 2λ2 * (yl - yr) - k * g,
        0 ~ xb * xl + yb * yl,
        0 ~ (xl - xr)^2 + (yl - yr)^2 - L^2,
    ]
    sys = mtkcompile(System(eqs, t; name = :caraxis, kwargs...))
    return sys, (; xl, yl, xr, yr, dxl, dyl, dxr, dyr, λ1, λ2, yb)
end

@testset "coefficient vanishing at t0 is not divided by (car axis)" begin
    sys, v = caraxis(; tspan = (0.0, 3.0))
    @test any(isequal(v.λ1), unknowns(sys))
    prob = ODEProblem(sys, [v.yl => 0.5, v.yr => 0.5, v.dyl => 0.0, v.dyr => 0.0], (0.0, 3.0))
    @test finite_at_t0(prob)
    sol = solve(prob, Rodas5P(); abstol = 1.0e-10, reltol = 1.0e-10)
    @test SciMLBase.successful_retcode(sol)
    # positions at t = 3 from the IVP test set reference solution
    ref = [0.493455784275402809e-1, 0.496989460230171153, 0.104174252488542151e1, 0.373911027265361256]
    @test sol[[v.xl, v.yl, v.xr, v.yr]][end] ≈ ref rtol = 1.0e-6
    # this is the same selection the manual `maybe_zeros` workaround gives
    sys2, v2 = caraxis(; tspan = (0.0, 3.0), maybe_zeros = [v.yb])
    @test issetequal(string.(unknowns(sys)), string.(unknowns(sys2)))
    # when the simulation starts where `yb ≠ 0`, or the initial time is unknown, the
    # original selection is kept
    sys3, v3 = caraxis(; tspan = (1.0, 3.0))
    @test !any(isequal(v3.λ1), unknowns(sys3))
    sys4, v4 = caraxis()
    @test !any(isequal(v4.λ1), unknowns(sys4))
end

@testset "non-finite initial values warn" begin
    @variables x(t) y(t)
    op = MTKBase.SymmapT([x => NaN, y => 1.0])
    @test_logs (:warn, r"given for the following unknowns is not finite.*x\(t\)") MTKBase.warn_nonfinite_u0([x, y], [NaN, 1.0], op)
    @test_logs (:warn, r"likely singular at the initial point") MTKBase.warn_nonfinite_u0([x, y], [1.0, Inf], op)
    @test_logs MTKBase.warn_nonfinite_u0([x, y], [1.0, 2.0], op)
    @mtkcompile sys = System([D(x) ~ -x], t)
    @test_logs (:warn, r"given for the following unknowns is not finite.*x\(t\)") match_mode = :any ODEProblem(sys, [x => NaN], (0.0, 1.0))
    # the car axis without `tspan` has no initial time, so the singular selection remains
    # and initialization derives a `NaN`
    sys, v = caraxis()
    @test_logs (:warn, r"likely singular at the initial point") match_mode = :any ODEProblem(
        sys, [v.yl => 0.5, v.yr => 0.5, v.dyl => 0.0, v.dyr => 0.0], (0.0, 3.0)
    )
end
