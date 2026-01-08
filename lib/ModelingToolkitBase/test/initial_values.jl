using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D, get_u0
using OrdinaryDiffEq
using DataInterpolations
using StaticArrays
using SymbolicIndexingInterface
using Test

missing_guess_value = if @isdefined(ModelingToolkit)
    MissingGuessValue.Error()
else
    MissingGuessValue.Constant(1.0)
end

@variables x(t)[1:3] = [1.0, 2.0, 3.0] y(t) z(t)[1:2]

@mtkcompile sys = System([D(x) ~ t * x], t) simplify = false
reorderer = getsym(sys, x)
@test reorderer(get_u0(sys, [])) == [1.0, 2.0, 3.0]
@test reorderer(get_u0(sys, [x => [2.0, 3.0, 4.0]])) == [2.0, 3.0, 4.0]
@test reorderer(get_u0(sys, [x[1] => 2.0, x[2] => 3.0, x[3] => 4.0])) == [2.0, 3.0, 4.0]

@mtkcompile sys = System(
    [
        D(x) ~ 3x,
        D(y) ~ t,
        D(z[1]) ~ z[2] + t,
        D(z[2]) ~ y + z[1],
    ], t
) simplify = false

@test_throws ModelingToolkitBase.MissingVariablesError get_u0(sys, [])
getter = getu(sys, [x..., y, z...])
@test getter(get_u0(sys, [y => 4.0, z => [5.0, 6.0]])) == collect(1.0:6.0)
@test getter(get_u0(sys, [y => 4.0, z => [3y, 4y]])) == [1.0, 2.0, 3.0, 4.0, 12.0, 16.0]
@test getter(get_u0(sys, [y => 3.0, z[1] => 3y, z[2] => 2x[1]])) ==
    [1.0, 2.0, 3.0, 3.0, 9.0, 2.0]

@variables w(t)
@parameters p1 p2

@test getter(get_u0(sys, [y => 2p1, z => [3y, 2p2], p1 => 5.0, p2 => 6.0])) ==
    [1.0, 2.0, 3.0, 10.0, 30.0, 12.0]
@test_throws Any getter(get_u0(sys, [y => 2w, w => 3.0, z[1] => 2p1, z[2] => 3p2]))
@test getter(
    get_u0(
        sys, [y => 2w, w => 3.0, z[1] => 2p1, z[2] => 3p2, p1 => 3.0, p2 => 4.0]
    )
) ==
    [1.0, 2.0, 3.0, 6.0, 6.0, 12.0]

# Issue#2566
@variables X(t)
@parameters p1 p2 p3

p_vals = [p1 => 1.0, p2 => 2.0]
u_vals = [X => 3.0]

var_vals = [p1 => 1.0, p2 => 2.0, X => 3.0]
desired_values = [p1, p2, p3]
defaults = Dict([p3 => X])
vals = ModelingToolkitBase.varmap_to_vars(merge(defaults, Dict(var_vals)), desired_values)
@test vals == [1.0, 2.0, 3.0]

# Issue#2565
# Create ODESystem.
@variables X1(t) X2(t)
@parameters k1 k2 Γ[1:1] = missing
eq = D(X1) ~ -k1 * X1 + k2 * (-X1 + Γ[1])
obs = X2 ~ Γ[1] - X1
@mtkcompile osys_m = System([eq, obs], t, [X1, X2], [k1, k2, Γ])

# Creates ODEProblem.
u0 = [X1 => 1.0, X2 => 2.0]
tspan = (0.0, 1.0)
ps = [k1 => 1.0, k2 => 5.0]
oprob = ODEProblem(osys_m, [u0; ps], tspan; guesses = [Γ[1] => 1.0])
@test SciMLBase.successful_retcode(solve(oprob))

# Initialization of ODEProblem with dummy derivatives of multidimensional arrays
# Issue#1283
@variables z(t)[1:2, 1:2]
eqs = [D(D(z)) ~ ones(2, 2)]
@mtkcompile sys = System(eqs, t)
@test_nowarn ODEProblem(sys, [z => zeros(2, 2), D(z) => ones(2, 2)], (0.0, 10.0))

@testset "binding=nothing is skipped" begin
    @parameters p = nothing
    @variables x(t) = nothing y(t)
    @named sys = System(Equation[], t, [x, y], [p]; bindings = [y => nothing])
    @test isempty(ModelingToolkitBase.bindings(sys))
end

# Using indepvar in initialization
# Issue#2799
@variables x(t)
@parameters p
@mtkcompile sys = System([D(x) ~ p], t; bindings = [x => t, p => 2t])
prob = ODEProblem(sys, [], (1.0, 2.0))
@test prob[x] == 1.0
@test prob.ps[p] == 2.0

@testset "Array of symbolics is unwrapped" begin
    @variables x(t)[1:2] y(t)
    @mtkcompile sys = System([D(x) ~ x, D(y) ~ t], t; initial_conditions = [x => [y, 3.0]])
    prob = ODEProblem(sys, [y => 1.0], (0.0, 1.0))
    @test eltype(prob.u0) <: Float64
    prob = ODEProblem(sys, [x => [y, 4.0], y => 2.0], (0.0, 1.0))
    @test eltype(prob.u0) <: Float64
end

@testset "split=false systems with all parameter initial conditions" begin
    @variables x(t) = 1.0
    @parameters p = 1.0 q = 2.0 r = 3.0
    @mtkcompile sys = System(D(x) ~ p * x + q * t + r, t) split = false
    prob = @test_nowarn ODEProblem(sys, [], (0.0, 1.0))
    @test prob.p isa Vector{Float64}
end

@testset "Issue#3153" begin
    @variables x(t) y(t)
    @parameters c1 c2 c3
    eqs = [
        D(x) ~ y,
        y ~ ifelse(t < c1, 0.0, (-c1 + t)^(c3)),
    ]
    sps = [x, y]
    ps = [c1, c2, c3]
    @mtkcompile osys = System(eqs, t, sps, ps)
    u0map = [x => 1.0]
    pmap = [c1 => 5.0, c2 => 1.0, c3 => 1.2]
    oprob = ODEProblem(osys, [u0map; pmap], (0.0, 10.0))
end

@testset "Cyclic dependency checking and substitution limits" begin
    @variables x(t) y(t)
    @mtkcompile sys = System(
        [D(x) ~ x, D(y) ~ y], t; initialization_eqs = [x ~ 2y + 3, y ~ 2x],
        guesses = [x => 2y, y => 2x]
    )
    @test_warn ["Cycle", "unknowns", "x", "y"] ODEProblem(sys, [], (0.0, 1.0), warn_cyclic_dependency = true)
    @test_throws ModelingToolkitBase.UnexpectedSymbolicValueInVarmap ODEProblem(
        sys, [x => 2y + 1, y => 2x], (0.0, 1.0); build_initializeprob = false,
        substitution_limit = 10
    )

    @parameters p q
    @mtkcompile sys = System(
        [D(x) ~ x * p, D(y) ~ y * q], t; guesses = [p => 1.0, q => 2.0]
    )
    # "unknowns" because they are initialization unknowns
    @test_warn ["Cycle", "parameters", "p", "q"] ODEProblem(
        sys, [x => 1, y => 2, p => 2q, q => 3p],
        (0.0, 1.0); warn_cyclic_dependency = true
    )
end

@testset "`add_fallbacks!` checks scalarized array parameters correctly" begin
    @variables x(t)[1:2]
    @parameters p[1:2, 1:2]
    @mtkcompile sys = System(D(x) ~ p * x, t)
    # used to throw a `MethodError` complaining about `getindex(::Nothing, ::CartesianIndex{2})`
    @test_throws "Could not evaluate" ODEProblem(
        sys, [x => ones(2)], (0.0, 1.0)
    )
end

@testset "Unscalarized default for scalarized observed variable" begin
    @parameters p[1:4] = rand(4)
    @variables x(t)[1:4] y(t)[1:2]
    eqs = [
        D(x) ~ x,
        y[1] ~ x[3],
        y[2] ~ x[4],
    ]
    @mtkcompile sys = System(eqs, t; bindings = [x => vcat(ones(2), y), y => x[1:2] ./ 2])
    prob = ODEProblem(sys, [], (0.0, 1.0))
    sol = solve(prob)
    @test SciMLBase.successful_retcode(sol)
    @test sol[x, 1] ≈ [1.0, 1.0, 0.5, 0.5]
end

# Manually do index reduction to allow testing with MTKBase
function index_reduced_pend(; name)
    @parameters g
    @variables x(t) y(t) [state_priority = 10] λ(t) yˍt(t) xˍt(t) xˍtt(t)
    eqs = [
        0 ~ 1 - y^2 - x^2
        D(y) ~ yˍt
        0 ~ -yˍt * y - xˍt * x
        D(yˍt) ~ -g + y * λ
        0 ~ -yˍt^2 - xˍt^2 - x * xˍtt - y * (-g + y * λ)
    ]
    return System(eqs, t, [xˍt, y, x, yˍt, λ], [g]; observed = [xˍtt ~ λ * x], name)
end

@testset "Missing/cyclic guesses throws error" begin
    @parameters g
    @variables x(t) y(t) [state_priority = 10] λ(t)
    @named pend = index_reduced_pend()
    pend = complete(pend)

    @test_throws ModelingToolkitBase.MissingGuessError ODEProblem(
        pend, [x => 1, g => 1], (0, 1); guesses = [y => λ, λ => y + 1]
    )
    ODEProblem(
        pend, [x => 1, g => 1], (0, 1);
        guesses = [y => λ, λ => 0.5], missing_guess_value
    )

    # Throw multiple if multiple are missing
    @variables a(t) b(t) c(t) d(t) e(t)
    eqs = [D(a) ~ b, D(b) ~ c, D(c) ~ d, D(d) ~ e, D(e) ~ 1]
    @mtkcompile sys = System(eqs, t)
    @test_throws ["d(t)", "c(t)"] ODEProblem(
        sys, [e => 2, a => b, b => a + 1, c => d, d => c + 1], (0, 1)
    )
end

@testset "Issue#3490: `remake` works with callable parameters" begin
    ts = collect(0.0:0.1:10.0)
    spline = LinearInterpolation(ts .^ 2, ts)
    Tspline = typeof(spline)
    @variables x(t)
    @parameters (interp::Tspline)(..)

    @mtkcompile sys = System(D(x) ~ interp(t), t)

    prob = ODEProblem(sys, [x => 0.0, interp => spline], (0.0, 1.0))
    spline2 = LinearInterpolation(ts .^ 2, ts .^ 2)
    p_new = [interp => spline2]
    prob2 = remake(prob; p = p_new)
    @test prob2.ps[interp] == spline2
end

@testset "Issue#3523: don't substitute inside initial in `build_operating_point!`" begin
    @variables (X(t))[1:2]
    @parameters p[1:2]
    eqs = [
        0 ~ p[1] - X[1],
        0 ~ p[2] - X[2],
    ]
    @named nlsys = System(eqs)
    nlsys = complete(nlsys)

    # Creates the `NonlinearProblem`.
    u0 = [X => [1.0, 2.0]]
    ps = [p => [4.0, 5.0]]
    @test_nowarn NonlinearProblem(nlsys, [u0; ps])
end

@testset "Issue#3553: Retain `Float32` initial values" begin
    @parameters p d
    @variables X(t)
    eqs = [D(X) ~ p - d * X]
    @mtkcompile osys = System(eqs, t)
    u0 = [X => 1.0f0]
    ps = [p => 1.0f0, d => 2.0f0]
    oprob = ODEProblem(osys, [u0; ps], (0.0f0, 1.0f0))
    sol = solve(oprob)
    @test eltype(oprob.u0) == Float32
    @test eltype(eltype(sol.u)) == Float32
end

@testset "Array initials and scalar parameters with `split = false`" begin
    @variables x(t)[1:2]
    @parameters p
    @mtkcompile sys = System([D(x[1]) ~ x[1], D(x[2]) ~ x[2] + p], t) split = false
    ps = Set(parameters(sys; initial_parameters = true))
    @test length(ps) == 5
    for i in 1:2
        @test Initial(x[i]) in ps
        @test Initial(D(x[i])) in ps
    end
    @test p in ps
    prob = ODEProblem(sys, [x => ones(2), p => 1.0], (0.0, 1.0))
    @test prob.p isa Vector{Float64}
    @test length(prob.p) == 5
end

@testset "Temporary values for solved variables are guesses" begin
    @parameters σ ρ β = missing [guess = 8 / 3]
    @variables x(t) y(t) z(t) w(t) w2(t)

    eqs = [
        D(D(x)) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z,
        w ~ x + y + z + 2 * β,
        0 ~ x^2 + y^2 - w2^2,
    ]

    @mtkcompile sys = System(eqs, t)

    u0 = [
        D(x) => 2.0,
        x => 1.0,
        y => 0.0,
        z => 0.0,
    ]

    p = [
        σ => 28.0,
        ρ => 10.0,
    ]

    tspan = (0.0, 100.0)
    prob = ODEProblem(
        sys, [u0; p], tspan, jac = true, guesses = [w2 => -1.0],
        warn_initialize_determined = false
    )
    @test prob[w2] ≈ -1.0
    @test prob.ps[β] ≈ 8 / 3
end

@testset "MTKParameters uses given `pType` for inner buffers" begin
    @parameters σ ρ β
    @variables x(t) y(t) z(t)

    eqs = [
        D(D(x)) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z,
    ]

    @mtkcompile sys = System(eqs, t)

    u0 = SA[
        D(x) => 2.0f0,
        x => 1.0f0,
        y => 0.0f0,
        z => 0.0f0,
    ]

    p = SA[
        σ => 28.0f0,
        ρ => 10.0f0,
        β => 8.0f0 / 3,
    ]

    tspan = (0.0f0, 100.0f0)
    prob = ODEProblem(sys, [u0; p], tspan)
    @test prob.p.tunable isa SVector
    @test prob.p.initials isa SVector
end

@testset "`p_constructor` keyword argument" begin
    @parameters g = 1.0
    @variables x(t) y(t) [state_priority = 10, guess = 1.0] λ(t) [guess = 1.0]
    @named pend = index_reduced_pend()
    pend = complete(pend)
    guesses(pend)[y] = 1.0
    guesses(pend)[λ] = 1.0
    initial_conditions(pend)[g] = 1.0

    u0 = [y => -1.0, D(y) => 0.0]
    u0_constructor = p_constructor = vals -> SVector{length(vals)}(vals...)
    tspan = (0.0, 5.0)
    prob = ODEProblem(pend, u0, tspan; u0_constructor, p_constructor, missing_guess_value)
    @test prob.u0 isa SVector
    @test prob.p.tunable isa SVector
    @test prob.p.initials isa SVector
    initdata = prob.f.initialization_data
    @test state_values(initdata.initializeprob) isa SVector
    @test parameter_values(initdata.initializeprob).tunable isa SVector

    pend = complete(pend; split = false)
    prob = ODEProblem(pend, u0, tspan; u0_constructor, p_constructor, missing_guess_value)
    @test prob.p isa SVector
    initdata = prob.f.initialization_data
    @test state_values(initdata.initializeprob) isa SVector
    @test parameter_values(initdata.initializeprob) isa SVector
end

@testset "Type promotion of `p` works with non-dual types" begin
    @variables x(t) y(t)
    @mtkcompile sys = System([D(x) ~ x + y, x^3 + y^3 ~ 5], t; guesses = [y => 1.0])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0))
    prob2 = remake(prob; u0 = BigFloat.(prob.u0))
    @test prob2.p.initials isa Vector{BigFloat}
    sol = solve(prob2)
    @test SciMLBase.successful_retcode(sol)
end

@testset "metadata defaults for variables removed by `mtkcompile` are removed" begin
    @parameters p d
    @variables X(t) Y(t) = 1.0
    eq = [
        D(X) ~ p - d * X,
    ]
    @named sys = System(eq, t, [X, Y], [p, d])
    @test haskey(initial_conditions(sys), Y)
    sys = mtkcompile(sys)
    @test !haskey(initial_conditions(sys), Y)
    @test_nowarn ODEProblem(sys, [X => 0.0, p => 2.0, d => 0.1], (0.0, 1.0))

    @variables X(t) Y(t) = X
    @named sys = System(eq, t, [X, Y], [p, d])
    @test haskey(bindings(sys), Y)
    sys = mtkcompile(sys)
    @test !haskey(bindings(sys), Y)
    @test_nowarn ODEProblem(sys, [X => 0.0, p => 2.0, d => 0.1], (0.0, 1.0))
end
