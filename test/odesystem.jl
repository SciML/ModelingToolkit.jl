using ModelingToolkit, StaticArrays, LinearAlgebra
using OrdinaryDiffEq, Sundials
using DiffEqBase, SparseArrays
using StaticArrays
using Test

using ModelingToolkit: value

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)

# Define a differential equation
eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

ModelingToolkit.toexpr.(eqs)[1]
@named de = ODESystem(eqs; defaults = Dict(x => 1))
@test eval(toexpr(de)) == de
@test hash(deepcopy(de)) == hash(de)

generate_function(de)

function test_diffeq_inference(name, sys, iv, dvs, ps)
    @testset "ODESystem construction: $name" begin
        @test isequal(independent_variables(sys)[1], value(iv))
        @test length(independent_variables(sys)) == 1
        @test isempty(setdiff(Set(states(sys)), Set(value.(dvs))))
        @test isempty(setdiff(Set(parameters(sys)), Set(value.(ps))))
    end
end

test_diffeq_inference("standard", de, t, [x, y, z], [ρ, σ, β])
generate_function(de, [x, y, z], [σ, ρ, β])
jac_expr = generate_jacobian(de)
jac = calculate_jacobian(de)
jacfun = eval(jac_expr[2])

for f in [
    ODEFunction(de, [x, y, z], [σ, ρ, β], tgrad = true, jac = true),
    eval(ODEFunctionExpr(de, [x, y, z], [σ, ρ, β], tgrad = true, jac = true)),
]
    # iip
    du = zeros(3)
    u = collect(1:3)
    p = collect(4:6)
    f.f(du, u, p, 0.1)
    @test du == [4, 0, -16]

    # oop
    du = @SArray zeros(3)
    u = SVector(1:3...)
    p = SVector(4:6...)
    @test f.f(u, p, 0.1) === @SArray [4, 0, -16]

    # iip vs oop
    du = zeros(3)
    g = similar(du)
    J = zeros(3, 3)
    u = collect(1:3)
    p = collect(4:6)
    f.f(du, u, p, 0.1)
    @test du == f(u, p, 0.1)
    f.tgrad(g, u, p, t)
    @test g == f.tgrad(u, p, t)
    f.jac(J, u, p, t)
    @test J == f.jac(u, p, t)
end

eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y * t,
    D(z) ~ x * y - β * z]
@named de = ODESystem(eqs)
ModelingToolkit.calculate_tgrad(de)

tgrad_oop, tgrad_iip = eval.(ModelingToolkit.generate_tgrad(de))

u = SVector(1:3...)
p = SVector(4:6...)
@test tgrad_oop(u, p, t) == [0.0, -u[2], 0.0]
du = zeros(3)
tgrad_iip(du, u, p, t)
@test du == [0.0, -u[2], 0.0]

@parameters σ′(t - 1)
eqs = [D(x) ~ σ′ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]
@named de = ODESystem(eqs)
test_diffeq_inference("global iv-varying", de, t, (x, y, z), (σ′, ρ, β))

f = eval(generate_function(de, [x, y, z], [σ′, ρ, β])[2])
du = [0.0, 0.0, 0.0]
f(du, [1.0, 2.0, 3.0], [x -> x + 7, 2, 3], 5.0)
@test du ≈ [11, -3, -7]

@parameters σ(..)
eqs = [D(x) ~ σ(t - 1) * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]
@named de = ODESystem(eqs)
test_diffeq_inference("single internal iv-varying", de, t, (x, y, z), (σ(t - 1), ρ, β))
f = eval(generate_function(de, [x, y, z], [σ, ρ, β])[2])
du = [0.0, 0.0, 0.0]
f(du, [1.0, 2.0, 3.0], [x -> x + 7, 2, 3], 5.0)
@test du ≈ [11, -3, -7]

eqs = [D(x) ~ x + 10σ(t - 1) + 100σ(t - 2) + 1000σ(t^2)]
@named de = ODESystem(eqs)
test_diffeq_inference("many internal iv-varying", de, t, (x,), (σ(t - 2), σ(t^2), σ(t - 1)))
f = eval(generate_function(de, [x], [σ])[2])
du = [0.0]
f(du, [1.0], [t -> t + 2], 5.0)
@test du ≈ [27561]

# Conversion to first-order ODEs #17
D3 = Differential(t)^3
D2 = Differential(t)^2
@variables u(t) uˍtt(t) uˍt(t) xˍt(t)
eqs = [D3(u) ~ 2(D2(u)) + D(u) + D(x) + 1
       D2(x) ~ D(x) + 2]
@named de = ODESystem(eqs)
de1 = ode_order_lowering(de)
lowered_eqs = [D(uˍtt) ~ 2uˍtt + uˍt + xˍt + 1
               D(xˍt) ~ xˍt + 2
               D(uˍt) ~ uˍtt
               D(u) ~ uˍt
               D(x) ~ xˍt]

#@test de1 == ODESystem(lowered_eqs)

# issue #219
@test all(isequal.([ModelingToolkit.var_from_nested_derivative(eq.lhs)[1]
                    for eq in equations(de1)],
                   states(@named lowered = ODESystem(lowered_eqs))))

test_diffeq_inference("first-order transform", de1, t, [uˍtt, xˍt, uˍt, u, x], [])
du = zeros(5)
ODEFunction(de1, [uˍtt, xˍt, uˍt, u, x], [])(du, ones(5), nothing, 0.1)
@test du == [5.0, 3.0, 1.0, 1.0, 1.0]

# Internal calculations
@parameters σ
a = y - x
eqs = [D(x) ~ σ * a,
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]
@named de = ODESystem(eqs)
generate_function(de, [x, y, z], [σ, ρ, β])
jac = calculate_jacobian(de)
@test ModelingToolkit.jacobian_sparsity(de).colptr == sparse(jac).colptr
@test ModelingToolkit.jacobian_sparsity(de).rowval == sparse(jac).rowval

f = ODEFunction(de, [x, y, z], [σ, ρ, β])

D = Differential(t)
@parameters A B C
_x = y / C
eqs = [D(x) ~ -A * x,
    D(y) ~ A * x - B * _x]
@named de = ODESystem(eqs)
@test begin
    local f, du
    f = eval(generate_function(de, [x, y], [A, B, C])[2])
    du = [0.0, 0.0]
    f(du, [1.0, 2.0], [1, 2, 3], 0.0)
    du ≈ [-1, -1 / 3]
    f = eval(generate_function(de, [x, y], [A, B, C])[1])
    du ≈ f([1.0, 2.0], [1, 2, 3], 0.0)
end

function lotka(u, p, t)
    x = u[1]
    y = u[2]
    [p[1] * x - p[2] * x * y,
        -p[3] * y + p[4] * x * y]
end

prob = ODEProblem(ODEFunction{false}(lotka), [1.0, 1.0], (0.0, 1.0), [1.5, 1.0, 3.0, 1.0])
de = modelingtoolkitize(prob)
ODEFunction(de)(similar(prob.u0), prob.u0, prob.p, 0.1)

function lotka(du, u, p, t)
    x = u[1]
    y = u[2]
    du[1] = p[1] * x - p[2] * x * y
    du[2] = -p[3] * y + p[4] * x * y
end

prob = ODEProblem(lotka, [1.0, 1.0], (0.0, 1.0), [1.5, 1.0, 3.0, 1.0])

de = modelingtoolkitize(prob)
ODEFunction(de)(similar(prob.u0), prob.u0, prob.p, 0.1)

# automatic state detection for DAEs
@parameters t k₁ k₂ k₃
@variables y₁(t) y₂(t) y₃(t)
D = Differential(t)
# reorder the system just to be a little spicier
eqs = [D(y₁) ~ -k₁ * y₁ + k₃ * y₂ * y₃,
    0 ~ y₁ + y₂ + y₃ - 1,
    D(y₂) ~ k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃]
@named sys = ODESystem(eqs, defaults = [k₁ => 100, k₂ => 3e7, y₁ => 1.0])
u0 = Pair[]
push!(u0, y₂ => 0.0)
push!(u0, y₃ => 0.0)
p = [k₁ => 0.04,
    k₃ => 1e4]
p2 = (k₁ => 0.04,
      k₂ => 3e7,
      k₃ => 1e4)
tspan = (0.0, 100000.0)
prob1 = ODEProblem(sys, u0, tspan, p)
prob12 = ODEProblem(sys, u0, tspan, [0.04, 3e7, 1e4])
prob13 = ODEProblem(sys, u0, tspan, (0.04, 3e7, 1e4))
prob14 = ODEProblem(sys, u0, tspan, p2)
for p in [prob1, prob14]
    @test Set(Num.(parameters(sys)) .=> p.p) == Set([k₁ => 0.04, k₂ => 3e7, k₃ => 1e4])
    @test Set(Num.(states(sys)) .=> p.u0) == Set([y₁ => 1, y₂ => 0, y₃ => 0])
end
prob2 = ODEProblem(sys, u0, tspan, p, jac = true)
prob3 = ODEProblem(sys, u0, tspan, p, jac = true, sparse = true)
@test prob3.f.jac_prototype isa SparseMatrixCSC
prob3 = ODEProblem(sys, u0, tspan, p, jac = true, sparsity = true)
@test prob3.f.sparsity isa SparseMatrixCSC
@test_throws ArgumentError ODEProblem(sys, zeros(5), tspan, p)
for (prob, atol) in [(prob1, 1e-12), (prob2, 1e-12), (prob3, 1e-12)]
    local sol
    sol = solve(prob, Rodas5())
    @test all(x -> ≈(sum(x), 1.0, atol = atol), sol.u)
end

du0 = [D(y₁) => -0.04
       D(y₂) => 0.04
       D(y₃) => 0.0]
prob4 = DAEProblem(sys, du0, u0, tspan, p2)
prob5 = eval(DAEProblemExpr(sys, du0, u0, tspan, p2))
for prob in [prob4, prob5]
    local sol
    @test prob.differential_vars == [true, true, false]
    sol = solve(prob, IDA())
    @test all(x -> ≈(sum(x), 1.0, atol = 1e-12), sol.u)
end

@parameters t σ β
@variables x(t) y(t) z(t)
D = Differential(t)
eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x - β * y,
    x + z ~ y]
@named sys = ODESystem(eqs)
@test all(isequal.(states(sys), [x, y, z]))
@test all(isequal.(parameters(sys), [σ, β]))
@test equations(sys) == eqs
@test ModelingToolkit.isautonomous(sys)

# issue 701
using ModelingToolkit
@parameters t a
@variables x(t)
D = Differential(t)
@named sys = ODESystem([D(x) ~ a])
@test equations(sys)[1].rhs isa Sym

# issue 708
@parameters t a
@variables x(t) y(t) z(t)
D = Differential(t)
@named sys = ODESystem([D(x) ~ y, 0 ~ x + z, 0 ~ x - y], t, [z, y, x], [])
sys2 = ode_order_lowering(sys)
M = ModelingToolkit.calculate_massmatrix(sys2)
@test M == Diagonal([1, 0, 0])

# issue #609
@variables t x1(t) x2(t)
D = Differential(t)

eqs = [
    D(x1) ~ -x1,
    0 ~ x1 - x2,
]
@named sys = ODESystem(eqs, t)
@test isequal(ModelingToolkit.get_iv(sys), t)
@test isequal(states(sys), [x1, x2])
@test isempty(parameters(sys))

# one equation ODESystem test
@parameters t r
@variables x(t)
D = Differential(t)
eq = D(x) ~ r * x
@named ode = ODESystem(eq)
@test equations(ode) == [eq]
# issue #808
@testset "Combined system name collisions" begin
    function makesys(name)
        @parameters t a
        @variables x(t) f(t)
        D = Differential(t)

        ODESystem([D(x) ~ -a * x + f]; name)
    end

    function issue808()
        sys1 = makesys(:sys1)
        sys2 = makesys(:sys1)

        @parameters t
        D = Differential(t)
        @test_throws ArgumentError ODESystem([sys2.f ~ sys1.x, D(sys1.f) ~ 0], t,
                                             systems = [sys1, sys2], name = :foo)
    end
    issue808()
end

#Issue 998
@parameters t
pars = []
vars = @variables((u1,))
der = Differential(t)
eqs = [
    der(u1) ~ 1,
]
@test_throws ArgumentError ODESystem(eqs, t, vars, pars, name = :foo)

#Issue 1063/998
pars = [t]
vars = @variables((u1(t),))
@test_throws ArgumentError ODESystem(eqs, t, vars, pars, name = :foo)

@parameters w
der = Differential(w)
eqs = [
    der(u1) ~ t,
]
@test_throws ArgumentError ModelingToolkit.ODESystem(eqs, t, vars, pars, name = :foo)

@variables x(t)
D = Differential(t)
@parameters M b k
eqs = [D(D(x)) ~ -b / M * D(x) - k / M * x]
ps = [M, b, k]
default_u0 = [D(x) => 0.0, x => 10.0]
default_p = [M => 1.0, b => 1.0, k => 1.0]
@named sys = ODESystem(eqs, t, [x], ps, defaults = [default_u0; default_p])
sys = ode_order_lowering(sys)
prob = ODEProblem(sys, [], tspan)
sol = solve(prob, Tsit5())
@test sum(abs, sol[end]) < 1

# check_eqs_u0 kwarg test
@parameters t
@variables x1(t) x2(t)
D = Differential(t)
eqs = [D(x1) ~ -x1]
@named sys = ODESystem(eqs, t, [x1, x2], [])
@test_throws ArgumentError ODEProblem(sys, [1.0, 1.0], (0.0, 1.0))
@test_nowarn ODEProblem(sys, [1.0, 1.0], (0.0, 1.0), check_length = false)

# check inputs
let
    @parameters t f k d
    @variables x(t) ẋ(t)
    δ = Differential(t)

    eqs = [δ(x) ~ ẋ, δ(ẋ) ~ f - k * x - d * ẋ]
    @named sys = ODESystem(eqs, t, [x, ẋ], [f, d, k]; controls = [f])

    calculate_control_jacobian(sys)

    @test isequal(calculate_control_jacobian(sys),
                  reshape(Num[0, 1], 2, 1))
end

# issue 1109
let
    @variables t x[1:3, 1:3](t)
    D = Differential(t)
    @named sys = ODESystem(D.(x) .~ x)
    @test_nowarn structural_simplify(sys)
end

# Array vars
using Symbolics: unwrap, wrap
using LinearAlgebra
@variables t
sts = @variables x[1:3](t)=[1, 2, 3.0] y(t)=1.0
ps = @parameters p[1:3] = [1, 2, 3]
D = Differential(t)
eqs = [collect(D.(x) .~ x)
       D(y) ~ norm(collect(x)) * y - x[1]]
@named sys = ODESystem(eqs, t, [sts...;], [ps...;])
sys = structural_simplify(sys)
@test isequal(@nonamespace(sys.x), x)
@test isequal(@nonamespace(sys.y), y)
@test isequal(@nonamespace(sys.p), p)
@test_nowarn sys.x, sys.y, sys.p
@test all(x -> x isa Symbolics.Arr, (sys.x, sys.p))
@test all(x -> x isa Symbolics.Arr, @nonamespace (sys.x, sys.p))
@test ModelingToolkit.isvariable(Symbolics.unwrap(x[1]))
prob = ODEProblem(sys, [], (0, 1.0))
sol = solve(prob, Tsit5())
@test sol[2x[1] + 3x[3] + norm(x)] ≈
      2sol[x[1]] + 3sol[x[3]] + vec(mapslices(norm, hcat(sol[x]...), dims = 2))
@test sol[x + [y, 2y, 3y]] ≈ sol[x] + [sol[y], 2sol[y], 3sol[y]]

# Mixed Difference Differential equations
@parameters t a b c d
@variables x(t) y(t)
δ = Differential(t)
Δ = Difference(t; dt = 0.1)
U = DiscreteUpdate(t; dt = 0.1)
eqs = [δ(x) ~ a * x - b * x * y
       δ(y) ~ -c * y + d * x * y
       Δ(x) ~ y
       U(y) ~ x + 1]
@named de = ODESystem(eqs, t, [x, y], [a, b, c, d])
@test generate_difference_cb(de) isa ModelingToolkit.DiffEqCallbacks.DiscreteCallback

# doesn't work with ODEFunction
# prob = ODEProblem(ODEFunction{false}(de),[1.0,1.0],(0.0,1.0),[1.5,1.0,3.0,1.0])

prob = ODEProblem(de, [1.0, 1.0], (0.0, 1.0), [1.5, 1.0, 3.0, 1.0], check_length = false)
@test prob.kwargs[:callback] isa ModelingToolkit.DiffEqCallbacks.DiscreteCallback

sol = solve(prob, Tsit5(); callback = prob.kwargs[:callback],
            tstops = prob.tspan[1]:0.1:prob.tspan[2], verbose = false)

# Direct implementation
function lotka(du, u, p, t)
    x = u[1]
    y = u[2]
    du[1] = p[1] * x - p[2] * x * y
    du[2] = -p[3] * y + p[4] * x * y
end

prob2 = ODEProblem(lotka, [1.0, 1.0], (0.0, 1.0), [1.5, 1.0, 3.0, 1.0])
function periodic_difference_affect!(int)
    int.u = [int.u[1] + int.u[2], int.u[1] + 1]
    return nothing
end

difference_cb = ModelingToolkit.PeriodicCallback(periodic_difference_affect!, 0.1)

sol2 = solve(prob2, Tsit5(); callback = difference_cb,
             tstops = collect(prob.tspan[1]:0.1:prob.tspan[2])[2:end], verbose = false)

@test sol(0:0.01:1)[x] ≈ sol2(0:0.01:1)[1, :]
@test sol(0:0.01:1)[y] ≈ sol2(0:0.01:1)[2, :]

using ModelingToolkit

function submodel(; name)
    @variables t y(t)
    @parameters A[1:5]
    A = collect(A)
    D = Differential(t)
    ODESystem(D(y) ~ sum(A) * y; name = name)
end

# Buid system
@named sys1 = submodel()
@named sys2 = submodel()

@variables t
@named sys = ODESystem([0 ~ sys1.y + sys2.y], t; systems = [sys1, sys2])

# DelayDiffEq
using ModelingToolkit: hist
@variables t x(t) y(t)
D = Differential(t)
xₜ₋₁ = hist(x, t - 1)
eqs = [D(x) ~ x * y
       D(y) ~ y * x - xₜ₋₁]
@named sys = ODESystem(eqs, t)

# register
using StaticArrays
using SymbolicUtils: term
using SymbolicUtils.Code
using Symbolics: unwrap, wrap
function foo(a::Num, ms::AbstractVector)
    a = unwrap(a)
    ms = map(unwrap, ms)
    wrap(term(foo, a, term(SVector, ms...)))
end
foo(a, ms::AbstractVector) = a + sum(ms)
@variables t x(t) ms[1:3](t)
D = Differential(t)
ms = collect(ms)
eqs = [D(x) ~ foo(x, ms); D.(ms) .~ 1]
@named sys = ODESystem(eqs, t, [x; ms], [])
@named emptysys = ODESystem(Equation[], t)
@named outersys = compose(emptysys, sys)
prob = ODEProblem(outersys, [sys.x => 1.0; collect(sys.ms) .=> 1:3], (0, 1.0))
@test_nowarn solve(prob, Tsit5())

# x/x
@variables t x(t)
@named sys = ODESystem([D(x) ~ x / x], t)
@test equations(alias_elimination(sys)) == [D(x) ~ 1]

# observed variable handling
@variables t x(t) RHS(t)
@parameters τ
D = Differential(t)
@named fol = ODESystem([D(x) ~ (1 - x) / τ]; observed = [RHS ~ (1 - x) / τ])
@test isequal(RHS, @nonamespace fol.RHS)
RHS2 = RHS
@unpack RHS = fol
@test isequal(RHS, RHS2)

#1413 and 1389
@parameters t α β
@variables x(t) y(t) z(t)
D = Differential(t)

eqs = [
    D(x) ~ 0.1x + 0.9y,
    D(y) ~ 0.5x + 0.5y,
    z ~ α * x - β * y,
]

@named sys = ODESystem(eqs, t, [x, y, z], [α, β])
@test_throws Any ODEFunction(sys)

eqs = copy(eqs)
eqs[end] = D(D(z)) ~ α * x - β * y
@named sys = ODESystem(eqs, t, [x, y, z], [α, β])
@test_throws Any ODEFunction(sys)

@testset "Preface tests" begin
    using OrdinaryDiffEq
    using Symbolics
    using DiffEqBase: isinplace
    using ModelingToolkit
    using SymbolicUtils.Code
    using SymbolicUtils: Sym

    c = [0]
    function f(c, du::AbstractVector{Float64}, u::AbstractVector{Float64}, p, t::Float64)
        c .= [c[1] + 1]
        du .= randn(length(u))
        nothing
    end

    dummy_identity(x, _) = x
    @register_symbolic dummy_identity(x, y)

    u0 = ones(5)
    p0 = Float64[]
    syms = [Symbol(:a, i) for i in 1:5]
    syms_p = Symbol[]

    @assert isinplace(f, 5)
    wf = let buffer = similar(u0), u = similar(u0), p = similar(p0), c = c
        t -> (f(c, buffer, u, p, t); buffer)
    end

    num = hash(f) ⊻ length(u0) ⊻ length(p0)
    buffername = Symbol(:fmi_buffer_, num)

    D = Differential(t)
    us = map(s -> (@variables $s(t))[1], syms)
    ps = map(s -> (@variables $s(t))[1], syms_p)
    buffer, = @variables $buffername[1:length(u0)]
    dummy_var = Sym{Any}(:_) # this is safe because _ cannot be a rvalue in Julia

    ss = Iterators.flatten((us, ps))
    vv = Iterators.flatten((u0, p0))
    defs = Dict{Any, Any}(s => v for (s, v) in zip(ss, vv))

    preface = [Assignment(dummy_var, SetArray(true, term(getfield, wf, Meta.quot(:u)), us))
               Assignment(dummy_var, SetArray(true, term(getfield, wf, Meta.quot(:p)), ps))
               Assignment(buffer, term(wf, t))]
    eqs = map(1:length(us)) do i
        D(us[i]) ~ dummy_identity(buffer[i], us[i])
    end

    @named sys = ODESystem(eqs, t, us, ps; defaults = defs, preface = preface)
    prob = ODEProblem(sys, [], (0.0, 1.0))
    sol = solve(prob, Euler(); dt = 0.1)

    @test c[1] == length(sol)
end

let
    @parameters t
    D = Differential(t)
    @variables x[1:2](t) = zeros(2)
    @variables y(t) = 0
    @parameters k = 1
    eqs = [D(x[1]) ~ x[2]
           D(x[2]) ~ -x[1] - 0.5 * x[2] + k
           y ~ 0.9 * x[1] + x[2]]
    @named sys = ODESystem(eqs, t, vcat(x, [y]), [k])
    sys = structural_simplify(sys)

    u0 = [0.5, 0]
    du0 = 0 .* copy(u0)
    prob = DAEProblem(sys, du0, u0, (0, 50))
    @test prob.u0 ≈ u0
    @test prob.du0 ≈ du0
    @test prob.p ≈ [1]
    sol = solve(prob, IDA())
    @test isapprox(sol[x[1]][end], 1, atol = 1e-3)

    prob = DAEProblem(sys, [D(y) => 0, D(x[1]) => 0, D(x[2]) => 0], Pair[x[1] => 0.5],
                      (0, 50))
    @test prob.u0 ≈ [0.5, 0]
    @test prob.du0 ≈ [0, 0]
    @test prob.p ≈ [1]
    sol = solve(prob, IDA())
    @test isapprox(sol[x[1]][end], 1, atol = 1e-3)

    prob = DAEProblem(sys, [D(y) => 0, D(x[1]) => 0, D(x[2]) => 0], Pair[x[1] => 0.5],
                      (0, 50), [k => 2])
    @test prob.u0 ≈ [0.5, 0]
    @test prob.du0 ≈ [0, 0]
    @test prob.p ≈ [2]
    sol = solve(prob, IDA())
    @test isapprox(sol[x[1]][end], 2, atol = 1e-3)

    # no initial conditions for D(x[1]) and D(x[2]) provided
    @test_throws ArgumentError prob=DAEProblem(sys, Pair[], Pair[], (0, 50))

    prob = ODEProblem(sys, Pair[x[1] => 0], (0, 50))
    sol = solve(prob, Rosenbrock23())
    @test isapprox(sol[x[1]][end], 1, atol = 1e-3)
end

#issue 1475 (mixed numeric type for parameters)
let
    @parameters k1 k2
    @variables t A(t)
    D = Differential(t)
    eqs = [D(A) ~ -k1 * k2 * A]
    @named sys = ODESystem(eqs, t)
    u0map = [A => 1.0]
    pmap = (k1 => 1.0, k2 => 1)
    tspan = (0.0, 1.0)
    prob = ODEProblem(sys, u0map, tspan, pmap)
    @test prob.p === Tuple([(Dict(pmap))[k] for k in values(parameters(sys))])

    pmap = [k1 => 1, k2 => 1]
    tspan = (0.0, 1.0)
    prob = ODEProblem(sys, u0map, tspan, pmap)
    @test eltype(prob.p) === Float64

    pmap = Pair{Any, Union{Int, Float64}}[k1 => 1, k2 => 1.0]
    tspan = (0.0, 1.0)
    prob = ODEProblem(sys, u0map, tspan, pmap, use_union = true)
    @test eltype(prob.p) === Union{Float64, Int}
end

let
    @variables t s(t) I(t) r(t)
    @parameters N
    @test_throws Any @named tmp = ODESystem([s + I + r ~ N])
end

let
    @parameters C L R
    @variables t q(t) p(t) F(t)
    D = Differential(t)

    eqs = [D(q) ~ -p / L - F
           D(p) ~ q / C
           0 ~ q / C - R * F]

    @named sys = ODESystem(eqs, t)
    @test length(equations(structural_simplify(sys))) == 2
end

let
    eq_to_lhs(eq) = eq.lhs - eq.rhs ~ 0
    eqs_to_lhs(eqs) = eq_to_lhs.(eqs)

    @parameters σ=10 ρ=28 β=8 / 3 sigma rho beta
    @variables t t2 x(t)=1 y(t)=0 z(t)=0 x2(t2)=1 y2(t2)=0 z2(t2)=0 u[1:3](t2)

    D = Differential(t)
    D2 = Differential(t2)

    eqs = [D(x) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z]

    eqs2 = [
        D2(y2) ~ x2 * (rho - z2) - y2,
        D2(x2) ~ sigma * (y2 - x2),
        D2(z2) ~ x2 * y2 - beta * z2,
    ]

    # array u
    eqs3 = [D2(u[1]) ~ sigma * (u[2] - u[1]),
        D2(u[2]) ~ u[1] * (rho - u[3]) - u[2],
        D2(u[3]) ~ u[1] * u[2] - beta * u[3]]
    eqs3 = eqs_to_lhs(eqs3)

    eqs4 = [
        D2(y2) ~ x2 * (rho - z2) - y2,
        D2(x2) ~ sigma * (y2 - x2),
        D2(z2) ~ y2 - beta * z2, # missing x2 term
    ]

    @named sys1 = ODESystem(eqs)
    @named sys2 = ODESystem(eqs2)
    @named sys3 = ODESystem(eqs3, t2)
    ssys3 = structural_simplify(sys3)
    @named sys4 = ODESystem(eqs4)

    @test ModelingToolkit.isisomorphic(sys1, sys2)
    @test !ModelingToolkit.isisomorphic(sys1, sys3)
    @test ModelingToolkit.isisomorphic(sys1, ssys3) # I don't call structural_simplify in isisomorphic
    @test !ModelingToolkit.isisomorphic(sys1, sys4)

    # 1281
    iv2 = only(independent_variables(sys2))
    @test isequal(only(independent_variables(convert_system(ODESystem, sys1, iv2))), iv2)
end

let
    @variables t
    vars = @variables sP(t) spP(t) spm(t) sph(t)
    pars = @parameters a b
    eqs = [sP ~ 1
           spP ~ sP
           spm ~ a
           sph ~ b
           spm ~ 0
           sph ~ a]
    @named sys = ODESystem(eqs, t, vars, pars)
    @test_throws ModelingToolkit.ExtraEquationsSystemException structural_simplify(sys)
end

# 1561
let
    vars = @variables x y
    arr = ModelingToolkit.varmap_to_vars([x => 0.0, y => [0.0, 1.0]], vars) #error
    sol = Union{Float64, Vector{Float64}}[0.0, [0.0, 1.0]]
    @test arr == sol
    @test typeof(arr) == typeof(sol)
end

let
    @parameters t

    u = collect(first(@variables u[1:4](t)))
    Dt = Differential(t)

    eqs = [Differential(t)(u[2]) - 1.1u[1] ~ 0
           Differential(t)(u[3]) - 1.1u[2] ~ 0
           u[1] ~ 0.0
           u[4] ~ 0.0]

    ps = []

    @named sys = ODESystem(eqs, t, u, ps)
    @test_nowarn simpsys = structural_simplify(sys)
end

# https://github.com/SciML/ModelingToolkit.jl/issues/1583
let
    @parameters k
    @variables t A(t)
    D = Differential(t)
    eqs = [D(A) ~ -k * A]
    @named osys = ODESystem(eqs, t)
    oprob = ODEProblem(osys, [A => 1.0], (0.0, 10.0), [k => 1.0]; check_length = false)
    @test_nowarn sol = solve(oprob, Tsit5())
end

let
    function sys1(; name)
        vars = @variables x(t)=0.0 dx(t)=0.0

        ODESystem([D(x) ~ dx], t, vars, []; name, defaults = [D(x) => x])
    end

    function sys2(; name)
        @named s1 = sys1()

        ODESystem(Equation[], t, [], []; systems = [s1], name)
    end

    s1′ = sys1(; name = :s1)
    @named s2 = sys2()
    @unpack s1 = s2
    @test isequal(s1, s1′)

    defs = Dict(s1.dx => 0.0, D(s1.x) => s1.x, s1.x => 0.0)
    @test isequal(ModelingToolkit.defaults(s2), defs)
end
