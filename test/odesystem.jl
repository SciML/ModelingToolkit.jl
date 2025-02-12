using ModelingToolkit, StaticArrays, LinearAlgebra
using ModelingToolkit: get_metadata, MTKParameters
using SymbolicIndexingInterface
using OrdinaryDiffEq, Sundials
using DiffEqBase, SparseArrays
using StaticArrays
using Test
using SymbolicUtils: issym
using ForwardDiff
using ModelingToolkit: value
using ModelingToolkit: t_nounits as t, D_nounits as D

# Define some variables
@parameters σ ρ β
@constants κ = 1
@variables x(t) y(t) z(t)
@parameters k

# Define a differential equation
eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z * κ]

ModelingToolkit.toexpr.(eqs)[1]
@named de = ODESystem(eqs, t; defaults = Dict(x => 1))
subed = substitute(de, [σ => k])
ssort(eqs) = sort(eqs, by = string)
@test isequal(ssort(parameters(subed)), [k, β, ρ])
@test isequal(equations(subed),
    [D(x) ~ k * (y - x)
     D(y) ~ (ρ - z) * x - y
     D(z) ~ x * y - β * κ * z])
@named des[1:3] = ODESystem(eqs, t)
@test length(unique(x -> ModelingToolkit.get_tag(x), des)) == 1

@test eval(toexpr(de)) == de
@test hash(deepcopy(de)) == hash(de)

generate_function(de)

function test_diffeq_inference(name, sys, iv, dvs, ps)
    @testset "ODESystem construction: $name" begin
        @test isequal(independent_variables(sys)[1], value(iv))
        @test length(independent_variables(sys)) == 1
        @test isempty(setdiff(Set(unknowns(sys)), Set(value.(dvs))))
        @test isempty(setdiff(Set(parameters(sys)), Set(value.(ps))))
    end
end

test_diffeq_inference("standard", de, t, [x, y, z], [ρ, σ, β])
generate_function(de, [x, y, z], [σ, ρ, β])
jac_expr = generate_jacobian(de)
jac = calculate_jacobian(de)
jacfun = eval(jac_expr[2])

de = complete(de)
for f in [
    ODEFunction(de, [x, y, z], [σ, ρ, β], tgrad = true, jac = true),
    eval(ODEFunctionExpr(de, [x, y, z], [σ, ρ, β], tgrad = true, jac = true))
]
    # system
    @test f.sys === de

    # iip
    du = zeros(3)
    u = collect(1:3)
    p = ModelingToolkit.MTKParameters(de, [σ, ρ, β] .=> 4.0:6.0)
    f.f(du, u, p, 0.1)
    @test du == [4, 0, -16]

    # oop
    du = @SArray zeros(3)
    u = SVector(1:3...)
    p = ModelingToolkit.MTKParameters(de, SVector{3}([σ, ρ, β] .=> 4.0:6.0))
    @test f.f(u, p, 0.1) === @SArray [4.0, 0.0, -16.0]

    # iip vs oop
    du = zeros(3)
    g = similar(du)
    J = zeros(3, 3)
    u = collect(1:3)
    p = ModelingToolkit.MTKParameters(de, [σ, ρ, β] .=> 4.0:6.0)
    f.f(du, u, p, 0.1)
    @test du == f(u, p, 0.1)
    f.tgrad(g, u, p, t)
    @test g == f.tgrad(u, p, t)
    f.jac(J, u, p, t)
    @test J == f.jac(u, p, t)
end

#check iip_config
f = eval(ODEFunctionExpr(de, [x, y, z], [σ, ρ, β], iip_config = (false, true)))
du = zeros(3)
u = collect(1:3)
p = ModelingToolkit.MTKParameters(de, [σ, ρ, β] .=> 4.0:6.0)
f.f(du, u, p, 0.1)
@test du == [4, 0, -16]
@test_throws ArgumentError f.f(u, p, 0.1)

#check iip
f = eval(ODEFunctionExpr(de, [x, y, z], [σ, ρ, β]))
f2 = ODEFunction(de, [x, y, z], [σ, ρ, β])
@test SciMLBase.isinplace(f) === SciMLBase.isinplace(f2)
@test SciMLBase.specialization(f) === SciMLBase.specialization(f2)
for iip in (true, false)
    f = eval(ODEFunctionExpr{iip}(de, [x, y, z], [σ, ρ, β]))
    f2 = ODEFunction{iip}(de, [x, y, z], [σ, ρ, β])
    @test SciMLBase.isinplace(f) === SciMLBase.isinplace(f2) === iip
    @test SciMLBase.specialization(f) === SciMLBase.specialization(f2)

    for specialize in (SciMLBase.AutoSpecialize, SciMLBase.FullSpecialize)
        f = eval(ODEFunctionExpr{iip, specialize}(de, [x, y, z], [σ, ρ, β]))
        f2 = ODEFunction{iip, specialize}(de, [x, y, z], [σ, ρ, β])
        @test SciMLBase.isinplace(f) === SciMLBase.isinplace(f2) === iip
        @test SciMLBase.specialization(f) === SciMLBase.specialization(f2) === specialize
    end
end

#check sparsity
f = eval(ODEFunctionExpr(de, [x, y, z], [σ, ρ, β], sparsity = true))
@test f.sparsity == ModelingToolkit.jacobian_sparsity(de)

f = eval(ODEFunctionExpr(de, [x, y, z], [σ, ρ, β], sparsity = false))
@test isnothing(f.sparsity)

eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y * t,
    D(z) ~ x * y - β * z * κ]
@named de = ODESystem(eqs, t)
de = complete(de)
ModelingToolkit.calculate_tgrad(de)

tgrad_oop, tgrad_iip = eval.(ModelingToolkit.generate_tgrad(de))

u = SVector(1:3...)
p = ModelingToolkit.MTKParameters(de, SVector{3}([σ, ρ, β] .=> 4.0:6.0))
@test tgrad_oop(u, p, t) == [0.0, -u[2], 0.0]
du = zeros(3)
tgrad_iip(du, u, p, t)
@test du == [0.0, -u[2], 0.0]

@parameters σ′(t - 1)
eqs = [D(x) ~ σ′ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z * κ]
@named de = ODESystem(eqs, t)
test_diffeq_inference("global iv-varying", de, t, (x, y, z), (σ′, ρ, β))

f = eval(generate_function(de, [x, y, z], [σ′, ρ, β])[2])
du = [0.0, 0.0, 0.0]
f(du, [1.0, 2.0, 3.0], [x -> x + 7, 2, 3], 5.0)
@test du ≈ [11, -3, -7]

@parameters σ(..)
eqs = [D(x) ~ σ(t - 1) * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z * κ]
@named de = ODESystem(eqs, t)
test_diffeq_inference("single internal iv-varying", de, t, (x, y, z), (σ, ρ, β))
f = eval(generate_function(de, [x, y, z], [σ, ρ, β])[2])
du = [0.0, 0.0, 0.0]
f(du, [1.0, 2.0, 3.0], [x -> x + 7, 2, 3], 5.0)
@test du ≈ [11, -3, -7]

eqs = [D(x) ~ x + 10σ(t - 1) + 100σ(t - 2) + 1000σ(t^2)]
@named de = ODESystem(eqs, t)
test_diffeq_inference("many internal iv-varying", de, t, (x,), (σ,))
f = eval(generate_function(de, [x], [σ])[2])
du = [0.0]
f(du, [1.0], [t -> t + 2], 5.0)
@test du ≈ [27561]

# Conversion to first-order ODEs #17
D3 = D^3
D2 = D^2
@variables u(t) uˍtt(t) uˍt(t) xˍt(t)
eqs = [D3(u) ~ 2(D2(u)) + D(u) + D(x) + 1
       D2(x) ~ D(x) + 2]
@named de = ODESystem(eqs, t)
de1 = ode_order_lowering(de)
lowered_eqs = [D(uˍtt) ~ 2uˍtt + uˍt + xˍt + 1
               D(xˍt) ~ xˍt + 2
               D(uˍt) ~ uˍtt
               D(u) ~ uˍt
               D(x) ~ xˍt]

#@test de1 == ODESystem(lowered_eqs)

# issue #219
@test all(isequal.(
    [ModelingToolkit.var_from_nested_derivative(eq.lhs)[1]
     for eq in equations(de1)],
    unknowns(@named lowered = ODESystem(lowered_eqs, t))))

test_diffeq_inference("first-order transform", de1, t, [uˍtt, xˍt, uˍt, u, x], [])
du = zeros(5)
ODEFunction(complete(de1), [uˍtt, xˍt, uˍt, u, x], [])(du, ones(5), nothing, 0.1)
@test du == [5.0, 3.0, 1.0, 1.0, 1.0]

# Internal calculations
@parameters σ
a = y - x
eqs = [D(x) ~ σ * a,
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z * κ]
@named de = ODESystem(eqs, t)
generate_function(de, [x, y, z], [σ, ρ, β])
jac = calculate_jacobian(de)
@test ModelingToolkit.jacobian_sparsity(de).colptr == sparse(jac).colptr
@test ModelingToolkit.jacobian_sparsity(de).rowval == sparse(jac).rowval

f = ODEFunction(complete(de), [x, y, z], [σ, ρ, β])

@parameters A B C
_x = y / C
eqs = [D(x) ~ -A * x,
    D(y) ~ A * x - B * _x]
@named de = ODESystem(eqs, t)
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
de = complete(modelingtoolkitize(prob))
ODEFunction(de)(similar(prob.u0), prob.u0, (prob.p,), 0.1)

function lotka(du, u, p, t)
    x = u[1]
    y = u[2]
    du[1] = p[1] * x - p[2] * x * y
    du[2] = -p[3] * y + p[4] * x * y
end

prob = ODEProblem(lotka, [1.0, 1.0], (0.0, 1.0), [1.5, 1.0, 3.0, 1.0])

de = complete(modelingtoolkitize(prob))
ODEFunction(de)(similar(prob.u0), prob.u0, (prob.p,), 0.1)

# automatic unknown detection for DAEs
@parameters k₁ k₂ k₃
@variables y₁(t) y₂(t) y₃(t)
# reorder the system just to be a little spicier
eqs = [D(y₁) ~ -k₁ * y₁ + k₃ * y₂ * y₃,
    0 ~ y₁ + y₂ + y₃ - 1,
    D(y₂) ~ k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃ * κ]
@named sys = ODESystem(eqs, t, defaults = [k₁ => 100, k₂ => 3e7, y₁ => 1.0])
sys = complete(sys)
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
@test prob1.f.sys == sys
prob12 = ODEProblem(sys, u0, tspan, [k₁ => 0.04, k₂ => 3e7, k₃ => 1e4])
prob13 = ODEProblem(sys, u0, tspan, (k₁ => 0.04, k₂ => 3e7, k₃ => 1e4))
prob14 = ODEProblem(sys, u0, tspan, p2)
for p in [prob1, prob14]
    @test p.p == MTKParameters(sys, [k₁ => 0.04, k₂ => 3e7, k₃ => 1e4])
    @test Set(Num.(unknowns(sys)) .=> p.u0) == Set([y₁ => 1, y₂ => 0, y₃ => 0])
end
# test remake with symbols
p3 = [k₁ => 0.05,
    k₂ => 2e7,
    k₃ => 1.1e4]
u01 = [y₁ => 1, y₂ => 1, y₃ => 1]
prob_pmap = remake(prob14; p = p3, u0 = u01)
prob_dpmap = remake(prob14; p = Dict(p3), u0 = Dict(u01))
for p in [prob_pmap, prob_dpmap]
    @test p.p == MTKParameters(sys, [k₁ => 0.05, k₂ => 2e7, k₃ => 1.1e4])
    @test Set(Num.(unknowns(sys)) .=> p.u0) == Set([y₁ => 1, y₂ => 1, y₃ => 1])
end
sol_pmap = solve(prob_pmap, Rodas5())
sol_dpmap = solve(prob_dpmap, Rodas5())
@test all(isequal(0.05), sol_pmap.(0:10:100, idxs = k₁))

@test sol_pmap.u ≈ sol_dpmap.u

@testset "symbolic remake with nested system" begin
    function makesys(name)
        @parameters a = 1.0
        @variables x(t) = 0.0
        ODESystem([D(x) ~ -a * x], t; name)
    end

    function makecombinedsys()
        sys1 = makesys(:sys1)
        sys2 = makesys(:sys2)
        @parameters b = 1.0
        complete(ODESystem(Equation[], t, [], [b]; systems = [sys1, sys2], name = :foo))
    end

    sys = makecombinedsys()
    @unpack sys1, b = sys
    prob = ODEProblem(sys, Pair[])
    prob_new = SciMLBase.remake(prob, p = Dict(sys1.a => 3.0, b => 4.0),
        u0 = Dict(sys1.x => 1.0))
    @test prob_new.p == MTKParameters(sys, [b => 4.0, sys1.a => 3.0, sys.sys2.a => 1.0])
    @test prob_new.u0 == [1.0, 0.0]
end

# test kwargs
prob2 = ODEProblem(sys, u0, tspan, p, jac = true)
prob3 = ODEProblem(sys, u0, tspan, p, jac = true, sparse = true) #SparseMatrixCSC need to handle
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

@parameters σ β
@variables x(t) y(t) z(t)
eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x - β * y,
    x + z ~ y]
@named sys = ODESystem(eqs, t)
@test all(isequal.(unknowns(sys), [x, y, z]))
@test all(isequal.(parameters(sys), [σ, β]))
@test equations(sys) == eqs
@test ModelingToolkit.isautonomous(sys)

# issue 701
using ModelingToolkit
@parameters a
@variables x(t)
@named sys = ODESystem([D(x) ~ a], t)
@test issym(equations(sys)[1].rhs)

# issue 708
@parameters a
@variables x(t) y(t) z(t)
@named sys = ODESystem([D(x) ~ y, 0 ~ x + z, 0 ~ x - y], t, [z, y, x], [])
asys = add_accumulations(sys)
@variables accumulation_x(t) accumulation_y(t) accumulation_z(t)
eqs = [0 ~ x + z
       0 ~ x - y
       D(accumulation_x) ~ x
       D(accumulation_y) ~ y
       D(accumulation_z) ~ z
       D(x) ~ y]
@test ssort(equations(asys)) == ssort(eqs)
@variables ac(t)
asys = add_accumulations(sys, [ac => (x + y)^2])
eqs = [0 ~ x + z
       0 ~ x - y
       D(ac) ~ (x + y)^2
       D(x) ~ y]
@test ssort(equations(asys)) == ssort(eqs)

sys2 = ode_order_lowering(sys)
M = ModelingToolkit.calculate_massmatrix(sys2)
@test M == Diagonal([1, 0, 0])

# issue #609
@variables x1(t) x2(t)

eqs = [
    D(x1) ~ -x1,
    0 ~ x1 - x2
]
@named sys = ODESystem(eqs, t)
@test isequal(ModelingToolkit.get_iv(sys), t)
@test isequal(unknowns(sys), [x1, x2])
@test isempty(parameters(sys))

# one equation ODESystem test
@parameters r
@variables x(t)
eq = D(x) ~ r * x
@named ode = ODESystem(eq, t)
@test equations(ode) == [eq]
# issue #808
@testset "Combined system name collisions" begin
    function makesys(name)
        @parameters a
        @variables x(t) f(t)

        ODESystem([D(x) ~ -a * x + f], t; name)
    end

    function issue808()
        sys1 = makesys(:sys1)
        sys2 = makesys(:sys1)

        @test_throws ArgumentError ODESystem([sys2.f ~ sys1.x, D(sys1.f) ~ 0], t,
            systems = [sys1, sys2], name = :foo)
    end
    issue808()
end

#Issue 998
pars = []
vars = @variables((u1,))
eqs = [
    D(u1) ~ 1
]
@test_throws ArgumentError ODESystem(eqs, t, vars, pars, name = :foo)

#Issue 1063/998
pars = [t]
vars = @variables((u1(t),))
@test_throws ArgumentError ODESystem(eqs, t, vars, pars, name = :foo)

@parameters w
der = Differential(w)
eqs = [
    der(u1) ~ t
]
@test_throws ArgumentError ModelingToolkit.ODESystem(eqs, t, vars, pars, name = :foo)

@variables x(t)
@parameters M b k
eqs = [D(D(x)) ~ -b / M * D(x) - k / M * x]
ps = [M, b, k]
default_u0 = [D(x) => 0.0, x => 10.0]
default_p = [M => 1.0, b => 1.0, k => 1.0]
@named sys = ODESystem(eqs, t, [x], ps; defaults = [default_u0; default_p], tspan)
sys = ode_order_lowering(sys)
sys = complete(sys)
prob = ODEProblem(sys)
sol = solve(prob, Tsit5())
@test sol.t[end] == tspan[end]
@test sum(abs, sol.u[end]) < 1
prob = ODEProblem{false}(sys; u0_constructor = x -> SVector(x...))
@test prob.u0 isa SVector

# check_eqs_u0 kwarg test
@variables x1(t) x2(t)
eqs = [D(x1) ~ -x1]
@named sys = ODESystem(eqs, t, [x1, x2], [])
sys = complete(sys)
@test_throws ArgumentError ODEProblem(sys, [1.0, 1.0], (0.0, 1.0))
@test_nowarn ODEProblem(sys, [1.0, 1.0], (0.0, 1.0), check_length = false)

# check inputs
let
    @parameters f k d
    @variables x(t) ẋ(t)
    δ = D

    eqs = [δ(x) ~ ẋ, δ(ẋ) ~ f - k * x - d * ẋ]
    @named sys = ODESystem(eqs, t, [x, ẋ], [f, d, k]; controls = [f])

    calculate_control_jacobian(sys)

    @test isequal(calculate_control_jacobian(sys),
        reshape(Num[0, 1], 2, 1))
end

# issue 1109
let
    @variables x(t)[1:3, 1:3]
    @named sys = ODESystem(D.(x) .~ x, t)
    @test_nowarn structural_simplify(sys)
end

# Array vars
using Symbolics: unwrap, wrap
using LinearAlgebra
sts = @variables x(t)[1:3]=[1, 2, 3.0] y(t)=1.0
ps = @parameters p[1:3] = [1, 2, 3]
eqs = [collect(D.(x) .~ x)
       D(y) ~ norm(collect(x)) * y - x[1]]
@named sys = ODESystem(eqs, t, sts, ps)
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
      2sol[x[1]] + 3sol[x[3]] + sol[norm(x)]
@test sol[x .+ [y, 2y, 3y]] ≈ map((x...) -> [x...],
    map((x, y) -> x[1] .+ y, sol[x], sol[y]),
    map((x, y) -> x[2] .+ 2y, sol[x], sol[y]),
    map((x, y) -> x[3] .+ 3y, sol[x], sol[y]))

using ModelingToolkit

function submodel(; name)
    @variables y(t)
    @parameters A[1:5]
    A = collect(A)
    ODESystem(D(y) ~ sum(A) * y, t; name = name)
end

# Build system
@named sys1 = submodel()
@named sys2 = submodel()

@named sys = ODESystem([0 ~ sys1.y + sys2.y], t; systems = [sys1, sys2])

# DelayDiffEq
using ModelingToolkit: hist
@variables x(t) y(t)
xₜ₋₁ = hist(x, t - 1)
eqs = [D(x) ~ x * y
       D(y) ~ y * x - xₜ₋₁]
@named sys = ODESystem(eqs, t)

# register
using StaticArrays
using SymbolicUtils: term
using SymbolicUtils.Code
using Symbolics: unwrap, wrap, @register_symbolic
foo(a, ms::AbstractVector) = a + sum(ms)
@register_symbolic foo(a, ms::AbstractVector)
@variables x(t) ms(t)[1:3]
eqs = [D(x) ~ foo(x, ms); D(ms) ~ ones(3)]
@named sys = ODESystem(eqs, t, [x; ms], [])
@named emptysys = ODESystem(Equation[], t)
@mtkbuild outersys = compose(emptysys, sys)
prob = ODEProblem(
    outersys, [outersys.sys.x => 1.0; collect(outersys.sys.ms .=> 1:3)], (0, 1.0))
@test_nowarn solve(prob, Tsit5())

# array equations
bar(x, p) = p * x
@register_array_symbolic bar(x::AbstractVector, p::AbstractMatrix) begin
    size = size(x)
    eltype = promote_type(eltype(x), eltype(p))
end
@parameters p[1:3, 1:3]
eqs = [D(x) ~ foo(x, ms); D(ms) ~ bar(ms, p)]
@named sys = ODESystem(eqs, t)
@named emptysys = ODESystem(Equation[], t)
@mtkbuild outersys = compose(emptysys, sys)
prob = ODEProblem(
    outersys, [sys.x => 1.0, sys.ms => 1:3], (0.0, 1.0), [sys.p => ones(3, 3)])
@test_nowarn solve(prob, Tsit5())
obsfn = ModelingToolkit.build_explicit_observed_function(
    outersys, bar(3outersys.sys.ms, 3outersys.sys.p))
@test_nowarn obsfn(sol.u[1], prob.p, sol.t[1])

# x/x
@variables x(t)
@named sys = ODESystem([D(x) ~ x / x], t)
@test equations(alias_elimination(sys)) == [D(x) ~ 1]

# observed variable handling
@variables x(t) RHS(t)
@parameters τ
@named fol = ODESystem([D(x) ~ (1 - x) / τ], t; observed = [RHS ~ (1 - x) / τ])
@test isequal(RHS, @nonamespace fol.RHS)
RHS2 = RHS
@unpack RHS = fol
@test isequal(RHS, RHS2)

#1413 and 1389
@parameters α β
@variables x(t) y(t) z(t)

eqs = [
    D(x) ~ 0.1x + 0.9y,
    D(y) ~ 0.5x + 0.5y,
    z ~ α * x - β * y
]

@named sys = ODESystem(eqs, t, [x, y, z], [α, β])
sys = complete(sys)
@test_throws Any ODEFunction(sys)

eqs = copy(eqs)
eqs[end] = D(D(z)) ~ α * x - β * y
@named sys = ODESystem(eqs, t, [x, y, z], [α, β])
sys = complete(sys)
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
    sys = complete(sys)
    prob = ODEProblem(sys, [], (0.0, 1.0))
    sol = solve(prob, Euler(); dt = 0.1)

    @test c[1] == length(sol)
end

let
    x = map(xx -> xx(t), Symbolics.variables(:x, 1:2, T = SymbolicUtils.FnType))
    @variables y(t) = 0
    @parameters k = 1
    eqs = [D(x[1]) ~ x[2]
           D(x[2]) ~ -x[1] - 0.5 * x[2] + k
           y ~ 0.9 * x[1] + x[2]]
    @named sys = ODESystem(eqs, t, vcat(x, [y]), [k], defaults = Dict(x .=> 0))
    sys = structural_simplify(sys)

    u0 = [0.5, 0]
    du0 = 0 .* copy(u0)
    prob = DAEProblem(sys, du0, u0, (0, 50))
    @test prob.u0 ≈ u0
    @test prob.du0 ≈ du0
    @test vcat(prob.p...) ≈ [1]
    sol = solve(prob, IDA())
    @test sol[y] ≈ 0.9 * sol[x[1]] + sol[x[2]]
    @test isapprox(sol[x[1]][end], 1, atol = 1e-3)

    prob = DAEProblem(sys, [D(y) => 0, D(x[1]) => 0, D(x[2]) => 0], Pair[x[1] => 0.5],
        (0, 50))
    @test prob.u0 ≈ [0.5, 0]
    @test prob.du0 ≈ [0, 0]
    @test vcat(prob.p...) ≈ [1]
    sol = solve(prob, IDA())
    @test isapprox(sol[x[1]][end], 1, atol = 1e-3)

    prob = DAEProblem(sys, [D(y) => 0, D(x[1]) => 0, D(x[2]) => 0], Pair[x[1] => 0.5],
        (0, 50), [k => 2])
    @test prob.u0 ≈ [0.5, 0]
    @test prob.du0 ≈ [0, 0]
    @test vcat(prob.p...) ≈ [2]
    sol = solve(prob, IDA())
    @test isapprox(sol[x[1]][end], 2, atol = 1e-3)

    # no initial conditions for D(x[1]) and D(x[2]) provided
    @test_throws ModelingToolkit.MissingVariablesError prob=DAEProblem(
        sys, Pair[], Pair[], (0, 50))

    prob = ODEProblem(sys, Pair[x[1] => 0], (0, 50))
    sol = solve(prob, Rosenbrock23())
    @test isapprox(sol[x[1]][end], 1, atol = 1e-3)
end

#issue 1475 (mixed numeric type for parameters)
let
    @parameters k1 k2::Int
    @variables A(t)
    eqs = [D(A) ~ -k1 * k2 * A]
    @named sys = ODESystem(eqs, t)
    sys = complete(sys)
    u0map = [A => 1.0]
    pmap = (k1 => 1.0, k2 => 1)
    tspan = (0.0, 1.0)
    prob = ODEProblem(sys, u0map, tspan, pmap; tofloat = false)
    @test (prob.p...,) == ([1], [1.0]) || (prob.p...,) == ([1.0], [1])

    prob = ODEProblem(sys, u0map, tspan, pmap)
    @test vcat(prob.p...) isa Vector{Float64}

    pmap = [k1 => 1, k2 => 1]
    tspan = (0.0, 1.0)
    prob = ODEProblem(sys, u0map, tspan, pmap)
    @test eltype(vcat(prob.p...)) === Float64

    prob = ODEProblem(sys, u0map, tspan, pmap)
    @test vcat(prob.p...) isa Vector{Float64}

    # No longer supported, Tuple used instead
    # pmap = Pair{Any, Union{Int, Float64}}[k1 => 1, k2 => 1.0]
    # tspan = (0.0, 1.0)
    # prob = ODEProblem(sys, u0map, tspan, pmap, use_union = true)
    # @test eltype(prob.p) === Union{Float64, Int}
end

let
    @parameters C L R
    @variables q(t) p(t) F(t)

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
    @variables x(t)=1 y(t)=0 z(t)=0 x2(t)=1 y2(t)=0 z2(t)=0 u(t)[1:3]

    eqs = [D(x) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z]

    eqs2 = [
        D(y2) ~ x2 * (rho - z2) - y2,
        D(x2) ~ sigma * (y2 - x2),
        D(z2) ~ x2 * y2 - beta * z2
    ]

    # array u
    eqs3 = [D(u[1]) ~ sigma * (u[2] - u[1]),
        D(u[2]) ~ u[1] * (rho - u[3]) - u[2],
        D(u[3]) ~ u[1] * u[2] - beta * u[3]]
    eqs3 = eqs_to_lhs(eqs3)

    eqs4 = [
        D(y2) ~ x2 * (rho - z2) - y2,
        D(x2) ~ sigma * (y2 - x2),
        D(z2) ~ y2 - beta * z2 # missing x2 term
    ]

    @named sys1 = ODESystem(eqs, t)
    @named sys2 = ODESystem(eqs2, t)
    @named sys3 = ODESystem(eqs3, t)
    ssys3 = structural_simplify(sys3)
    @named sys4 = ODESystem(eqs4, t)

    @test ModelingToolkit.isisomorphic(sys1, sys2)
    @test !ModelingToolkit.isisomorphic(sys1, sys3)
    @test ModelingToolkit.isisomorphic(sys1, ssys3) # I don't call structural_simplify in isisomorphic
    @test !ModelingToolkit.isisomorphic(sys1, sys4)

    # 1281
    iv2 = only(independent_variables(sys2))
    @test isequal(only(independent_variables(convert_system(ODESystem, sys1, iv2))), iv2)
end

let
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
    u = collect(first(@variables u(t)[1:4]))
    Dt = D

    eqs = [Differential(t)(u[2]) - 1.1u[1] ~ 0
           Differential(t)(u[3]) - 1.1u[2] ~ 0
           u[1] ~ 0.0
           u[4] ~ 0.0]

    ps = []

    @named sys = ODESystem(eqs, t, u, ps)
    @test_nowarn simpsys = structural_simplify(sys)

    sys = structural_simplify(sys)

    u0 = ModelingToolkit.missing_variable_defaults(sys)
    u0_expected = Pair[s => 0.0 for s in unknowns(sys)]
    @test string(u0) == string(u0_expected)

    u0 = ModelingToolkit.missing_variable_defaults(sys, [1, 2])
    u0_expected = Pair[s => i for (i, s) in enumerate(unknowns(sys))]
    @test string(u0) == string(u0_expected)

    @test_nowarn ODEProblem(sys, u0, (0, 1))
end

# https://github.com/SciML/ModelingToolkit.jl/issues/1583
let
    @parameters k
    @variables A(t)
    eqs = [D(A) ~ -k * A]
    @named osys = ODESystem(eqs, t)
    osys = complete(osys)
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

# https://github.com/SciML/ModelingToolkit.jl/issues/1705
let
    x0 = 0.0
    v0 = 1.0

    kx = -1.0
    kv = -1.0

    tf = 10.0

    ## controller

    function pd_ctrl(; name)
        @parameters kx kv
        @variables u(t) x(t) v(t)

        eqs = [u ~ kx * x + kv * v]
        ODESystem(eqs, t; name)
    end

    @named ctrl = pd_ctrl()

    ## double integrator

    function double_int(; name)
        @variables u(t) x(t) v(t)

        eqs = [D(x) ~ v, D(v) ~ u]
        ODESystem(eqs, t; name)
    end

    @named sys = double_int()

    ## connections

    connections = [sys.u ~ ctrl.u, ctrl.x ~ sys.x, ctrl.v ~ sys.v]

    @named connected = ODESystem(connections, t)
    @named sys_con = compose(connected, sys, ctrl)

    sys_simp = structural_simplify(sys_con)
    true_eqs = [D(sys.x) ~ sys.v
                D(sys.v) ~ ctrl.kv * sys.v + ctrl.kx * sys.x]
    @test isequal(full_equations(sys_simp), true_eqs)
end

let
    @variables x(t) = 1
    @variables y(t) = 1
    @parameters pp = -1
    @named sys4 = ODESystem([D(x) ~ -y; D(y) ~ 1 + pp * y + x], t)
    sys4s = structural_simplify(sys4)
    prob = ODEProblem(sys4s, [x => 1.0, D(x) => 1.0], (0, 1.0))
    @test string.(unknowns(prob.f.sys)) == ["x(t)", "y(t)"]
    @test string.(parameters(prob.f.sys)) == ["pp"]
    @test string.(independent_variables(prob.f.sys)) == ["t"]
end

let
    @parameters P(t) Q(t)
    ∂t = D
    eqs = [∂t(Q) ~ 0.2P
           ∂t(P) ~ -80.0sin(Q)]
    @test_throws ArgumentError @named sys = ODESystem(eqs, t)
end

@parameters C L R
@variables q(t) p(t) F(t)

eqs = [D(q) ~ -p / L - F
       D(p) ~ q / C
       0 ~ q / C - R * F]
testdict = Dict([:name => "test"])
@named sys = ODESystem(eqs, t, metadata = testdict)
@test get_metadata(sys) == testdict

@variables P(t)=NaN Q(t)=NaN
eqs = [D(Q) ~ 1 / sin(P), D(P) ~ log(-cos(Q))]
@named sys = ODESystem(eqs, t, [P, Q], [])
sys = complete(debug_system(sys))
prob = ODEProblem(sys, [], (0.0, 1.0))
@test_throws "log(-cos(Q(t))) errors" prob.f([1, 0], prob.p, 0.0)
@test_throws "/(1, sin(P(t))) output non-finite value" prob.f([0, 2], prob.p, 0.0)

let
    @variables x(t) = 1
    @variables y(t) = 1
    @parameters pp = -1
    der = Differential(t)
    @named sys4 = ODESystem([der(x) ~ -y; der(y) ~ 1 + pp * y + x], t)
    sys4s = structural_simplify(sys4)
    prob = ODEProblem(sys4s, [x => 1.0, D(x) => 1.0], (0, 1.0))
    @test !isnothing(prob.f.sys)
end

# SYS 1:
vars_sub1 = @variables s1(t)
@named sub = ODESystem(Equation[], t, vars_sub1, [])

vars1 = @variables x1(t)
@named sys1 = ODESystem(Equation[], t, vars1, [], systems = [sub])
@named sys2 = ODESystem(Equation[], t, vars1, [], systems = [sys1, sub])

# SYS 2: Extension to SYS 1
vars_sub2 = @variables s2(t)
@named partial_sub = ODESystem(Equation[], t, vars_sub2, [])
@named sub = extend(partial_sub, sub)

# no warnings for systems without events
new_sys2 = @test_nowarn complete(substitute(sys2, Dict(:sub => sub)))
Set(unknowns(new_sys2)) == Set([new_sys2.x1, new_sys2.sys1.x1,
    new_sys2.sys1.sub.s1, new_sys2.sys1.sub.s2,
    new_sys2.sub.s1, new_sys2.sub.s2])

let # Issue https://github.com/SciML/ModelingToolkit.jl/issues/2322
    @parameters a=10 b=a / 10 c=a / 20

    Dt = D

    @variables x(t)=1 z(t)

    eqs = [Dt(x) ~ -b * (x - z),
        0 ~ z - c * x]

    sys = ODESystem(eqs, t; name = :kjshdf)

    sys_simp = structural_simplify(sys)

    @test a ∈ keys(ModelingToolkit.defaults(sys_simp))

    tspan = (0.0, 1)
    prob = ODEProblem(sys_simp, [], tspan)
    sol = solve(prob, Rodas4())
    @test sol(1)[]≈0.6065307685451087 rtol=1e-4
end

# Issue#2599
@variables x(t) y(t)
eqs = [D(x) ~ x * t, y ~ 2x]
@mtkbuild sys = ODESystem(eqs, t; continuous_events = [[y ~ 3] => [x ~ 2]])
prob = ODEProblem(sys, [x => 1.0], (0.0, 10.0))
@test_nowarn solve(prob, Tsit5())

# Issue#2383
@variables x(t)[1:3]
@parameters p[1:3, 1:3]
eqs = [
    D(x) ~ p * x
]
@mtkbuild sys = ODESystem(eqs, t; continuous_events = [[norm(x) ~ 3.0] => [x ~ ones(3)]])
# array affect equations used to not work
prob1 = @test_nowarn ODEProblem(sys, [x => ones(3)], (0.0, 10.0), [p => ones(3, 3)])
sol1 = @test_nowarn solve(prob1, Tsit5())

# array condition equations also used to not work
@mtkbuild sys = ODESystem(
    eqs, t; continuous_events = [[x ~ sqrt(3) * ones(3)] => [x ~ ones(3)]])
# array affect equations used to not work
prob2 = @test_nowarn ODEProblem(sys, [x => ones(3)], (0.0, 10.0), [p => ones(3, 3)])
sol2 = @test_nowarn solve(prob2, Tsit5())

@test sol1 ≈ sol2

# Requires fix in symbolics for `linear_expansion(p * x, D(y))`
@test_skip begin
    @variables x(t)[1:3] y(t)
    @parameters p[1:3, 1:3]
    @test_nowarn @mtkbuild sys = ODESystem([D(x) ~ p * x, D(y) ~ x' * p * x], t)
    @test_nowarn ODEProblem(sys, [x => ones(3), y => 2], (0.0, 10.0), [p => ones(3, 3)])
end

@parameters g L
@variables q₁(t) q₂(t) λ(t) θ(t)

eqs = [D(D(q₁)) ~ -λ * q₁,
    D(D(q₂)) ~ -λ * q₂ - g,
    q₁ ~ L * sin(θ),
    q₂ ~ L * cos(θ)]

@named pend = ODESystem(eqs, t)
@test_nowarn generate_initializesystem(
    pend, u0map = [q₁ => 1.0, q₂ => 0.0], guesses = [λ => 1])

# https://github.com/SciML/ModelingToolkit.jl/issues/2618
@parameters σ ρ β
@variables x(t) y(t) z(t)

eqs = [D(D(x)) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@mtkbuild sys = ODESystem(eqs, t)

u0 = [D(x) => 2.0,
    x => 1.0,
    y => 0.0,
    z => 0.0]

p = [σ => 28.0,
    ρ => 10.0,
    β => 8 / 3]

prob = SteadyStateProblem(sys, u0, p)
@test prob isa SteadyStateProblem
prob = SteadyStateProblem(ODEProblem(sys, u0, (0.0, 10.0), p))
@test prob isa SteadyStateProblem

# Issue#2344
using ModelingToolkitStandardLibrary.Blocks

function FML2(; name)
    @parameters begin
        k2[1:1] = [1.0]
    end
    systems = @named begin
        constant = Constant(k = k2[1])
    end
    @variables begin
        x(t) = 0
    end
    eqs = [
        D(x) ~ constant.output.u + k2[1]
    ]
    ODESystem(eqs, t; systems, name)
end

@mtkbuild model = FML2()

@test isequal(ModelingToolkit.defaults(model)[model.constant.k], model.k2[1])
@test_nowarn ODEProblem(model, [], (0.0, 10.0))

# Issue#2477
function RealExpression(; name, y)
    vars = @variables begin
        u(t)
    end
    eqns = [
        u ~ y
    ]
    sys = ODESystem(eqns, t, vars, []; name)
end

function RealExpressionSystem(; name)
    vars = @variables begin
        x(t)
        z(t)[1:1]
    end # doing a collect on z doesn't work either. 
    @named e1 = RealExpression(y = x) # This works perfectly. 
    @named e2 = RealExpression(y = z[1]) # This bugs. However, `full_equations(e2)` works as expected. 
    systems = [e1, e2]
    ODESystem(Equation[], t, Iterators.flatten(vars), []; systems, name)
end

@named sys = RealExpressionSystem()
sys = complete(sys)
@test Set(equations(sys)) == Set([sys.e1.u ~ sys.x, sys.e2.u ~ sys.z[1]])
tearing_state = TearingState(expand_connections(sys))
ts_vars = tearing_state.fullvars
orig_vars = unknowns(sys)
@test isempty(setdiff(ts_vars, orig_vars))

# Guesses in hierarchical systems
@variables x(t) y(t)
@named sys = ODESystem(Equation[], t, [x], []; guesses = [x => 1.0])
@named outer = ODESystem(
    [D(y) ~ sys.x + t, 0 ~ t + y - sys.x * y], t, [y], []; systems = [sys])
@test ModelingToolkit.guesses(outer)[sys.x] == 1.0
outer = structural_simplify(outer)
@test ModelingToolkit.get_guesses(outer)[sys.x] == 1.0
prob = ODEProblem(outer, [outer.y => 2.0], (0.0, 10.0))
int = init(prob, Rodas4())
@test int[outer.sys.x] == 1.0

# Ensure indexes of array symbolics are cached appropriately
@variables x(t)[1:2]
@named sys = ODESystem(Equation[], t, [x], [])
sys1 = complete(sys)
@named sys = ODESystem(Equation[], t, [x...], [])
sys2 = complete(sys)
for sys in [sys1, sys2]
    for (sym, idx) in [(x, 1:2), (x[1], 1), (x[2], 2)]
        @test is_variable(sys, sym)
        @test variable_index(sys, sym) == idx
    end
end

@variables x(t)[1:2, 1:2]
@named sys = ODESystem(Equation[], t, [x], [])
sys1 = complete(sys)
@named sys = ODESystem(Equation[], t, [x...], [])
sys2 = complete(sys)
for sys in [sys1, sys2]
    @test is_variable(sys, x)
    @test variable_index(sys, x) == [1 3; 2 4]
    for i in eachindex(x)
        @test is_variable(sys, x[i])
        @test variable_index(sys, x[i]) == variable_index(sys, x)[i]
    end
end

@testset "Non-1-indexed variable array (issue #2670)" begin
    @variables x(t)[0:1] # 0-indexed variable array
    @named sys = ODESystem([x[0] ~ 0.0, D(x[1]) ~ x[0]], t, [x], [])
    @test_nowarn sys = structural_simplify(sys)
    @test equations(sys) == [D(x[1]) ~ 0.0]
end

# Namespacing of array variables
@variables x(t)[1:2]
@named sys = ODESystem(Equation[], t)
@test getname(unknowns(sys, x)) == :sys₊x
@test size(unknowns(sys, x)) == size(x)

# Issue#2667 and Issue#2953
@testset "ForwardDiff through ODEProblem constructor" begin
    @parameters P
    @variables x(t)
    sys = structural_simplify(ODESystem([D(x) ~ P], t, [x], [P]; name = :sys))

    function x_at_1(P)
        prob = ODEProblem(sys, [x => P], (0.0, 1.0), [sys.P => P], use_union = false)
        return solve(prob, Tsit5())(1.0)
    end

    @test_nowarn ForwardDiff.derivative(P -> x_at_1(P), 1.0)
end

@testset "Inplace observed functions" begin
    @parameters P
    @variables x(t)
    sys = structural_simplify(ODESystem([D(x) ~ P], t, [x], [P]; name = :sys))
    obsfn = ModelingToolkit.build_explicit_observed_function(
        sys, [x + 1, x + P, x + t], return_inplace = true)[2]
    ps = ModelingToolkit.MTKParameters(sys, [P => 2.0])
    buffer = zeros(3)
    @test_nowarn obsfn(buffer, [1.0], ps..., 3.0)
    @test buffer ≈ [2.0, 3.0, 4.0]
end

# https://github.com/SciML/ModelingToolkit.jl/issues/2818
@testset "Custom independent variable" begin
    @independent_variables x
    @variables y(x)
    @test_nowarn @named sys = ODESystem([y ~ 0], x)

    # the same, but with @mtkmodel
    @independent_variables x
    @mtkmodel MyModel begin
        @variables begin
            y(x)
        end
        @equations begin
            y ~ 0
        end
    end
    @test_nowarn @mtkbuild sys = MyModel()

    @variables x y(x)
    @test_logs (:warn,) @named sys = ODESystem([y ~ 0], x)

    @parameters T
    D = Differential(T)
    @variables x(T)
    eqs = [D(x) ~ 0.0]
    initialization_eqs = [x ~ T]
    guesses = [x => 0.0]
    @named sys2 = ODESystem(eqs, T; initialization_eqs, guesses)
    prob2 = ODEProblem(structural_simplify(sys2), [], (1.0, 2.0), [])
    sol2 = solve(prob2)
    @test all(sol2[x] .== 1.0)
end

# https://github.com/SciML/ModelingToolkit.jl/issues/2502
@testset "Extend systems with a field that can be nothing" begin
    A = Dict(:a => 1)
    B = Dict(:b => 2)
    @named A1 = ODESystem(Equation[], t, [], [])
    @named B1 = ODESystem(Equation[], t, [], [])
    @named A2 = ODESystem(Equation[], t, [], []; metadata = A)
    @named B2 = ODESystem(Equation[], t, [], []; metadata = B)
    @test ModelingToolkit.get_metadata(extend(A1, B1)) == nothing
    @test ModelingToolkit.get_metadata(extend(A1, B2)) == B
    @test ModelingToolkit.get_metadata(extend(A2, B1)) == A
    @test Set(ModelingToolkit.get_metadata(extend(A2, B2))) == Set(A ∪ B)
end

# https://github.com/SciML/ModelingToolkit.jl/issues/2859
@testset "Initialization with defaults from observed equations (edge case)" begin
    @variables x(t) y(t) z(t)
    eqs = [D(x) ~ 0, y ~ x, D(z) ~ 0]
    defaults = [x => 1, z => y]
    @named sys = ODESystem(eqs, t; defaults)
    ssys = structural_simplify(sys)
    prob = ODEProblem(ssys, [], (0.0, 1.0), [])
    @test prob[x] == prob[y] == prob[z] == 1.0

    @parameters y0
    @variables x(t) y(t) z(t)
    eqs = [D(x) ~ 0, y ~ y0 / x, D(z) ~ y]
    defaults = [y0 => 1, x => 1, z => y]
    @named sys = ODESystem(eqs, t; defaults)
    ssys = structural_simplify(sys)
    prob = ODEProblem(ssys, [], (0.0, 1.0), [])
    @test prob[x] == prob[y] == prob[z] == 1.0
end

@testset "Scalarized parameters in array functions" begin
    @variables u(t)[1:2] x(t)[1:2] o(t)[1:2]
    @parameters p[1:2, 1:2] [tunable = false]
    @named sys = ODESystem(
        [D(u) ~ (sum(u) + sum(x) + sum(p) + sum(o)) * x, o ~ prod(u) * x],
        t, [u..., x..., o...], [p...])
    sys1, = structural_simplify(sys, ([x...], []))
    fn1, = ModelingToolkit.generate_function(sys1; expression = Val{false})
    @test_nowarn fn1(ones(4), (2ones(2), 3ones(2, 2)), 4.0)
    sys2, = structural_simplify(sys, ([x...], []); split = false)
    fn2, = ModelingToolkit.generate_function(sys2; expression = Val{false})
    @test_nowarn fn2(ones(4), 2ones(6), 4.0)
end

# https://github.com/SciML/ModelingToolkit.jl/issues/2969
@testset "Constant substitution" begin
    make_model = function (c_a, c_b; name = nothing)
        @mtkmodel ModelA begin
            @constants begin
                a = c_a
            end
            @variables begin
                x(t)
            end
            @equations begin
                D(x) ~ -a * x
            end
        end

        @mtkmodel ModelB begin
            @constants begin
                b = c_b
            end
            @variables begin
                y(t)
            end
            @components begin
                modela = ModelA()
            end
            @equations begin
                D(y) ~ -b * y
            end
        end
        return ModelB(; name = name)
    end
    c_a, c_b = 1.234, 5.578
    @named sys = make_model(c_a, c_b)
    sys = complete(sys)

    u0 = [sys.y => -1.0, sys.modela.x => -1.0]
    p = defaults(sys)
    prob = ODEProblem(sys, u0, (0.0, 1.0), p)

    # evaluate
    u0_v, p_v, _ = ModelingToolkit.get_u0_p(sys, u0, p)
    @test prob.f(u0_v, p_v, 0.0) == [c_b, c_a]
end

@testset "Independent variable as system property" begin
    @variables x(t)
    @named sys = ODESystem([x ~ t], t)
    @named sys = compose(sys, sys) # nest into a hierarchical system
    @test t === sys.t === sys.sys.t
end

@testset "Substituting preserves parameter dependencies, defaults, guesses" begin
    @parameters p1 p2
    @variables x(t) y(t)
    @named sys = ODESystem([D(x) ~ y + p2], t; parameter_dependencies = [p2 ~ 2p1],
        defaults = [p1 => 1.0, p2 => 2.0], guesses = [p1 => 2.0, p2 => 3.0])
    @parameters p3
    sys2 = substitute(sys, [p1 => p3])
    @test length(parameters(sys2)) == 1
    @test is_parameter(sys2, p3)
    @test !is_parameter(sys2, p1)
    @test length(ModelingToolkit.defaults(sys2)) == 2
    @test ModelingToolkit.defaults(sys2)[p3] == 1.0
    @test length(ModelingToolkit.guesses(sys2)) == 2
    @test ModelingToolkit.guesses(sys2)[p3] == 2.0
end

@testset "Substituting with nested systems" begin
    @parameters p1 p2
    @variables x(t) y(t)
    @named innersys = ODESystem([D(x) ~ y + p2], t; parameter_dependencies = [p2 ~ 2p1],
        defaults = [p1 => 1.0, p2 => 2.0], guesses = [p1 => 2.0, p2 => 3.0])
    @parameters p3 p4
    @named outersys = ODESystem(
        [D(innersys.y) ~ innersys.y + p4], t; parameter_dependencies = [p4 ~ 3p3],
        defaults = [p3 => 3.0, p4 => 9.0], guesses = [p4 => 10.0], systems = [innersys])
    @test_nowarn structural_simplify(outersys)
    @parameters p5
    sys2 = substitute(outersys, [p4 => p5])
    @test_nowarn structural_simplify(sys2)
    @test length(equations(sys2)) == 2
    @test length(parameters(sys2)) == 2
    @test length(full_parameters(sys2)) == 4
    @test all(!isequal(p4), full_parameters(sys2))
    @test any(isequal(p5), full_parameters(sys2))
    @test length(ModelingToolkit.defaults(sys2)) == 4
    @test ModelingToolkit.defaults(sys2)[p5] == 9.0
    @test length(ModelingToolkit.guesses(sys2)) == 3
    @test ModelingToolkit.guesses(sys2)[p5] == 10.0
end

@testset "Observed with inputs" begin
    @variables u(t)[1:2] x(t)[1:2] o(t)[1:2]
    @parameters p[1:4]

    eqs = [D(u[1]) ~ p[1] * u[1] - p[2] * u[1] * u[2] + x[1] + 0.1
           D(u[2]) ~ p[4] * u[1] * u[2] - p[3] * u[2] - x[2]
           o[1] ~ sum(p) * sum(u)
           o[2] ~ sum(p) * sum(x)]

    @named sys = ODESystem(eqs, t, [u..., x..., o], [p...])
    sys1, = structural_simplify(sys, ([x...], [o...]), split = false)

    @test_nowarn ModelingToolkit.build_explicit_observed_function(sys1, u; inputs = [x...])

    obsfn = ModelingToolkit.build_explicit_observed_function(
        sys1, u + x + p[1:2]; inputs = [x...])

    @test obsfn(ones(2), 2ones(2), 3ones(4), 4.0) == 6ones(2)
end

@testset "Passing `nothing` to `u0`" begin
    @variables x(t) = 1
    @mtkbuild sys = ODESystem(D(x) ~ t, t)
    prob = @test_nowarn ODEProblem(sys, nothing, (0.0, 1.0))
    @test_nowarn solve(prob)
end

@testset "ODEs are not DDEs" begin
    @variables x(t)
    @named sys = ODESystem(D(x) ~ x, t)
    @test !ModelingToolkit.is_dde(sys)
    @test is_markovian(sys)
    @named sys2 = ODESystem(Equation[], t; systems = [sys])
    @test !ModelingToolkit.is_dde(sys)
    @test is_markovian(sys)
end

@testset "Issue #2597" begin
    @variables x(t)[1:2]=ones(2) y(t)=1.0

    for eqs in [D(x) ~ x, collect(D(x) .~ x)]
        for dvs in [[x], collect(x)]
            @named sys = ODESystem(eqs, t, dvs, [])
            sys = complete(sys)
            if eqs isa Vector && length(eqs) == 2 && length(dvs) == 2
                @test_nowarn ODEProblem(sys, [], (0.0, 1.0))
            else
                @test_throws [
                    r"array (equations|unknowns)", "structural_simplify", "scalarize"] ODEProblem(
                    sys, [], (0.0, 1.0))
            end
        end
    end
    for eqs in [[D(x) ~ x, D(y) ~ y], [collect(D(x) .~ x); D(y) ~ y]]
        for dvs in [[x, y], [x..., y]]
            @named sys = ODESystem(eqs, t, dvs, [])
            sys = complete(sys)
            if eqs isa Vector && length(eqs) == 3 && length(dvs) == 3
                @test_nowarn ODEProblem(sys, [], (0.0, 1.0))
            else
                @test_throws [
                    r"array (equations|unknowns)", "structural_simplify", "scalarize"] ODEProblem(
                    sys, [], (0.0, 1.0))
            end
        end
    end
end

@testset "Parameter dependencies with constant RHS" begin
    @parameters p
    @test_nowarn ODESystem(Equation[], t; parameter_dependencies = [p ~ 1.0], name = :a)
end

@testset "Variable discovery in arrays of `Num` inside callable symbolic" begin
    @variables x(t) y(t)
    @parameters foo(::AbstractVector)
    sys = @test_nowarn ODESystem(D(x) ~ foo([x, 2y]), t; name = :sys)
    @test length(unknowns(sys)) == 2
    @test any(isequal(y), unknowns(sys))
end

@testset "Inplace observed" begin
    @variables x(t)
    @parameters p[1:2] q
    @mtkbuild sys = ODESystem(D(x) ~ sum(p) * x + q * t, t)
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0), [p => ones(2), q => 2])
    obsfn = ModelingToolkit.build_explicit_observed_function(
        sys, [p..., q], return_inplace = true)[2]
    buf = zeros(3)
    obsfn(buf, prob.u0, prob.p, 0.0)
    @test buf ≈ [1.0, 1.0, 2.0]
end

@testset "`complete` expands connections" begin
    using ModelingToolkitStandardLibrary.Electrical
    @mtkmodel RC begin
        @parameters begin
            R = 1.0
            C = 1.0
            V = 1.0
        end
        @components begin
            resistor = Resistor(R = R)
            capacitor = Capacitor(C = C, v = 0.0)
            source = Voltage()
            constant = Constant(k = V)
            ground = Ground()
        end
        @equations begin
            connect(constant.output, source.V)
            connect(source.p, resistor.p)
            connect(resistor.n, capacitor.p)
            connect(capacitor.n, source.n, ground.g)
        end
    end
    @named sys = RC()
    total_eqs = length(equations(expand_connections(sys)))
    sys2 = complete(sys)
    @test length(equations(sys2)) == total_eqs
end

@testset "`complete` with `split = false` removes the index cache" begin
    @variables x(t)
    @parameters p
    @mtkbuild sys = ODESystem(D(x) ~ p * t, t)
    @test ModelingToolkit.get_index_cache(sys) !== nothing
    sys2 = complete(sys; split = false)
    @test ModelingToolkit.get_index_cache(sys2) === nothing
end

# https://github.com/SciML/SciMLBase.jl/issues/786
@testset "Observed variables dependent on discrete parameters" begin
    @variables x(t) obs(t)
    @parameters c(t)
    @mtkbuild sys = ODESystem(
        [D(x) ~ c * cos(x), obs ~ c], t, [x], [c]; discrete_events = [1.0 => [c ~ c + 1]])
    prob = ODEProblem(sys, [x => 0.0], (0.0, 2pi), [c => 1.0])
    sol = solve(prob, Tsit5())
    @test sol[obs] ≈ 1:7
end

@testset "DAEProblem with array parameters" begin
    @variables x(t)=1.0 y(t) [guess = 1.0]
    @parameters p[1:2] = [1.0, 2.0]
    @mtkbuild sys = ODESystem([D(x) ~ x, y^2 ~ x + sum(p)], t)
    prob = DAEProblem(sys, [D(x) => x, D(y) => D(x) / 2y], [], (0.0, 1.0))
    sol = solve(prob, DFBDF(), abstol = 1e-8, reltol = 1e-8)
    @test sol[x]≈sol[y^2 - sum(p)] atol=1e-5
end

@testset "Symbolic tstops" begin
    @variables x(t) = 1.0
    @parameters p=0.15 q=0.25 r[1:2]=[0.35, 0.45]
    @mtkbuild sys = ODESystem(
        [D(x) ~ p * x + q * t + sum(r)], t; tstops = [0.5p, [0.1, 0.2], [p + 2q], r])
    prob = ODEProblem(sys, [], (0.0, 5.0))
    sol = solve(prob)
    expected_tstops = unique!(sort!(vcat(0.0:0.075:5.0, 0.1, 0.2, 0.65, 0.35, 0.45)))
    @test all(x -> any(isapprox(x, atol = 1e-6), sol.t), expected_tstops)
    prob2 = remake(prob; tspan = (0.0, 10.0))
    sol2 = solve(prob2)
    expected_tstops = unique!(sort!(vcat(0.0:0.075:10.0, 0.1, 0.2, 0.65, 0.35, 0.45)))
    @test all(x -> any(isapprox(x, atol = 1e-6), sol2.t), expected_tstops)

    @variables y(t) [guess = 1.0]
    @mtkbuild sys = ODESystem([D(x) ~ p * x + q * t + sum(r), y^3 ~ 2x + 1],
        t; tstops = [0.5p, [0.1, 0.2], [p + 2q], r])
    prob = DAEProblem(
        sys, [D(y) => 2D(x) / 3y^2, D(x) => p * x + q * t + sum(r)], [], (0.0, 5.0))
    sol = solve(prob, DImplicitEuler())
    expected_tstops = unique!(sort!(vcat(0.0:0.075:5.0, 0.1, 0.2, 0.65, 0.35, 0.45)))
    @test all(x -> any(isapprox(x, atol = 1e-6), sol.t), expected_tstops)
    prob2 = remake(prob; tspan = (0.0, 10.0))
    sol2 = solve(prob2, DImplicitEuler())
    expected_tstops = unique!(sort!(vcat(0.0:0.075:10.0, 0.1, 0.2, 0.65, 0.35, 0.45)))
    @test all(x -> any(isapprox(x, atol = 1e-6), sol2.t), expected_tstops)
end

@testset "Validate input types" begin
    @parameters p d
    @variables X(t)::Int64
    eq = D(X) ~ p - d * X
    @test_throws ArgumentError @mtkbuild osys = ODESystem([eq], t)
    @variables Y(t)[1:3]::String
    eq = D(Y) ~ [p, p, p]
    @test_throws ArgumentError @mtkbuild osys = ODESystem([eq], t)
end

# Test `isequal`
@testset "`isequal`" begin
    @variables X(t)
    @parameters p d
    eq = D(X) ~ p - d * X

    osys1 = complete(ODESystem([eq], t; name = :osys))
    osys2 = complete(ODESystem([eq], t; name = :osys))
    @test osys1 == osys2 # true

    continuous_events = [[X ~ 1.0] => [X ~ X + 5.0]]
    discrete_events = [5.0 => [d ~ d / 2.0]]

    osys1 = complete(ODESystem([eq], t; name = :osys, continuous_events))
    osys2 = complete(ODESystem([eq], t; name = :osys))
    @test osys1 !== osys2

    osys1 = complete(ODESystem([eq], t; name = :osys, discrete_events))
    osys2 = complete(ODESystem([eq], t; name = :osys))
    @test osys1 !== osys2

    osys1 = complete(ODESystem([eq], t; name = :osys, continuous_events))
    osys2 = complete(ODESystem([eq], t; name = :osys, discrete_events))
    @test osys1 !== osys2
end

@testset "dae_order_lowering basic test" begin
    @parameters a
    @variables x(t) y(t) z(t)
    @named dae_sys = ODESystem([
            D(x) ~ y,
            0 ~ x + z,
            0 ~ x - y + z
        ], t, [z, y, x], [])

    lowered_dae_sys = dae_order_lowering(dae_sys)
    @variables x1(t) y1(t) z1(t)
    expected_eqs = [
        0 ~ x + z,
        0 ~ x - y + z,
        Differential(t)(x) ~ y
    ]
    lowered_eqs = equations(lowered_dae_sys)
    sorted_lowered_eqs = sort(lowered_eqs, by = string)
    sorted_expected_eqs = sort(expected_eqs, by = string)
    @test sorted_lowered_eqs == sorted_expected_eqs

    expected_vars = Set([z, y, x])
    lowered_vars = Set(unknowns(lowered_dae_sys))
    @test lowered_vars == expected_vars
end

@testset "dae_order_lowering test with structural_simplify" begin
    @variables x(t) y(t) z(t)
    @parameters M b k
    eqs = [
        D(D(x)) ~ -b / M * D(x) - k / M * x,
        0 ~ y - D(x),
        0 ~ z - x
    ]
    ps = [M, b, k]
    default_u0 = [
        D(x) => 0.0, x => 10.0, y => 0.0, z => 10.0
    ]
    default_p = [M => 1.0, b => 1.0, k => 1.0]
    @named dae_sys = ODESystem(eqs, t, [x, y, z], ps; defaults = [default_u0; default_p])

    simplified_dae_sys = structural_simplify(dae_sys)

    lowered_dae_sys = dae_order_lowering(simplified_dae_sys)
    lowered_dae_sys = complete(lowered_dae_sys)

    tspan = (0.0, 10.0)
    prob = ODEProblem(lowered_dae_sys, nothing, tspan)
    sol = solve(prob, Tsit5())

    @test sol.t[end] == tspan[end]
    @test sum(abs, sol.u[end]) < 1

    prob = ODEProblem{false}(lowered_dae_sys; u0_constructor = x -> SVector(x...))
    @test prob.u0 isa SVector
end
