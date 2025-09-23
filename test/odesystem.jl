using ModelingToolkit, StaticArrays, LinearAlgebra
using ModelingToolkit: get_metadata, MTKParameters, SymbolicDiscreteCallback,
                       SymbolicContinuousCallback
using SymbolicIndexingInterface
using OrdinaryDiffEq, Sundials
using DiffEqBase, SparseArrays
using StaticArrays
using Test
using SymbolicUtils.Code
using SymbolicUtils: Sym, issym
using ForwardDiff
using ModelingToolkit: value
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics
using Symbolics: unwrap
using DiffEqBase: isinplace

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
@named de = System(eqs, t; defaults = Dict(x => 1))
subed = substitute(de, [σ => k])
ssort(eqs) = sort(eqs, by = string)
@test isequal(ssort(parameters(subed)), [k, β, κ, ρ])
@test isequal(equations(subed),
    [D(x) ~ k * (y - x)
     D(y) ~ (ρ - z) * x - y
     D(z) ~ x * y - β * κ * z])
@named des[1:3] = System(eqs, t)
@test length(unique(x -> ModelingToolkit.get_tag(x), des)) == 1

de2 = eval(toexpr(de))
@test issetequal(equations(de2), eqs)
@test issetequal(unknowns(de2), unknowns(de))
@test issetequal(parameters(de2), parameters(de))

function test_diffeq_inference(name, sys, iv, dvs, ps)
    @testset "System construction: $name" begin
        @test isequal(independent_variables(sys)[1], value(iv))
        @test length(independent_variables(sys)) == 1
        @test isempty(setdiff(Set(unknowns(sys)), Set(value.(dvs))))
        @test isempty(setdiff(Set(parameters(sys)), Set(value.(ps))))
    end
end

test_diffeq_inference("standard", de, t, [x, y, z], [ρ, σ, β, κ])
jac_expr = generate_jacobian(de)
jac = calculate_jacobian(de)
jacfun = eval(jac_expr[2])

de = complete(de)
f = ODEFunction(de, tgrad = true, jac = true)
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

#check iip_config
f = ODEFunction(de; iip_config = (false, true))
du = zeros(3)
u = collect(1:3)
p = ModelingToolkit.MTKParameters(de, [σ, ρ, β] .=> 4.0:6.0)
f.f(du, u, p, 0.1)
@test du == [4, 0, -16]
@test_throws ArgumentError f.f(u, p, 0.1)

#check iip
f = eval(ODEFunction(de; expression = Val{true}))
f2 = ODEFunction(de)
@test SciMLBase.isinplace(f) === SciMLBase.isinplace(f2)
@test SciMLBase.specialization(f) === SciMLBase.specialization(f2)
for iip in (true, false)
    f = eval(ODEFunction{iip}(de; expression = Val{true}))
    f2 = ODEFunction{iip}(de)
    @test SciMLBase.isinplace(f) === SciMLBase.isinplace(f2) === iip
    @test SciMLBase.specialization(f) === SciMLBase.specialization(f2)

    for specialize in (SciMLBase.AutoSpecialize, SciMLBase.FullSpecialize)
        f = eval(ODEFunction{iip, specialize}(de; expression = Val{true}))
        f2 = ODEFunction{iip, specialize}(de)
        @test SciMLBase.isinplace(f) === SciMLBase.isinplace(f2) === iip
        @test SciMLBase.specialization(f) === SciMLBase.specialization(f2) === specialize
    end
end

#check sparsity
f = eval(ODEFunction(de, sparsity = true, expression = Val{true}))
@test f.sparsity == ModelingToolkit.jacobian_sparsity(de)

f = eval(ODEFunction(de, sparsity = false, expression = Val{true}))
@test isnothing(f.sparsity)

eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y * t,
    D(z) ~ x * y - β * z * κ]
@named de = System(eqs, t)
de = complete(de)
ModelingToolkit.calculate_tgrad(de)

tgrad_oop, tgrad_iip = eval.(ModelingToolkit.generate_tgrad(de))

u = SVector(1:3...)
p = ModelingToolkit.MTKParameters(de, SVector{3}([σ, ρ, β] .=> 4.0:6.0))
@test tgrad_oop(u, p, t) == [0.0, -u[2], 0.0]
du = zeros(3)
tgrad_iip(du, u, p, t)
@test du == [0.0, -u[2], 0.0]

@parameters (σ::Function)(..)
eqs = [D(x) ~ σ(t - 1) * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z * κ]
@named de = System(eqs, t, [x, y, z], [σ, ρ, β, κ])
test_diffeq_inference("single internal iv-varying", de, t, (x, y, z), (σ, ρ, β, κ))
f = generate_rhs(de, expression = Val{false}, wrap_gfw = Val{true})
du = [0.0, 0.0, 0.0]
f(du, [1.0, 2.0, 3.0], [x -> x + 7, 2, 3, 1], 5.0)
@test du ≈ [11, -3, -7]

eqs = [D(x) ~ x + 10σ(t - 1) + 100σ(t - 2) + 1000σ(t^2)]
@named de = System(eqs, t)
test_diffeq_inference("many internal iv-varying", de, t, (x,), (σ,))
f = generate_rhs(de, expression = Val{false}, wrap_gfw = Val{true})
du = [0.0]
f(du, [1.0], [t -> t + 2], 5.0)
@test du ≈ [27561]

# Internal calculations
@parameters σ
a = y - x
eqs = [D(x) ~ σ * a,
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z * κ]
@named de = System(eqs, t)
jac = calculate_jacobian(de)
@test ModelingToolkit.jacobian_sparsity(de).colptr == sparse(jac).colptr
@test ModelingToolkit.jacobian_sparsity(de).rowval == sparse(jac).rowval

f = ODEFunction(complete(de))

@parameters A B C
_x = y / C
eqs = [D(x) ~ -A * x,
    D(y) ~ A * x - B * _x]
@named de = System(eqs, t)
@test begin
    local f, du
    f = generate_rhs(de, expression = Val{false}, wrap_gfw = Val{true})
    du = [0.0, 0.0]
    f(du, [1.0, 2.0], [1, 2, 3], 0.0)
    du ≈ [-1, -1 / 3]
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
ODEFunction(de)(similar(prob.u0), prob.u0, prob.p, 0.1)

function lotka(du, u, p, t)
    x = u[1]
    y = u[2]
    du[1] = p[1] * x - p[2] * x * y
    du[2] = -p[3] * y + p[4] * x * y
end

prob = ODEProblem(lotka, [1.0, 1.0], (0.0, 1.0), [1.5, 1.0, 3.0, 1.0])

de = complete(modelingtoolkitize(prob))
ODEFunction(de)(similar(prob.u0), prob.u0, prob.p, 0.1)

# automatic unknown detection for DAEs
@parameters k₁ k₂ k₃
@variables y₁(t) y₂(t) y₃(t)
# reorder the system just to be a little spicier
eqs = [D(y₁) ~ -k₁ * y₁ + k₃ * y₂ * y₃,
    0 ~ y₁ + y₂ + y₃ - 1,
    D(y₂) ~ k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃ * κ]
@named sys = System(eqs, t, defaults = [k₁ => 100, k₂ => 3e7, y₁ => 1.0])
sys = complete(sys)
u0 = Pair[]
push!(u0, y₂ => 0.0)
push!(u0, y₃ => 0.0)
p = [k₁ => 0.04,
    k₃ => 1e4]
p2 = [k₁ => 0.04,
    k₂ => 3e7,
    k₃ => 1e4]
tspan = (0.0, 100000.0)
prob1 = ODEProblem(sys, [u0; p], tspan)
@test prob1.f.sys == sys
prob12 = ODEProblem(sys, [u0; [k₁ => 0.04, k₂ => 3e7, k₃ => 1e4]], tspan)
prob13 = ODEProblem(sys, [u0; [k₁ => 0.04, k₂ => 3e7, k₃ => 1e4]], tspan)
prob14 = ODEProblem(sys, [u0; p2], tspan)
for p in [prob1, prob14]
    @test p.p isa MTKParameters
    p.ps[k₁] ≈ 0.04
    p.ps[k₂] ≈ 3e7
    p.ps[k₃] ≈ 1e-4
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
    @test p.p isa MTKParameters
    p.ps[k₁] ≈ 0.05
    p.ps[k₂] ≈ 2e7
    p.ps[k₃] ≈ 1.1e-4
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
        System([D(x) ~ -a * x], t; name)
    end

    function makecombinedsys()
        sys1 = makesys(:sys1)
        sys2 = makesys(:sys2)
        @parameters b = 1.0
        complete(System(Equation[], t, [], [b]; systems = [sys1, sys2], name = :foo))
    end

    sys = makecombinedsys()
    @unpack sys1, b = sys
    prob = ODEProblem(sys, Pair[], (0.0, 1.0))
    prob_new = SciMLBase.remake(prob, p = Dict(sys1.a => 3.0, b => 4.0),
        u0 = Dict(sys1.x => 1.0))
    @test prob_new.p isa MTKParameters
    @test prob_new.ps[b] ≈ 4.0
    @test prob_new.ps[sys1.a] ≈ 3.0
    @test prob_new.ps[sys.sys2.a] ≈ 1.0
    @test prob_new.u0 == [1.0, 0.0]
end

# test kwargs
prob2 = ODEProblem(sys, [u0; p], tspan, jac = true)
prob3 = ODEProblem(sys, [u0; p], tspan, jac = true, sparse = true) #SparseMatrixCSC need to handle
@test prob3.f.jac_prototype isa SparseMatrixCSC
prob3 = ODEProblem(sys, [u0; p], tspan, jac = true, sparsity = true)
@test prob3.f.sparsity isa SparseMatrixCSC
@test_throws ArgumentError ODEProblem(sys, zeros(5), tspan)
for (prob, atol) in [(prob1, 1e-12), (prob2, 1e-12), (prob3, 1e-12)]
    local sol
    sol = solve(prob, Rodas5())
    @test all(x -> ≈(sum(x), 1.0, atol = atol), sol.u)
end

du0 = [D(y₁) => -0.04
       D(y₂) => 0.04
       D(y₃) => 0.0]
prob4 = DAEProblem(sys, [du0; u0; p2], tspan)
prob5 = eval(DAEProblem(sys, [du0; u0; p2], tspan; expression = Val{true}))
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
@named sys = System(eqs, t)
@test issetequal(unknowns(sys), [x, y, z])
@test issetequal(parameters(sys), [σ, β])
@test equations(sys) == eqs
@test ModelingToolkit.isautonomous(sys)

@testset "Issue#701: `collect_vars!` handles non-call symbolics" begin
    @parameters a
    @variables x(t)
    @named sys = System([D(x) ~ a], t)
    @test issym(equations(sys)[1].rhs)
end

# issue #609
@variables x1(t) x2(t)

eqs = [
    D(x1) ~ -x1,
    0 ~ x1 - x2
]
@named sys = System(eqs, t)
@test isequal(ModelingToolkit.get_iv(sys), t)
@test isequal(unknowns(sys), [x1, x2])
@test isempty(parameters(sys))

# one equation System test
@parameters r
@variables x(t)
eq = D(x) ~ r * x
@named ode = System(eq, t)
@test equations(ode) == [eq]
# issue #808
@testset "Combined system name collisions" begin
    function makesys(name)
        @parameters a
        @variables x(t) f(t)

        System([D(x) ~ -a * x + f], t; name)
    end

    function issue808()
        sys1 = makesys(:sys1)
        sys2 = makesys(:sys1)

        @test_throws ModelingToolkit.NonUniqueSubsystemsError System(
            [sys2.f ~ sys1.x, D(sys1.f) ~ 0], t,
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
@test_throws ArgumentError System(eqs, t, vars, pars, name = :foo)

#Issue 1063/998
pars = [t]
vars = @variables((u1(t),))
@test_throws ArgumentError System(eqs, t, vars, pars, name = :foo)

@parameters w
der = Differential(w)
eqs = [
    der(u1) ~ t
]
@test_throws ArgumentError ModelingToolkit.System(eqs, t, vars, pars, name = :foo)

# check_eqs_u0 kwarg test
@variables x1(t) x2(t)
eqs = [D(x1) ~ -x1]
@named sys = System(eqs, t, [x1, x2], [])
sys = complete(sys)
@test_throws ArgumentError ODEProblem(sys, [1.0, 1.0], (0.0, 1.0))
@test_nowarn ODEProblem(sys, [1.0, 1.0], (0.0, 1.0), check_length = false)

@testset "Issue#1109" begin
    @variables x(t)[1:3, 1:3]
    @named sys = System(D(x) ~ x, t)
    @test_nowarn mtkcompile(sys)
end

# Array vars
using Symbolics: unwrap, wrap
using LinearAlgebra
sts = @variables x(t)[1:3]=[1, 2, 3.0] y(t)=1.0
ps = @parameters p[1:3] = [1, 2, 3]
eqs = [collect(D.(x) .~ x)
       D(y) ~ norm(collect(x)) * y - x[1]]
@named sys = System(eqs, t, sts, ps)
sys = mtkcompile(sys)
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
    System(D(y) ~ sum(A) * y, t; name = name)
end

# Build system
@named sys1 = submodel()
@named sys2 = submodel()

@named sys = System([0 ~ sys1.y + sys2.y], t; systems = [sys1, sys2])

# register
using StaticArrays
using SymbolicUtils: term
using SymbolicUtils.Code
using Symbolics: unwrap, wrap, @register_symbolic
foo(a, ms::AbstractVector) = a + sum(ms)
@register_symbolic foo(a, ms::AbstractVector)
@variables x(t) ms(t)[1:3]
eqs = [D(x) ~ foo(x, ms); D(ms) ~ ones(3)]
@named sys = System(eqs, t, [x; ms], [])
@named emptysys = System(Equation[], t)
@mtkcompile outersys = compose(emptysys, sys)
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
@named sys = System(eqs, t)
@named emptysys = System(Equation[], t)
@mtkcompile outersys = compose(emptysys, sys)
prob = ODEProblem(
    outersys, [sys.x => 1.0, sys.ms => 1:3, sys.p => ones(3, 3)], (0.0, 1.0))
@test_nowarn solve(prob, Tsit5())
obsfn = ModelingToolkit.build_explicit_observed_function(
    outersys, bar(3outersys.sys.ms, 3outersys.sys.p))
@test_nowarn obsfn(sol.u[1], prob.p, sol.t[1])

# x/x
@variables x(t)
@named sys = System([D(x) ~ x / x], t)
@test equations(alias_elimination(sys)) == [D(x) ~ 1]

# observed variable handling
@variables x(t) RHS(t)
@parameters τ
@named fol = System([D(x) ~ (1 - x) / τ], t; observed = [RHS ~ (1 - x) / τ])
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

@named sys = System(eqs, t, [x, y, z], [α, β])
sys = complete(sys)
@test_throws Any ODEFunction(sys)

eqs = copy(eqs)
eqs[end] = D(D(z)) ~ α * x - β * y
@named sys = System(eqs, t, [x, y, z], [α, β])
sys = complete(sys)
@test_throws Any ODEFunction(sys)

@testset "Preface tests" begin
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

    @named sys = System(eqs, t, us, ps; defaults = defs, preface = preface)
    sys = complete(sys)
    # don't build initializeprob because it will use preface in other functions and
    # affect `c`
    prob = ODEProblem(sys, [], (0.0, 1.0); build_initializeprob = false)
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
    @named sys = System(eqs, t, vcat(x, [y]), [k])
    sys = mtkcompile(sys)

    u0 = x .=> [0.5, 0]
    du0 = D.(x) .=> 0.0
    prob = DAEProblem(sys, du0, (0, 50); guesses = u0)
    @test prob[x] ≈ [0.5, 1.0]
    @test prob.du0 ≈ [0.0, 0.0]
    @test prob.p isa MTKParameters
    @test prob.ps[k] ≈ 1
    sol = solve(prob, IDA())
    @test sol[y] ≈ 0.9 * sol[x[1]] + sol[x[2]]
    @test isapprox(sol[x[1]][end], 1, atol = 1e-3)

    prob = DAEProblem(sys, [D(y) => 0, D(x[1]) => 0, D(x[2]) => 0], (0, 50); guesses = u0)

    @test prob[x] ≈ [0.5, 1]
    @test prob.du0 ≈ [0, 0]
    @test prob.p isa MTKParameters
    @test prob.ps[k] ≈ 1
    sol = solve(prob, IDA())
    @test isapprox(sol[x[1]][end], 1, atol = 1e-3)

    prob = DAEProblem(sys, [D(y) => 0, D(x[1]) => 0, D(x[2]) => 0, k => 2],
        (0, 50); guesses = u0)
    @test prob[x] ≈ [0.5, 3]
    @test prob.du0 ≈ [0, 0]
    @test prob.p isa MTKParameters
    @test prob.ps[k] ≈ 2
    sol = solve(prob, IDA())
    @test isapprox(sol[x[1]][end], 2, atol = 1e-3)

    # no initial conditions for D(x[1]) and D(x[2]) provided
    @test_throws ModelingToolkit.MissingVariablesError prob=DAEProblem(
        sys, Pair[], (0, 50); guesses = u0)

    prob = ODEProblem(sys, Pair[x[1] => 0], (0, 50))
    sol = solve(prob, Rosenbrock23())
    @test isapprox(sol[x[1]][end], 1, atol = 1e-3)
end

#issue 1475 (mixed numeric type for parameters)
let
    @parameters k1 k2::Int
    @variables A(t)
    eqs = [D(A) ~ -k1 * k2 * A]
    @named sys = System(eqs, t)
    sys = complete(sys)
    ivmap = [A => 1.0, k1 => 1.0, k2 => 1.0]
    tspan = (0.0, 1.0)
    prob = ODEProblem(sys, ivmap, tspan; tofloat = false)
    @test prob.p isa MTKParameters
    @test prob.ps[k1] ≈ 1.0
    @test prob.ps[k2] == 1 && prob.ps[k2] isa Int
end

let
    @parameters C L R
    @variables q(t) p(t) F(t)

    eqs = [D(q) ~ -p / L - F
           D(p) ~ q / C
           0 ~ q / C - R * F]

    @named sys = System(eqs, t)
    @test length(equations(mtkcompile(sys))) == 2
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
    @named sys = System(eqs, t, vars, pars)
    @test_throws ModelingToolkit.ExtraEquationsSystemException mtkcompile(sys)
end

# 1561
let
    vars = @variables x y
    arr = ModelingToolkit.varmap_to_vars(
        Dict([x => 0.0, y => [0.0, 1.0]]), vars; use_union = true) #error
    sol = Union{Float64, Vector{Float64}}[0.0, [
        0.0, 1.0]]
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

    @named sys = System(eqs, t, u, ps)
    @test_nowarn simpsys = mtkcompile(sys)

    sys = mtkcompile(sys)

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
    @named osys = System(eqs, t)
    osys = complete(osys)
    oprob = ODEProblem(osys, [A => 1.0, k => 1.0], (0.0, 10.0); check_length = false)
    @test_nowarn sol = solve(oprob, Tsit5())
end

let
    function sys1(; name)
        vars = @variables x(t)=0.0 dx(t)=0.0

        System([D(x) ~ dx], t, vars, []; name, defaults = [D(x) => x])
    end

    function sys2(; name)
        @named s1 = sys1()

        System(Equation[], t, [], []; systems = [s1], name)
    end

    s1′ = sys1(; name = :s1)
    @named s2 = sys2()
    @unpack s1 = s2
    @test isequal(unknowns(s1), unknowns(s1′))
    @test isequal(parameters(s1), parameters(s1′))
    @test isequal(equations(s1), equations(s1′))

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
        System(eqs, t; name)
    end

    @named ctrl = pd_ctrl()

    ## double integrator

    function double_int(; name)
        @variables u(t) x(t) v(t)

        eqs = [D(x) ~ v, D(v) ~ u]
        System(eqs, t; name)
    end

    @named sys = double_int()

    ## connections

    connections = [sys.u ~ ctrl.u, ctrl.x ~ sys.x, ctrl.v ~ sys.v]

    @named connected = System(connections, t)
    @named sys_con = compose(connected, sys, ctrl)

    sys_simp = mtkcompile(sys_con)
    true_eqs = [D(sys.x) ~ sys.v
                D(sys.v) ~ ctrl.kv * sys.v + ctrl.kx * sys.x]
    @test issetequal(full_equations(sys_simp), true_eqs)
end

let
    @variables x(t) = 1
    @variables y(t) = 1
    @parameters pp = -1
    @named sys4 = System([D(x) ~ -y; D(y) ~ 1 + pp * y + x], t)
    sys4s = mtkcompile(sys4)
    prob = ODEProblem(sys4s, [x => 1.0, D(x) => 1.0], (0, 1.0))
    @test issetequal(string.(unknowns(prob.f.sys)), ["x(t)", "y(t)"])
    @test string.(parameters(prob.f.sys)) == ["pp"]
    @test string.(independent_variables(prob.f.sys)) == ["t"]
end

@variables P(t)=NaN Q(t)=NaN
eqs = [D(Q) ~ 1 / sin(P), D(P) ~ log(-cos(Q))]
@named sys = System(eqs, t, [P, Q], [])
sys = complete(debug_system(sys))
prob = ODEProblem(sys, [], (0.0, 1.0))
@test_throws "log(-cos(Q(t))) errors" prob.f([1, 0], prob.p, 0.0)
@test_throws "/(1, sin(P(t))) output non-finite value" prob.f([0, 2], prob.p, 0.0)

let
    @variables x(t) = 1
    @variables y(t) = 1
    @parameters pp = -1
    der = Differential(t)
    @named sys4 = System([der(x) ~ -y; der(y) ~ 1 + pp * y + x], t)
    sys4s = mtkcompile(sys4)
    prob = ODEProblem(sys4s, [x => 1.0, D(x) => 1.0], (0, 1.0))
    @test !isnothing(prob.f.sys)
end

# SYS 1:
vars_sub1 = @variables s1(t)
@named sub = System(Equation[], t, vars_sub1, [])

vars1 = @variables x1(t)
@named sys1 = System(Equation[], t, vars1, [], systems = [sub])
@named sys2 = System(Equation[], t, vars1, [], systems = [sys1, sub])

# SYS 2: Extension to SYS 1
vars_sub2 = @variables s2(t)
@named partial_sub = System(Equation[], t, vars_sub2, [])
@named sub = extend(partial_sub, sub)

# no warnings for systems without events
new_sys2 = @test_nowarn complete(substitute(sys2, Dict(:sub => sub)))
Set(unknowns(new_sys2)) == Set([new_sys2.x1, new_sys2.sys1.x1,
    new_sys2.sys1.sub.s1, new_sys2.sys1.sub.s2,
    new_sys2.sub.s1, new_sys2.sub.s2])

let # Issue https://github.com/SciML/ModelingToolkit.jl/issues/2322
    @parameters a=10 b=a/10 c=a/20

    Dt = D

    @variables x(t)=1 z(t)

    eqs = [Dt(x) ~ -b * (x - z),
        0 ~ z - c * x]

    sys = System(eqs, t; name = :kjshdf)

    sys_simp = mtkcompile(sys)

    @test a ∈ keys(ModelingToolkit.defaults(sys_simp))

    tspan = (0.0, 1)
    prob = ODEProblem(sys_simp, [], tspan)
    sol = solve(prob, Rodas4())
    @test sol(1)[]≈0.6065307685451087 rtol=1e-4
end

# Issue#2599
@variables x(t) y(t)
eqs = [D(x) ~ x * t, y ~ 2x]
@mtkcompile sys = System(eqs, t; continuous_events = [[y ~ 3] => [x ~ 2]])
prob = ODEProblem(sys, [x => 1.0], (0.0, 10.0))
@test_nowarn solve(prob, Tsit5())

# Issue#2383
@testset "Arrays in affect/condition equations" begin
    @variables x(t)[1:3]
    @parameters p[1:3, 1:3]
    eqs = [
        D(x) ~ p * x
    ]
    @mtkcompile sys = System(
        eqs, t; continuous_events = [[norm(x) ~ 3.0] => [x ~ ones(3)]])
    # array affect equations used to not work
    prob1 = @test_nowarn ODEProblem(sys, [x => ones(3), p => ones(3, 3)], (0.0, 10.0))
    sol1 = @test_nowarn solve(prob1, Tsit5())

    # array condition equations also used to not work
    @mtkcompile sys = System(
        eqs, t; continuous_events = [[x ~ sqrt(3) * ones(3)] => [x ~ ones(3)]])
    # array affect equations used to not work
    prob2 = @test_nowarn ODEProblem(sys, [x => ones(3), p => ones(3, 3)], (0.0, 10.0))
    sol2 = @test_nowarn solve(prob2, Tsit5())

    @test sol1.u ≈ sol2.u
end

# Requires fix in symbolics for `linear_expansion(p * x, D(y))`
@test_skip begin
    @variables x(t)[1:3] y(t)
    @parameters p[1:3, 1:3]
    @test_nowarn @mtkcompile sys = System([D(x) ~ p * x, D(y) ~ x' * p * x], t)
    @test_nowarn ODEProblem(sys, [x => ones(3), y => 2, p => ones(3, 3)], (0.0, 10.0))
end

@parameters g L
@variables q₁(t) q₂(t) λ(t) θ(t)

eqs = [D(D(q₁)) ~ -λ * q₁,
    D(D(q₂)) ~ -λ * q₂ - g,
    q₁ ~ L * sin(θ),
    q₂ ~ L * cos(θ)]

@named pend = System(eqs, t)
pend = complete(pend)
@test_nowarn generate_initializesystem(
    pend; op = [q₁ => 1.0, q₂ => 0.0], guesses = [λ => 1])

# https://github.com/SciML/ModelingToolkit.jl/issues/2618
@parameters σ ρ β
@variables x(t) y(t) z(t)

eqs = [D(D(x)) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@mtkcompile sys = System(eqs, t)

u0 = [D(x) => 2.0,
    x => 1.0,
    y => 0.0,
    z => 0.0]

p = [σ => 28.0,
    ρ => 10.0,
    β => 8 / 3]

prob = SteadyStateProblem(sys, [u0; p])
@test prob isa SteadyStateProblem
prob = SteadyStateProblem(ODEProblem(sys, [u0; p], (0.0, 10.0)))
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
    System(eqs, t; systems, name)
end

@mtkcompile model = FML2()

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
    sys = System(eqns, t, vars, []; name)
end

function RealExpressionSystem(; name)
    vars = @variables begin
        x(t)
        z(t)[1:1]
    end # doing a collect on z doesn't work either. 
    @named e1 = RealExpression(y = x) # This works perfectly. 
    @named e2 = RealExpression(y = z[1]) # This bugs. However, `full_equations(e2)` works as expected. 
    systems = [e1, e2]
    System(Equation[], t, Iterators.flatten(vars), []; systems, name)
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
@named sys = System(Equation[], t, [x], []; guesses = [x => 1.0])
@named outer = System(
    [D(y) ~ sys.x + t, 0 ~ t + y - sys.x * y], t, [y], []; systems = [sys])
@test ModelingToolkit.guesses(outer)[sys.x] == 1.0
outer = mtkcompile(outer)
@test ModelingToolkit.get_guesses(outer)[sys.x] == 1.0
prob = ODEProblem(outer, [outer.y => 2.0], (0.0, 10.0))
int = init(prob, Rodas4())
@test int[outer.sys.x] == 1.0

# Ensure indexes of array symbolics are cached appropriately
@variables x(t)[1:2]
@named sys = System(Equation[], t, [x], [])
sys1 = complete(sys)
@named sys = System(Equation[], t, [x...], [])
sys2 = complete(sys)
for sys in [sys1, sys2]
    for (sym, idx) in [(x, 1:2), (x[1], 1), (x[2], 2)]
        @test is_variable(sys, sym)
        @test variable_index(sys, sym) == idx
    end
end

@variables x(t)[1:2, 1:2]
@named sys = System(Equation[], t, [x], [])
sys1 = complete(sys)
@named sys = System(Equation[], t, [x...], [])
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
    @named sys = System([x[0] ~ 0.0, D(x[1]) ~ x[0]], t, [x], [])
    @test_nowarn sys = mtkcompile(sys)
    @test equations(sys) == [D(x[1]) ~ 0.0]
end

# Namespacing of array variables
@testset "Namespacing of array variables" begin
    @variables x(t)[1:2]
    @named sys = System(Equation[], t)
    @test getname(unknowns(sys, x)) == :sys₊x
    @test size(unknowns(sys, x)) == size(x)
end

# Issue#2667 and Issue#2953
@testset "ForwardDiff through ODEProblem constructor" begin
    @parameters P
    @variables x(t)
    sys = mtkcompile(System([D(x) ~ P], t, [x], [P]; name = :sys))

    function x_at_1(P)
        prob = ODEProblem(sys, [x => P, sys.P => P], (0.0, 1.0))
        return solve(prob, Tsit5())(1.0)
    end

    @test_nowarn ForwardDiff.derivative(P -> x_at_1(P), 1.0)
end

@testset "Inplace observed functions" begin
    @parameters P
    @variables x(t)
    sys = mtkcompile(System([D(x) ~ P], t, [x], [P]; name = :sys))
    obsfn = ModelingToolkit.build_explicit_observed_function(
        sys, [x + 1, x + P, x + t], return_inplace = true)[2]
    ps = ModelingToolkit.MTKParameters(sys, [P => 2.0])
    buffer = zeros(3)
    @test_nowarn obsfn(buffer, [1.0], ps, 3.0)
    @test buffer ≈ [2.0, 3.0, 4.0]
end

# https://github.com/SciML/ModelingToolkit.jl/issues/2818
@testset "Custom independent variable" begin
    @independent_variables x
    @variables y(x)
    @test_nowarn @named sys = System([y ~ 0], x)

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
    @test_nowarn @mtkcompile sys = MyModel()

    @variables x y(x)
    @test_logs (:warn,) @named sys = System([y ~ 0], x)

    @parameters T
    D = Differential(T)
    @variables x(T)
    eqs = [D(x) ~ 0.0]
    initialization_eqs = [x ~ T]
    guesses = [x => 0.0]
    @named sys2 = System(eqs, T; initialization_eqs, guesses)
    prob2 = ODEProblem(mtkcompile(sys2), [], (1.0, 2.0))
    sol2 = solve(prob2)
    @test all(sol2[x] .== 1.0)
end

# https://github.com/SciML/ModelingToolkit.jl/issues/2502
@testset "Extend systems with a field that can be nothing" begin
    A = Dict(Int => 1)
    B = Dict(String => 2)
    @named A1 = System(Equation[], t, [], [])
    @named B1 = System(Equation[], t, [], [])
    @named A2 = System(Equation[], t, [], []; metadata = A)
    @named B2 = System(Equation[], t, [], []; metadata = B)
    n_core_metadata = length(ModelingToolkit.get_metadata(A1))
    @test length(ModelingToolkit.get_metadata(extend(A1, B1))) == n_core_metadata
    meta = ModelingToolkit.get_metadata(extend(A1, B2))
    @test length(meta) == n_core_metadata + 1
    @test meta[String] == 2
    meta = ModelingToolkit.get_metadata(extend(A2, B1))
    @test length(meta) == n_core_metadata + 1
    @test meta[Int] == 1
    meta = ModelingToolkit.get_metadata(extend(A2, B2))
    @test length(meta) == n_core_metadata + 2
    @test meta[Int] == 1
    @test meta[String] == 2
end

# https://github.com/SciML/ModelingToolkit.jl/issues/2859
@testset "Initialization with defaults from observed equations (edge case)" begin
    @variables x(t) y(t) z(t)
    eqs = [D(x) ~ 0, y ~ x, D(z) ~ 0]
    defaults = [x => 1, z => y]
    @named sys = System(eqs, t; defaults)
    ssys = mtkcompile(sys)
    prob = ODEProblem(ssys, [], (0.0, 1.0))
    @test prob[x] == prob[y] == prob[z] == 1.0

    @parameters y0
    @variables x(t) y(t) z(t)
    eqs = [D(x) ~ 0, y ~ y0 / x, D(z) ~ y]
    defaults = [y0 => 1, x => 1, z => y]
    @named sys = System(eqs, t; defaults)
    ssys = mtkcompile(sys)
    prob = ODEProblem(ssys, [], (0.0, 1.0))
    @test prob[x] == prob[y] == prob[z] == 1.0
end

@testset "Scalarized parameters in array functions" begin
    @variables u(t)[1:2] x(t)[1:2] o(t)[1:2]
    @parameters p[1:2, 1:2] [tunable = false]
    @named sys = System(
        [D(u) ~ (sum(u) + sum(x) + sum(p) + sum(o)) * x, o ~ prod(u) * x],
        t, [u..., x..., o...], [p...])
    sys1 = mtkcompile(sys, inputs = [x...], outputs = [])
    fn1, = ModelingToolkit.generate_rhs(sys1; expression = Val{false})
    ps = MTKParameters(sys1, [x => 2ones(2), p => 3ones(2, 2)])
    @test_nowarn fn1(ones(4), ps, 4.0)
    sys2 = mtkcompile(sys, inputs = [x...], outputs = [], split = false)
    fn2, = ModelingToolkit.generate_rhs(sys2; expression = Val{false})
    ps = zeros(8)
    setp(sys2, x)(ps, 2ones(2))
    setp(sys2, p)(ps, 2ones(2, 2))
    @test_nowarn fn2(ones(4), 2ones(14), 4.0)
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
    prob = ODEProblem(sys, merge(p, Dict(u0)), (0.0, 1.0))

    # evaluate
    u0_v = prob.u0
    p_v = prob.p
    @test prob.f(u0_v, p_v, 0.0) == [c_b, c_a]
end

@testset "Independent variable as system property" begin
    @variables x(t)
    @named sys = System([x ~ t], t)
    @named sys = compose(sys, sys) # nest into a hierarchical system
    @test t === sys.t === sys.sys.t
end

@testset "Substituting preserves parameter dependencies, defaults, guesses" begin
    @parameters p1 p2
    @variables x(t) y(t)
    @named sys = System([D(x) ~ y + p2, p2 ~ 2p1], t;
        defaults = [p1 => 1.0, p2 => 2.0], guesses = [p1 => 2.0, p2 => 3.0])
    @parameters p3
    sys2 = substitute(sys, [p1 => p3])
    sys2 = complete(sys2)
    @test length(parameters(sys2)) == 1
    @test is_parameter(sys2, p3)
    @test !is_parameter(sys2, p1)
    @test length(ModelingToolkit.defaults(sys2)) == 7
    @test ModelingToolkit.defaults(sys2)[p3] == 1.0
    @test length(ModelingToolkit.guesses(sys2)) == 2
    @test ModelingToolkit.guesses(sys2)[p3] == 2.0
end

@testset "Substituting with nested systems" begin
    @parameters p1 p2
    @variables x(t) y(t)
    @named innersys = System([D(x) ~ y + p2; p2 ~ 2p1], t;
        defaults = [p1 => 1.0, p2 => 2.0], guesses = [p1 => 2.0, p2 => 3.0])
    @parameters p3 p4
    @named outersys = System(
        [D(innersys.y) ~ innersys.y + p4, p4 ~ 3p3], t;
        defaults = [p3 => 3.0, p4 => 9.0], guesses = [p4 => 10.0], systems = [innersys])
    @test_nowarn mtkcompile(outersys)
    @parameters p5
    sys2 = substitute(outersys, [p4 => p5])
    sys2 = complete(sys2)
    @test_nowarn mtkcompile(sys2)
    @test length(equations(sys2)) == 2
    @test length(parameters(sys2)) == 2
    @test length(full_parameters(sys2)) == 10
    @test all(!isequal(p4), full_parameters(sys2))
    @test any(isequal(p5), full_parameters(sys2))
    @test length(ModelingToolkit.defaults(sys2)) == 10
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

    @named sys = System(eqs, t, [u..., x..., o], [p...])
    sys1 = mtkcompile(sys, inputs = [x...], outputs = [o...], split = false)

    @test_nowarn ModelingToolkit.build_explicit_observed_function(sys1, u; inputs = [x...])

    obsfn = ModelingToolkit.build_explicit_observed_function(
        sys1, u + x + p[1:2]; inputs = [x...])

    @test obsfn(ones(2), 2ones(2), 3ones(12), 4.0) == 6ones(2)
end

@testset "Passing `nothing` to `u0`" begin
    @variables x(t) = 1
    @mtkcompile sys = System(D(x) ~ t, t)
    prob = @test_nowarn ODEProblem(sys, nothing, (0.0, 1.0))
    @test_nowarn solve(prob)
end

@testset "ODEs are not DDEs" begin
    @variables x(t)
    @named sys = System(D(x) ~ x, t)
    @test !ModelingToolkit.is_dde(sys)
    @test is_markovian(sys)
    @named sys2 = System(Equation[], t; systems = [sys])
    @test !ModelingToolkit.is_dde(sys)
    @test is_markovian(sys)
end

@testset "Issue #2597" begin
    @variables x(t)[1:2]=ones(2) y(t)=1.0

    for eqs in [D(x) ~ x, collect(D(x) .~ x)]
        for dvs in [[x], collect(x)]
            @named sys = System(eqs, t, dvs, [])
            sys = complete(sys)
            if eqs isa Vector && length(eqs) == 2 && length(dvs) == 2
                @test_nowarn ODEProblem(sys, [], (0.0, 1.0))
            else
                @test_throws [
                    r"array (equations|unknowns)", "mtkcompile", "scalarize"] ODEProblem(
                    sys, [], (0.0, 1.0))
            end
        end
    end
    for eqs in [[D(x) ~ x, D(y) ~ y], [collect(D(x) .~ x); D(y) ~ y]]
        for dvs in [[x, y], [x..., y]]
            @named sys = System(eqs, t, dvs, [])
            sys = complete(sys)
            if eqs isa Vector && length(eqs) == 3 && length(dvs) == 3
                @test_nowarn ODEProblem(sys, [], (0.0, 1.0))
            else
                @test_throws [
                    r"array (equations|unknowns)", "mtkcompile", "scalarize"] ODEProblem(
                    sys, [], (0.0, 1.0))
            end
        end
    end
end

@testset "Parameter dependencies with constant RHS" begin
    @parameters p
    @test_nowarn System([p ~ 1.0], t; name = :a)
end

@testset "Variable discovery in arrays of `Num` inside callable symbolic" begin
    @variables x(t) y(t)
    @parameters foo(::AbstractVector)
    sys = @test_nowarn System(D(x) ~ foo([x, 2y]), t; name = :sys)
    @test length(unknowns(sys)) == 2
    @test any(isequal(y), unknowns(sys))
end

@testset "Inplace observed" begin
    @variables x(t)
    @parameters p[1:2] q
    @mtkcompile sys = System(D(x) ~ sum(p) * x + q * t, t)
    prob = ODEProblem(sys, [x => 1.0, p => ones(2), q => 2], (0.0, 1.0))
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
    @mtkcompile sys = System(D(x) ~ p * t, t)
    @test ModelingToolkit.get_index_cache(sys) !== nothing
    sys2 = complete(sys; split = false)
    @test ModelingToolkit.get_index_cache(sys2) === nothing
end

# https://github.com/SciML/SciMLBase.jl/issues/786
@testset "Observed variables dependent on discrete parameters" begin
    @variables x(t) obs(t)
    @parameters c(t)
    @mtkcompile sys = System([D(x) ~ c * cos(x), obs ~ c],
        t,
        [x, obs],
        [c];
        discrete_events = [SymbolicDiscreteCallback(
            1.0 => [c ~ Pre(c) + 1], discrete_parameters = [c])])
    prob = ODEProblem(sys, [x => 0.0, c => 1.0], (0.0, 2pi))
    sol = solve(prob, Tsit5())
    @test sol[obs] ≈ 1:7
end

@testset "DAEProblem with array parameters" begin
    @variables x(t)=1.0 y(t) [guess = 1.0]
    @parameters p[1:2] = [1.0, 2.0]
    @mtkcompile sys = System([D(x) ~ x, y^2 ~ x + sum(p)], t)
    prob = DAEProblem(sys, [D(x) => x, D(y) => D(x) / 2y], (0.0, 1.0))
    sol = solve(prob, DFBDF(), abstol = 1e-8, reltol = 1e-8)
    @test sol[x]≈sol[y ^ 2 - sum(p)] atol=1e-5
end

@testset "Symbolic tstops" begin
    @variables x(t) = 1.0
    @parameters p=0.15 q=0.25 r[1:2]=[0.35, 0.45]
    @mtkcompile sys = System(
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
    @mtkcompile sys = System([D(x) ~ p * x + q * t + sum(r), y^3 ~ 2x + 1],
        t; tstops = [0.5p, [0.1, 0.2], [p + 2q], r])
    prob = DAEProblem(
        sys, [D(y) => 2D(x) / 3y^2, D(x) => p * x + q * t + sum(r)], (0.0, 5.0))
    sol = solve(prob, DImplicitEuler())
    expected_tstops = unique!(sort!(vcat(0.0:0.075:5.0, 0.1, 0.2, 0.65, 0.35, 0.45)))
    @test all(x -> any(isapprox(x, atol = 1e-6), sol.t), expected_tstops)
    prob2 = remake(prob; tspan = (0.0, 10.0))
    sol2 = solve(prob2, DImplicitEuler())
    expected_tstops = unique!(sort!(vcat(0.0:0.075:10.0, 0.1, 0.2, 0.65, 0.35, 0.45)))
    @test all(x -> any(isapprox(x, atol = 1e-6), sol2.t), expected_tstops)

    @mtkcompile sys = System([D(x) ~ x + p], t; tstops = [[p]])
    prob = ODEProblem(sys, [], (0.0, 1.0))
    @test prob.kwargs[:tstops](prob.p, prob.tspan) ≈ [0.15]
end

@testset "Validate input types" begin
    @parameters p d
    @variables X(t)::Int64
    eq = D(X) ~ p - d * X
    @test_throws ModelingToolkit.ContinuousOperatorDiscreteArgumentError @mtkcompile osys = System(
        [eq], t)
    @variables Y(t)[1:3]::String
    eq = D(Y) ~ [p, p, p]
    @test_throws ModelingToolkit.ContinuousOperatorDiscreteArgumentError @mtkcompile osys = System(
        [eq], t)

    @variables X(t)::Complex
    eqs = D(X) ~ p - d * X
    @test_nowarn @named osys = System(eqs, t)
end

@testset "Constraint system construction" begin
    @variables x(..) y(..) z(..)
    @parameters a b c d e
    eqs = [D(x(t)) ~ 3 * a * y(t), D(y(t)) ~ x(t) - z(t), D(z(t)) ~ e * x(t)^2]
    cons = [x(0.3) ~ c * d, y(0.7) ~ 3]

    # Test variables + parameters infer correctly.
    @mtkcompile sys = System(eqs, t; constraints = cons)
    @test issetequal(parameters(sys), [a, c, d, e])
    @test issetequal(unknowns(sys), [x(t), y(t), z(t)])

    @parameters t_c
    cons = [x(t_c) ~ 3]
    @mtkcompile sys = System(eqs, t; constraints = cons)
    @test issetequal(parameters(sys), [a, e, t_c])

    @parameters g(..) h i
    cons = [g(h, i) * x(3) ~ c]
    @mtkcompile sys = System(eqs, t; constraints = cons)
    @test issetequal(parameters(sys), [g, h, i, a, e, c])

    # Test that bad constraints throw errors.
    cons = [x(3, 4) ~ 3] # unknowns cannot have multiple args.
    @test_throws ArgumentError @mtkcompile sys = System(eqs, t; constraints = cons)

    cons = [x(y(t)) ~ 2] # unknown arg must be parameter, value, or t
    @test_throws ArgumentError @mtkcompile sys = System(eqs, t; constraints = cons)

    @variables u(t) v
    cons = [x(t) * u ~ 3]
    @test_throws ArgumentError @mtkcompile sys = System(eqs, t; constraints = cons)
    cons = [x(t) * v ~ 3]
    @test_throws ArgumentError @mtkcompile sys = System(eqs, t; constraints = cons) # Need time argument.

    # Test array variables
    @variables x(..)[1:5]
    mat = [1 2 0 3 2
           0 0 3 2 0
           0 1 3 0 4
           2 0 0 2 1
           0 0 2 0 5]
    eqs = D(x(t)) ~ mat * x(t)
    cons = [x(3) ~ [2, 3, 3, 5, 4]]
    @mtkcompile ode = System(D(x(t)) ~ mat * x(t), t; constraints = cons)
    @test length(constraints(ode)) == 1
end

@testset "`build_explicit_observed_function` with `expression = true` returns `Expr`" begin
    @variables x(t)
    @mtkcompile sys = System(D(x) ~ 2x, t)
    obsfn_expr = ModelingToolkit.build_explicit_observed_function(
        sys, 2x + 1, expression = true)
    @test obsfn_expr isa Expr
    obsfn_expr_oop,
    obsfn_expr_iip = ModelingToolkit.build_explicit_observed_function(
        sys, [x + 1, x + 2, x + t], return_inplace = true, expression = true)
    @test obsfn_expr_oop isa Expr
    @test obsfn_expr_iip isa Expr
end

@testset "Solve with `split=false` static arrays" begin
    @parameters σ ρ β
    @variables x(t) y(t) z(t)

    eqs = [D(D(x)) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z]

    @mtkcompile sys=System(eqs, t) split=false

    u0 = SA[D(x) => 2.0f0,
    x => 1.0f0,
    y => 0.0f0,
    z => 0.0f0]

    p = SA[σ => 28.0f0,
    ρ => 10.0f0,
    β => 8.0f0 / 3.0f0]

    tspan = (0.0f0, 100.0f0)
    prob = ODEProblem{false}(sys, [u0; p], tspan)
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
end

@testset "`@named` always wraps in `ParentScope`" begin
    function SysA(; name, var1)
        @variables x(t)
        scope = ModelingToolkit.getmetadata(unwrap(var1), ModelingToolkit.SymScope, nothing)
        @test scope isa ParentScope
        @test scope.parent isa ParentScope
        @test scope.parent.parent isa LocalScope
        return System(D(x) ~ var1, t; name)
    end
    function SysB(; name, var1)
        @variables x(t)
        @named subsys = SysA(; var1)
        return System(D(x) ~ x, t; systems = [subsys], name)
    end
    function SysC(; name)
        @variables x(t)
        @named subsys = SysB(; var1 = x)
        return System(D(x) ~ x, t; systems = [subsys], name)
    end
    @mtkcompile sys = SysC()
    @test length(unknowns(sys)) == 3
end

@testset "`full_equations` doesn't recurse infinitely" begin
    code = """
    using ModelingToolkit
    using ModelingToolkit: t_nounits as t, D_nounits as D
    @variables x(t)[1:3]=[0,0,1]
    @variables u1(t)=0 u2(t)=0 
    y₁, y₂, y₃ = x
    k₁, k₂, k₃ = 1,1,1
    eqs = [
        D(y₁) ~ -k₁*y₁ + k₃*y₂*y₃ + u1
        D(y₂) ~ k₁*y₁ - k₃*y₂*y₃ - k₂*y₂^2 + u2
        y₁ + y₂ + y₃ ~ 1
    ]

    @named sys = System(eqs, t)

    inputs = [u1, u2]
    outputs = [y₁, y₂, y₃]
    ss = mtkcompile(sys; inputs)
    full_equations(ss)
    """

    cmd = `$(Base.julia_cmd()) --project=$(@__DIR__) -e $code`
    proc = run(cmd, stdin, stdout, stderr; wait = false)
    sleep(180)
    @test !process_running(proc)
    kill(proc, Base.SIGKILL)
end

@testset "`ProblemTypeCtx`" begin
    @variables x(t)
    @mtkcompile sys = System(
        [D(x) ~ x], t; metadata = [ModelingToolkit.ProblemTypeCtx => "A"])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0))
    @test prob.problem_type == "A"
end

@testset "`substitute` retains events and metadata" begin
    @parameters p(t) = 1.0
    @variables x(t) = 0.0
    event = [0.5] => [p ~ Pre(t)]
    event2 = [x ~ 0.75] => [p ~ 2 * Pre(t)]

    struct TestMeta end

    eq = [
        D(x) ~ p
    ]
    @named sys = System(eq, t, [x], [p], discrete_events = [event],
        continuous_events = [event2], metadata = Dict(TestMeta => "test"))

    @variables x2(t) = 0.0
    sys2 = substitute(sys, [x => x2])

    @test length(ModelingToolkit.get_discrete_events(sys)) == 1
    @test length(ModelingToolkit.get_discrete_events(sys2)) == 1
    @test length(ModelingToolkit.get_continuous_events(sys)) == 1
    @test length(ModelingToolkit.get_continuous_events(sys2)) == 1
    @test getmetadata(sys, TestMeta, nothing) == "test"
    @test getmetadata(sys2, TestMeta, nothing) == "test"
end

struct TestWrapper
    sys::ODESystem
end

@testset "`ODESystem` is a type" begin
    @variables x(t)
    @named sys = ODESystem(D(x) ~ x, t)
    @test sys isa ODESystem
    @test sys isa System
    arr = ODESystem[]
    @test_nowarn push!(arr, sys)
    @test_nowarn TestWrapper(sys)
end

# ensure `@mtkbuild` works when `@mtkcompile` is not imported
module MtkbuildTestModule
import ModelingToolkit: @variables, System, t_nounits as t, D_nounits as D, @mtkbuild
import Test: @test
@variables x(t)
@mtkbuild sys = System(D(x) ~ t, t)
@test sys isa System
end

@testset "Empty system can be simplified" begin
    @named sys = System(Equation[], t)
    ss = mtkcompile(sys)
    @test length(equations(ss)) == length(unknowns(ss)) == 0
end
