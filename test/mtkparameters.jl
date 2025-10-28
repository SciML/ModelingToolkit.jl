using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D, MTKParameters
using SymbolicIndexingInterface, StaticArrays
using SciMLStructures: SciMLStructures, canonicalize, Tunable, Discrete, Constants
using ModelingToolkitStandardLibrary.Electrical, ModelingToolkitStandardLibrary.Blocks
using BlockArrays: BlockedArray, BlockedVector, Block
using OrdinaryDiffEq
using ForwardDiff
using JET

@parameters a b c(t) d::Integer e[1:3] f[1:3, 1:3]::Int g::Vector{AbstractFloat} h::String
@named sys = System(
    [b ~ 2a], t, [], [a, b, c, d, e, f, g, h];
    continuous_events = [ModelingToolkit.SymbolicContinuousCallback(
        [a ~ 0] => [c ~ 0], discrete_parameters = c)], defaults = Dict(a => 0.0))
sys = complete(sys)

ivs = Dict(c => 3a, d => 4, e => [5.0, 6.0, 7.0],
    f => ones(Int, 3, 3), g => [0.1, 0.2, 0.3], h => "foo")

ps = MTKParameters(sys, ivs)
@test_nowarn copy(ps)
ps_copy = copy(ps)
ps_field_equals = map(fieldnames(typeof(ps))) do f
    getfield(ps, f) == getfield(ps_copy, f)
end
@test all(ps_field_equals)
# dependent initialization, also using defaults
@test getp(sys, a)(ps) == getp(sys, b)(ps) == getp(sys, c)(ps) == 0.0
@test getp(sys, d)(ps) isa Int

@testset "`p_constructor`" begin
    ps2 = MTKParameters(sys, ivs; p_constructor = x -> SArray{Tuple{size(x)...}}(x))
    @test ps2.tunable isa SVector
    @test ps2.initials isa SVector
    @test ps2.discrete isa Tuple{<:BlockedVector{Float64, <:SVector}}
    @test ps2.constant isa Tuple{<:SVector, <:SVector, <:SVector{1, <:SMatrix}}
    @test ps2.nonnumeric isa Tuple{<:SVector}
end

ivs[a] = 1.0
ps = MTKParameters(sys, ivs)
for (p, val) in ivs
    if isequal(p, c)
        val = 3ivs[a]
    end
    idx = parameter_index(sys, p)
    # ensure getindex with `ParameterIndex` works
    @test ps[idx] == getp(sys, p)(ps) == val
end

# ensure setindex! with `ParameterIndex` works
ps[parameter_index(sys, a)] = 3.0
@test getp(sys, a)(ps) == 3.0
setp(sys, a)(ps, 1.0)

@test getp(sys, a)(ps) == getp(sys, b)(ps) / 2 == getp(sys, c)(ps) / 3 == 1.0

for (portion, values) in [(Tunable(), [1.0, 5.0, 6.0, 7.0])
     (Discrete(), [3.0])
     (Constants(), vcat([0.1, 0.2, 0.3], ones(9), [4.0]))]
    buffer, repack, alias = canonicalize(portion, ps)
    @test alias
    @test sort(collect(buffer)) == values
    @test all(isone,
        canonicalize(portion, SciMLStructures.replace(portion, ps, ones(length(buffer))))[1])
    # make sure it is out-of-place
    @test sort(collect(buffer)) == values
    SciMLStructures.replace!(portion, ps, ones(length(buffer)))
    # make sure it is in-place
    @test all(isone, canonicalize(portion, ps)[1])
    global ps = repack(zeros(length(buffer)))
    @test all(iszero, canonicalize(portion, ps)[1])
end

setp(sys, a)(ps, 2.0) # test set_parameter!
@test getp(sys, a)(ps) == 2.0

setp(sys, e)(ps, 5ones(3)) # with an array
@test getp(sys, e)(ps) == 5ones(3)

setp(sys, f[2, 2])(ps, 42) # with a sub-index
@test getp(sys, f[2, 2])(ps) == 42

setp(sys, g)(ps, ones(100)) # with non-fixed-length array
@test getp(sys, g)(ps) == ones(100)

setp(sys, h)(ps, "bar") # with a non-numeric
@test getp(sys, h)(ps) == "bar"

varmap = Dict(a => 1.0f0, b => 5.0f0, c => 2.0, d => 0x5, e => Float32[0.4, 0.5, 0.6],
    f => 3ones(UInt, 3, 3), g => ones(Float32, 4), h => "bar")
@test_deprecated remake_buffer(sys, ps, varmap)
@test_warn ["Symbolic variable b", "non-dependent", "parameter"] remake_buffer(
    sys, ps, keys(varmap), values(varmap))
newps = remake_buffer(sys, ps, keys(varmap), values(varmap))

for fname in (:tunable, :discrete, :constant)
    # ensure same number of sub-buffers
    @test length(getfield(ps, fname)) == length(getfield(newps, fname))
end

@test getp(sys, a)(newps) isa Float32
@test getp(sys, b)(newps) == 2.0f0 # ensure dependent update still happened, despite explicit value
@test getp(sys, c)(newps) isa Float64
@test getp(sys, d)(newps) isa UInt8
@test getp(sys, f)(newps) isa Matrix{UInt}
@test getp(sys, g)(newps) isa Vector{Float32}

@testset "Type-stability of `remake_buffer`" begin
    prob = ODEProblem(sys, ivs, (0.0, 1.0))

    idxs = (a, c, d, e, f, g, h)
    vals = (1.0, 2.0, 3, ones(3), ones(Int, 3, 3), ones(2), "a")

    setter = setsym_oop(prob, idxs)
    @test_nowarn @inferred setter(prob, vals)
    @test_throws ErrorException @inferred setter(prob, collect(vals))

    idxs = (a, c, e...)
    vals = Float16[1.0, 2.0, 3.0, 4.0, 5.0]
    setter = setsym_oop(prob, idxs)
    @test_nowarn @inferred setter(prob, vals)

    idxs = [a, e]
    vals = (Float16(1.0), ForwardDiff.Dual{Nothing, Float16, 0}[1.0, 2.0, 3.0])
    setter = setsym_oop(prob, idxs)
    @test_nowarn @inferred setter(prob, vals)
end

ps = MTKParameters(sys, ivs)
function loss(value, sys, ps)
    @test value isa ForwardDiff.Dual
    ps = remake_buffer(sys, ps, (a,), (value,))
    getp(sys, a)(ps) + getp(sys, b)(ps)
end

@test ForwardDiff.derivative(x -> loss(x, sys, ps), 1.5) == 3.0

# Issue#2615
@parameters p::Vector{Float64}
@variables X(t)
eq = D(X) ~ p[1] - p[2] * X
@mtkcompile osys = System([eq], t)

u0 = [X => 1.0]
ps = [p => [2.0, 0.1]]
p = MTKParameters(osys, [ps; u0])
@test p.tunable == [2.0, 0.1]

# Ensure partial update promotes the buffer
@parameters p q r
@named sys = System(Equation[], t, [], [p, q, r])
sys = complete(sys)
ps = MTKParameters(sys, [p => 1.0, q => 2.0, r => 3.0])
newps = remake_buffer(sys, ps, (p,), (1.0f0,))
@test newps.tunable isa Vector{Float32}
@test newps.tunable == [1.0f0, 2.0f0, 3.0f0]

# Issue#2624
@parameters p d
@variables X(t)
eqs = [D(X) ~ p - d * X]
@mtkcompile sys = System(eqs, t)

u0 = [X => 1.0]
tspan = (0.0, 100.0)
ps = [p => 1.0] # Value for `d` is missing

@test_throws ModelingToolkit.MissingParametersError ODEProblem(sys, [u0; ps], tspan)
@test_nowarn ODEProblem(sys, [u0; ps; [d => 1.0]], tspan)

# JET tests

# scalar parameters only
function level1()
    @parameters p1=0.5 [tunable=true] p2=1 [tunable=true] p3=3 [tunable=false] p4=3 [tunable=true] y0=1
    @variables x(t)=2 y(t)=y0
    D = Differential(t)

    eqs = [D(x) ~ p1 * x - p2 * x * y
           D(y) ~ -p3 * y + p4 * x * y
           y0 ~ 2p4]

    sys = mtkcompile(complete(System(
        eqs, t, name = :sys)))
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(sys, [], (0.0, 3.0))
end

# scalar and vector parameters
function level2()
    @parameters p1=0.5 [tunable=true] (p23[1:2]=[1, 3.0]) [tunable=true] p4=3 [tunable=false] y0=1
    @variables x(t)=2 y(t)=y0
    D = Differential(t)

    eqs = [D(x) ~ p1 * x - p23[1] * x * y
           D(y) ~ -p23[2] * y + p4 * x * y
           y0 ~ 2p4]

    sys = mtkcompile(complete(System(
        eqs, t, name = :sys)))
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(sys, [], (0.0, 3.0))
end

# scalar and vector parameters with different scalar types
function level3()
    @parameters p1=0.5 [tunable=true] (p23[1:2]=[1, 3.0]) [tunable=true] p4::Int=3 [tunable=true] y0::Int=1
    @variables x(t)=2 y(t)=y0
    D = Differential(t)

    eqs = [D(x) ~ p1 * x - p23[1] * x * y
           D(y) ~ -p23[2] * y + p4 * x * y
           y0 ~ 2p4]

    sys = mtkcompile(complete(System(
        eqs, t, name = :sys)))
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(sys, [], (0.0, 3.0))
end

@testset "level$i" for (i, prob) in enumerate([level1(), level2(), level3()])
    ps = prob.p
    @testset "Type stability of $portion" for portion in [
        Tunable(), Discrete(), Constants()]
        @test_call canonicalize(portion, ps)
        @inferred canonicalize(portion, ps)

        # broken because the size of a vector of vectors can't be determined at compile time
        @test_opt target_modules=(ModelingToolkit,) canonicalize(
            portion, ps)

        buffer, repack, alias = canonicalize(portion, ps)

        # broken because dependent update functions break inference
        @test_call target_modules=(ModelingToolkit,) SciMLStructures.replace(
            portion, ps, ones(length(buffer)))
        @inferred SciMLStructures.replace(
            portion, ps, ones(length(buffer)))
        @inferred MTKParameters SciMLStructures.replace(portion, ps, ones(length(buffer)))
        @test_opt target_modules=(ModelingToolkit,) SciMLStructures.replace(
            portion, ps, ones(length(buffer)))

        @test_call target_modules=(ModelingToolkit,) SciMLStructures.replace!(
            portion, ps, ones(length(buffer)))
        @inferred SciMLStructures.replace!(portion, ps, ones(length(buffer)))
        @test_opt target_modules=(ModelingToolkit,) SciMLStructures.replace!(
            portion, ps, ones(length(buffer)))
    end
end

# Issue#2642
@parameters α β γ δ
@variables x(t) y(t)
eqs = [D(x) ~ (α - β * y) * x
       D(y) ~ (δ * x - γ) * y]
@mtkcompile odesys = System(eqs, t)
odeprob = ODEProblem(
    odesys, [x => 1.0, y => 1.0, α => 1.5, β => 1.0, γ => 3.0, δ => 1.0], (0.0, 10.0))
tunables, _... = canonicalize(Tunable(), odeprob.p)
@test tunables isa AbstractVector{Float64}

function loss(x)
    ps = odeprob.p
    newps = SciMLStructures.replace(Tunable(), ps, x)
    newprob = remake(odeprob, p = newps)
    sol = solve(newprob, Tsit5())
    return sum(sol)
end

@test_nowarn ForwardDiff.gradient(loss, collect(tunables))

VDual = Vector{<:ForwardDiff.Dual}
VVDual = Vector{<:Vector{<:ForwardDiff.Dual}}

@testset "Parameter type validation" begin
    struct Foo{T}
        x::T
    end

    @parameters a b::Int c::Vector{Float64} d[1:2, 1:2]::Int e::Foo{Int} f::Foo
    @named sys = System(Equation[], t, [], [a, b, c, d, e, f])
    sys = complete(sys)
    ps = MTKParameters(sys,
        Dict(a => 1.0, b => 2, c => 3ones(2),
            d => 3ones(Int, 2, 2), e => Foo(1), f => Foo("a")))
    @test_nowarn setp(sys, c)(ps, ones(4)) # so this is fixed when SII is fixed
    @test_throws DimensionMismatch set_parameter!(
        ps, 4ones(Int, 3, 2), parameter_index(sys, d))
    @test_throws DimensionMismatch set_parameter!(
        ps, 4ones(Int, 4), parameter_index(sys, d)) # size has to match, not just length
    @test_nowarn setp(sys, f)(ps, Foo(:a)) # can change non-concrete type

    # Same flexibility is afforded to `b::Int` to allow for ForwardDiff
    for sym in [a, b]
        @test_nowarn remake_buffer(sys, ps, (sym,), (1,))
        newps = @test_nowarn remake_buffer(sys, ps, (sym,), (1.0f0,)) # Can change type if it's numeric
        @test getp(sys, sym)(newps) isa Float32
        newps = @test_nowarn remake_buffer(sys, ps, sym, ForwardDiff.Dual(1.0))
        @test getp(sys, sym)(newps) isa ForwardDiff.Dual
        @test_throws TypeError remake_buffer(sys, ps, (sym,), (:a,)) # still has to be numeric
    end

    newps = @test_nowarn remake_buffer(sys, ps, (c,), (view(1.0:4.0, 2:4),)) # can change type of array
    @test getp(sys, c)(newps) == 2.0:4.0
    @test parameter_values(newps, parameter_index(sys, c)) ≈ [2.0, 3.0, 4.0]
    @test_throws TypeError remake_buffer(sys, ps, (c,), ([:a, :b, :c],)) # can't arbitrarily change eltype
    @test_throws TypeError remake_buffer(sys, ps, (c,), (:a,)) # can't arbitrarily change type

    newps = @test_nowarn remake_buffer(sys, ps, (d,), (ForwardDiff.Dual.(ones(2, 2)),)) # can change eltype
    @test_throws TypeError remake_buffer(sys, ps, (d,), ([:a :b; :c :d],)) # eltype still has to be numeric
    @test getp(sys, d)(newps) isa Matrix{<:ForwardDiff.Dual}

    @test_throws TypeError remake_buffer(sys, ps, (e,), (Foo(2.0),)) # need exact same type for nonnumeric
    @test_nowarn remake_buffer(sys, ps, (f,), (Foo(:a),))
end

@testset "Error on missing parameter defaults" begin
    @parameters a b c
    @named sys = System(Equation[], t, [], [a, b]; defaults = Dict(b => 2c))
    sys = complete(sys)
    @test_throws ["Could not evaluate", "b", "Missing", "2c"] MTKParameters(sys, [a => 1.0])
end

@testset "Issue#2804" begin
    @parameters k[1:4]
    @variables (V(t))[1:2]
    eqs = [
        D(V[1]) ~ k[1] - k[2] * V[1],
        D(V[2]) ~ k[3] - k[4] * V[2]
    ]
    @mtkcompile osys_scal = System(eqs, t, [V[1], V[2]], [k[1], k[2], k[3], k[4]])

    u0 = [V => [10.0, 20.0]]
    ps_vec = [k => [2.0, 3.0, 4.0, 5.0]]
    ps_scal = [k[1] => 1.0, k[2] => 2.0, k[3] => 3.0, k[4] => 4.0]
    oprob_scal_scal = ODEProblem(osys_scal, [u0; ps_scal], 1.0)
    newoprob = remake(oprob_scal_scal; p = ps_vec, build_initializeprob = false)
    @test newoprob.ps[k] == [2.0, 3.0, 4.0, 5.0]
end

# Parameter timeseries
ps = MTKParameters([1.0, 1.0], (), (BlockedArray(zeros(4), [2, 2]),),
    (), (), ())
ps2 = SciMLStructures.replace(Discrete(), ps, ones(4))
@test typeof(ps2.discrete) == typeof(ps.discrete)
with_updated_parameter_timeseries_values(
    sys, ps, 1 => ModelingToolkit.NestedGetIndex(([5.0, 10.0],)))
@test ps.discrete[1][Block(1)] == [5.0, 10.0]
with_updated_parameter_timeseries_values(
    sys, ps, 1 => ModelingToolkit.NestedGetIndex(([3.0, 30.0],)),
    2 => ModelingToolkit.NestedGetIndex(([4.0, 40.0],)))
@test ps.discrete[1][Block(1)] == [3.0, 30.0]
@test ps.discrete[1][Block(2)] == [4.0, 40.0]
@test SciMLBase.get_saveable_values(sys, ps, 1).x == (ps.discrete[1][Block(1)],)

# With multiple types and clocks
ps = MTKParameters(
    (), (),
    (BlockedArray([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [3, 3]),
        BlockedArray(falses(1), [1, 0])),
    (), (), ())
@test SciMLBase.get_saveable_values(sys, ps, 1).x isa Tuple{Vector{Float64}, BitVector}
tsidx1 = 1
tsidx2 = 2
@test length(ps.discrete[1][Block(tsidx1)]) == 3
@test length(ps.discrete[2][Block(tsidx1)]) == 1
@test length(ps.discrete[1][Block(tsidx2)]) == 3
@test length(ps.discrete[2][Block(tsidx2)]) == 0
with_updated_parameter_timeseries_values(
    sys, ps, tsidx1 => ModelingToolkit.NestedGetIndex(([10.0, 11.0, 12.0], [false])))
@test ps.discrete[1][Block(tsidx1)] == [10.0, 11.0, 12.0]
@test ps.discrete[2][Block(tsidx1)][] == false

@testset "Avoid specialization of nonnumeric parameters on `remake_buffer`" begin
    @variables x(t)
    @parameters p::Any
    @named sys = System(D(x) ~ x, t, [x], [p])
    sys = complete(sys)
    ps = MTKParameters(sys, [p => 1.0])
    @test ps.nonnumeric isa Tuple{Vector{Any}}
    ps2 = remake_buffer(sys, ps, [p], [:a])
    @test ps2.nonnumeric isa Tuple{Vector{Any}}
end

@testset "Issue#3925: Autodiff after `subset_tunables`" begin
    function circuit_model()
        @named resistor1 = Resistor(R=5.0)
        @named resistor2 = Resistor(R=2.0)
        @named capacitor1 = Capacitor(C=2.4)
        @named capacitor2 = Capacitor(C=60.0)
        @named source = Voltage()
        @named input_signal = Sine(frequency=1.0)
        @named ground = Ground()
        @named ampermeter = CurrentSensor()

        eqs = [connect(input_signal.output, source.V)
            connect(source.p, capacitor1.n, capacitor2.n)
            connect(source.n, resistor1.p, resistor2.p, ground.g)
            connect(resistor1.n, capacitor1.p, ampermeter.n)
            connect(resistor2.n, capacitor2.p, ampermeter.p)]

        @named circuit_model = System(eqs, t,
            systems=[
                resistor1, resistor2, capacitor1, capacitor2,
                source, input_signal, ground, ampermeter
            ])
    end

    model = circuit_model()
    sys = mtkcompile(model)

    tunable_parameters(sys)

    sub_sys = subset_tunables(sys, [sys.capacitor2.C])

    tunable_parameters(sub_sys)

    prob = ODEProblem(sub_sys, [sys.capacitor2.v => 0.0], (0, 3.))

    setter = setsym_oop(prob, [sys.capacitor2.C]);

    function loss(x, ps)
        setter, prob = ps
        u0, p = setter(prob, x)
        new_prob = remake(prob; u0, p)
        sol = solve(new_prob, Rodas5P())
        sum(sol)
    end

    grad = ForwardDiff.gradient(Base.Fix2(loss, (setter, prob)), [3.0])
    @test grad ≈ [0.14882627068752538] atol=1e-10
end
