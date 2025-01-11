using ModelingToolkit, SymbolicIndexingInterface, SciMLBase
using ModelingToolkit: t_nounits as t, D_nounits as D, ParameterIndex
using SciMLStructures: Tunable

@testset "ODESystem" begin
    @parameters a b
    @variables x(t)=1.0 y(t)=2.0 xy(t)
    eqs = [D(x) ~ a * y + t, D(y) ~ b * t]
    @named odesys = ODESystem(eqs, t, [x, y], [a, b]; observed = [xy ~ x + y])
    odesys = complete(odesys)
    @test SymbolicIndexingInterface.supports_tuple_observed(odesys)
    @test all(is_variable.((odesys,), [x, y, 1, 2, :x, :y]))
    @test all(.!is_variable.((odesys,), [a, b, t, 3, 0, :a, :b]))
    @test variable_index.((odesys,), [x, y, a, b, t, 1, 2, :x, :y, :a, :b]) ==
          [1, 2, nothing, nothing, nothing, 1, 2, 1, 2, nothing, nothing]
    @test isequal(variable_symbols(odesys), [x, y])
    @test all(is_parameter.((odesys,), [a, b, ParameterIndex(Tunable(), 1), :a, :b]))
    @test all(.!is_parameter.((odesys,), [x, y, t, 3, 0, :x, :y]))
    @test parameter_index(odesys, a) == parameter_index(odesys, :a)
    @test parameter_index(odesys, a) isa ParameterIndex{Tunable, Int}
    @test parameter_index(odesys, b) == parameter_index(odesys, :b)
    @test parameter_index(odesys, b) isa ParameterIndex{Tunable, Int}
    @test parameter_index.(
        (odesys,), [x, y, t, ParameterIndex(Tunable(), 1), :x, :y]) ==
          [nothing, nothing, nothing, ParameterIndex(Tunable(), 1), nothing, nothing]
    @test isequal(parameter_symbols(odesys), [a, b])
    @test all(is_independent_variable.((odesys,), [t, :t]))
    @test all(.!is_independent_variable.((odesys,), [x, y, a, :x, :y, :a]))
    @test isequal(independent_variable_symbols(odesys), [t])
    @test is_time_dependent(odesys)
    @test constant_structure(odesys)
    @test !isempty(default_values(odesys))
    @test default_values(odesys)[x] == 1.0
    @test default_values(odesys)[y] == 2.0
    @test isequal(default_values(odesys)[xy], x + y)

    prob = ODEProblem(odesys, [], (0.0, 1.0), [a => 1.0, b => 2.0])
    getter = getu(odesys, (x + 1, x + 2))
    @test getter(prob) isa Tuple
    @test_nowarn @inferred getter(prob)
    getter = getp(odesys, (a + 1, a + 2))
    @test getter(prob) isa Tuple
    @test_nowarn @inferred getter(prob)

    @named odesys = ODESystem(
        eqs, t, [x, y], [a, b]; defaults = [xy => 3.0], observed = [xy ~ x + y])
    odesys = complete(odesys)
    @test default_values(odesys)[xy] == 3.0
    pobs = parameter_observed(odesys, a + b)
    @test isempty(get_all_timeseries_indexes(odesys, a + b))
    @test pobs(
        ModelingToolkit.MTKParameters(odesys, [a => 1.0, b => 2.0]), 0.0) ≈ 3.0
    pobs = parameter_observed(odesys, [a + b, a - b])
    @test isempty(get_all_timeseries_indexes(odesys, [a + b, a - b]))
    @test pobs(
        ModelingToolkit.MTKParameters(odesys, [a => 1.0, b => 2.0]), 0.0) ≈ [3.0, -1.0]
end

# @testset "Clock system" begin
#     dt = 0.1
#     dt2 = 0.2
#     @variables x(t)=0 y(t)=0 u(t)=0 yd1(t)=0 ud1(t)=0 yd2(t)=0 ud2(t)=0
#     @parameters kp=1 r=1

#     eqs = [
#            # controller (time discrete part `dt=0.1`)
#            yd1 ~ Sample(t, dt)(y)
#            ud1 ~ kp * (r - yd1)
#            # controller (time discrete part `dt=0.2`)
#            yd2 ~ Sample(t, dt2)(y)
#            ud2 ~ kp * (r - yd2)

#            # plant (time continuous part)
#            u ~ Hold(ud1) + Hold(ud2)
#            D(x) ~ -x + u
#            y ~ x]

#     @mtkbuild cl = ODESystem(eqs, t)
#     partition1_params = [Hold(ud1), Sample(t, dt)(y), ud1, yd1]
#     partition2_params = [Hold(ud2), Sample(t, dt2)(y), ud2, yd2]
#     @test all(
#         Base.Fix1(is_timeseries_parameter, cl), vcat(partition1_params, partition2_params))
#     @test allequal(timeseries_parameter_index(cl, p).timeseries_idx
#     for p in partition1_params)
#     @test allequal(timeseries_parameter_index(cl, p).timeseries_idx
#     for p in partition2_params)
#     tsidx1 = timeseries_parameter_index(cl, partition1_params[1]).timeseries_idx
#     tsidx2 = timeseries_parameter_index(cl, partition2_params[1]).timeseries_idx
#     @test tsidx1 != tsidx2
#     ps = ModelingToolkit.MTKParameters(cl, [kp => 1.0, Sample(t, dt)(y) => 1.0])
#     pobs = parameter_observed(cl, Shift(t, 1)(yd1))
#     @test pobs.timeseries_idx == tsidx1
#     @test pobs.observed_fn(ps, 0.0) == 1.0
#     pobs = parameter_observed(cl, [Shift(t, 1)(yd1), Shift(t, 1)(ud1)])
#     @test pobs.timeseries_idx == tsidx1
#     @test pobs.observed_fn(ps, 0.0) == [1.0, 0.0]
#     pobs = parameter_observed(cl, [Shift(t, 1)(yd1), Shift(t, 1)(ud2)])
#     @test pobs.timeseries_idx === nothing
#     @test pobs.observed_fn(ps, 0.0) == [1.0, 1.0]
# end

@testset "Nonlinear system" begin
    @variables x y z
    @parameters σ ρ β

    eqs = [0 ~ σ * (y - x),
        0 ~ x * (ρ - z) - y,
        0 ~ x * y - β * z]
    @named ns = NonlinearSystem(eqs, [x, y, z], [σ, ρ, β])
    ns = complete(ns)
    @test SymbolicIndexingInterface.supports_tuple_observed(ns)
    @test !is_time_dependent(ns)
    ps = ModelingToolkit.MTKParameters(ns, [σ => 1.0, ρ => 2.0, β => 3.0])
    pobs = parameter_observed(ns, σ + ρ)
    @test isempty(get_all_timeseries_indexes(ns, σ + ρ))
    @test pobs(ps) == 3.0
    pobs = parameter_observed(ns, [σ + ρ, ρ + β])
    @test isempty(get_all_timeseries_indexes(ns, [σ + ρ, ρ + β]))
    @test pobs(ps) == [3.0, 5.0]

    prob = NonlinearProblem(
        ns, [x => 1.0, y => 2.0, z => 3.0], [σ => 1.0, ρ => 2.0, β => 3.0])
    getter = getu(ns, (x + 1, x + 2))
    @test getter(prob) isa Tuple
    @test_nowarn @inferred getter(prob)
    getter = getp(ns, (σ + 1, σ + 2))
    @test getter(prob) isa Tuple
    @test_nowarn @inferred getter(prob)
end

@testset "PDESystem" begin
    @parameters x
    @variables u(..)
    Dxx = Differential(x)^2
    Dtt = Differential(t)^2
    Dt = D

    #2D PDE
    C = 1
    eq = Dtt(u(t, x)) ~ C^2 * Dxx(u(t, x))

    # Initial and boundary conditions
    bcs = [u(t, 0) ~ 0.0,# for all t > 0
        u(t, 1) ~ 0.0,# for all t > 0
        u(0, x) ~ x * (1.0 - x), #for all 0 < x < 1
        Dt(u(0, x)) ~ 0.0] #for all  0 < x < 1]

    # Space and time domains
    domains = [t ∈ (0.0, 1.0),
        x ∈ (0.0, 1.0)]

    @named pde_system = PDESystem(eq, bcs, domains, [t, x], [u])

    @test pde_system.ps == SciMLBase.NullParameters()
    @test parameter_symbols(pde_system) == []

    @parameters x
    @constants h = 1
    @variables u(..)
    Dt = D
    Dxx = Differential(x)^2
    eq = Dt(u(t, x)) ~ h * Dxx(u(t, x))
    bcs = [u(0, x) ~ -h * x * (x - 1) * sin(x),
        u(t, 0) ~ 0, u(t, 1) ~ 0]

    domains = [t ∈ (0.0, 1.0),
        x ∈ (0.0, 1.0)]

    analytic = [u(t, x) ~ -h * x * (x - 1) * sin(x) * exp(-2 * h * t)]
    analytic_function = (ps, t, x) -> -ps[1] * x * (x - 1) * sin(x) * exp(-2 * ps[1] * t)

    @named pdesys = PDESystem(eq, bcs, domains, [t, x], [u], [h], analytic = analytic)

    @test isequal(pdesys.ps, [h])
    @test isequal(parameter_symbols(pdesys), [h])
    @test isequal(parameters(pdesys), [h])
end

# Issue#2767
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SymbolicIndexingInterface

@parameters p1[1:2]=[1.0, 2.0] p2[1:2]=[0.0, 0.0]
@variables x(t) = 0

@named sys = ODESystem(
    [D(x) ~ sum(p1) * t + sum(p2)],
    t;
)
prob = ODEProblem(complete(sys))
get_dep = @test_nowarn getu(prob, 2p1)
@test get_dep(prob) == [2.0, 4.0]

@testset "Observed functions with variables as `Symbol`s" begin
    @variables x(t) y(t) z(t)[1:2]
    @parameters p1 p2[1:2, 1:2]
    @mtkbuild sys = ODESystem([D(x) ~ x * t + p1, y ~ 2x, D(z) ~ p2 * z], t)
    prob = ODEProblem(
        sys, [x => 1.0, z => ones(2)], (0.0, 1.0), [p1 => 2.0, p2 => ones(2, 2)])
    @test getu(prob, x)(prob) == getu(prob, :x)(prob)
    @test getu(prob, [x, y])(prob) == getu(prob, [:x, :y])(prob)
    @test getu(prob, z)(prob) == getu(prob, :z)(prob)
    @test getu(prob, p1)(prob) == getu(prob, :p1)(prob)
    @test getu(prob, p2)(prob) == getu(prob, :p2)(prob)
end

@testset "Parameter dependencies as symbols" begin
    @variables x(t) = 1.0
    @parameters a=1 b
    @named model = ODESystem(D(x) ~ x + a - b, t, parameter_dependencies = [b ~ a + 1])
    sys = complete(model)
    prob = ODEProblem(sys, [], (0.0, 1.0))
    @test prob.ps[b] == prob.ps[:b]
end

@testset "`get_all_timeseries_indexes` with non-split systems" begin
    @variables x(t) y(t) z(t)
    @parameters a
    @named sys = ODESystem([D(x) ~ a * x, y ~ 2x, z ~ 0.0], t)
    sys = structural_simplify(sys, split = false)
    for sym in [x, y, z, x + y, x + a, y / x]
        @test only(get_all_timeseries_indexes(sys, sym)) == ContinuousTimeseries()
    end
    @test isempty(get_all_timeseries_indexes(sys, a))
end

@testset "`timeseries_parameter_index` on unwrapped scalarized timeseries parameter" begin
    @variables x(t)[1:2]
    @parameters p(t)[1:2, 1:2]
    ev = [x[1] ~ 2.0] => [p ~ -ones(2, 2)]
    @mtkbuild sys = ODESystem(D(x) ~ p * x, t; continuous_events = [ev])
    p = ModelingToolkit.unwrap(p)
    @test timeseries_parameter_index(sys, p) === ParameterTimeseriesIndex(1, (1, 1))
    @test timeseries_parameter_index(sys, p[1, 1]) ===
          ParameterTimeseriesIndex(1, (1, 1, 1, 1))
end
