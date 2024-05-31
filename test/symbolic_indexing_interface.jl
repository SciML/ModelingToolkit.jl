using ModelingToolkit, SymbolicIndexingInterface, SciMLBase
using ModelingToolkit: t_nounits as t, D_nounits as D, ParameterIndex
using SciMLStructures: Tunable

@testset "ODESystem" begin
    @parameters a b
    @variables x(t)=1.0 y(t)=2.0 xy(t)
    eqs = [D(x) ~ a * y + t, D(y) ~ b * t]
    @named odesys = ODESystem(eqs, t, [x, y], [a, b]; observed = [xy ~ x + y])
    odesys = complete(odesys)
    @test all(is_variable.((odesys,), [x, y, 1, 2, :x, :y]))
    @test all(.!is_variable.((odesys,), [a, b, t, 3, 0, :a, :b]))
    @test variable_index.((odesys,), [x, y, a, b, t, 1, 2, :x, :y, :a, :b]) ==
          [1, 2, nothing, nothing, nothing, 1, 2, 1, 2, nothing, nothing]
    @test isequal(variable_symbols(odesys), [x, y])
    @test all(is_parameter.((odesys,), [a, b, ParameterIndex(Tunable(), (1, 1)), :a, :b]))
    @test all(.!is_parameter.((odesys,), [x, y, t, 3, 0, :x, :y]))
    @test parameter_index(odesys, a) == parameter_index(odesys, :a)
    @test parameter_index(odesys, a) isa ParameterIndex{Tunable, Tuple{Int, Int}}
    @test parameter_index(odesys, b) == parameter_index(odesys, :b)
    @test parameter_index(odesys, b) isa ParameterIndex{Tunable, Tuple{Int, Int}}
    @test parameter_index.(
        (odesys,), [x, y, t, ParameterIndex(Tunable(), (1, 1)), :x, :y]) ==
          [nothing, nothing, nothing, ParameterIndex(Tunable(), (1, 1)), nothing, nothing]
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

    @named odesys = ODESystem(
        eqs, t, [x, y], [a, b]; defaults = [xy => 3.0], observed = [xy ~ x + y])
    odesys = complete(odesys)
    @test default_values(odesys)[xy] == 3.0
    pobs = parameter_observed(odesys, a + b)
    @test pobs.timeseries_idx === nothing
    @test pobs.observed_fn(
        ModelingToolkit.MTKParameters(odesys, [a => 1.0, b => 2.0]), 0.0) ≈ 3.0
    pobs = parameter_observed(odesys, [a + b, a - b])
    @test pobs.timeseries_idx === nothing
    @test pobs.observed_fn(
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
    @test !is_time_dependent(ns)
    ps = ModelingToolkit.MTKParameters(ns, [σ => 1.0, ρ => 2.0, β => 3.0])
    pobs = parameter_observed(ns, σ + ρ)
    @test pobs.timeseries_idx === nothing
    @test pobs.observed_fn(ps) == 3.0
    pobs = parameter_observed(ns, [σ + ρ, ρ + β])
    @test pobs.timeseries_idx === nothing
    @test pobs.observed_fn(ps) == [3.0, 5.0]
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
    analytic = [u(t, x) ~ -h * x * (x - 1) * sin(x) * exp(-2 * h * t)]
    analytic_function = (ps, t, x) -> -ps[1] * x * (x - 1) * sin(x) * exp(-2 * ps[1] * t)

    @named pdesys = PDESystem(eq, bcs, domains, [t, x], [u], [h], analytic = analytic)

    @test isequal(pdesys.ps, [h])
    @test isequal(parameter_symbols(pdesys), [h])
    @test isequal(parameters(pdesys), [h])
end
