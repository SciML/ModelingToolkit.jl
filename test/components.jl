using Test
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: get_component_type
using ModelingToolkit.BipartiteGraphs
using ModelingToolkit.StructuralTransformations
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks
using LinearAlgebra
using ModelingToolkitStandardLibrary.Thermal
using SymbolicUtils: getmetadata
include("common/rc_model.jl")

@testset "Basics" begin
    @unpack resistor, capacitor, source = rc_model
    function check_contract(sys)
        state = ModelingToolkit.get_tearing_state(sys)
        graph = state.structure.graph
        fullvars = state.fullvars
        sys = tearing_substitution(sys)

        eqs = equations(sys)
        for (i, eq) in enumerate(eqs)
            actual = union(ModelingToolkit.vars(eq.lhs), ModelingToolkit.vars(eq.rhs))
            actual = filter(!ModelingToolkit.isparameter, collect(actual))
            current = Set(fullvars[𝑠neighbors(graph, i)])
            @test isempty(setdiff(actual, current))
        end
    end

    function check_rc_sol(sol)
        rpi = sol[rc_model.resistor.p.i]
        rpifun = sol.prob.f.observed(rc_model.resistor.p.i)
        @test rpifun.(sol.u, (sol.prob.p,), sol.t) == rpi
        @test any(!isequal(rpi[1]), rpi) # test that we don't have a constant system
        @test sol[rc_model.resistor.p.i] == sol[resistor.p.i] == sol[capacitor.p.i]
        @test sol[rc_model.resistor.n.i] == sol[resistor.n.i] == -sol[capacitor.p.i]
        @test sol[rc_model.capacitor.n.i] == sol[capacitor.n.i] == -sol[capacitor.p.i]
        @test iszero(sol[rc_model.ground.g.i])
        @test iszero(sol[rc_model.ground.g.v])
        @test sol[rc_model.resistor.v] == sol[resistor.v] ==
              sol[source.p.v] - sol[capacitor.p.v]
    end

    @named pin = Pin()
    @test get_component_type(pin).name == :Pin
    @test get_component_type(rc_model.resistor).name == :Resistor

    @mtkcompile rc_model_compile = RCModel()
    @test get_component_type(rc_model).name == :RCModel
    @test get_component_type(rc_model_compile).name == :RCModel

    completed_rc_model = complete(rc_model)
    @test isequal(completed_rc_model.resistor.n.i, resistor.n.i)
    @test ModelingToolkit.n_expanded_connection_equations(capacitor) == 2
    @test length(equations(mtkcompile(rc_model, allow_parameter = false))) == 2
    sys = mtkcompile(rc_model)
    @test_throws ModelingToolkit.RepeatedStructuralSimplificationError mtkcompile(sys)
    @test length(equations(sys)) == 1
    check_contract(sys)
    @test !isempty(ModelingToolkit.defaults(sys))
    u0 = [capacitor.v => 0.0]
    prob = ODEProblem(sys, u0, (0, 10.0))
    sol = solve(prob, Rodas4())
    check_rc_sol(sol)
end

@testset "Outer/inner connections" begin
    sys = mtkcompile(rc_model)

    prob = ODEProblem(sys, [sys.capacitor.v => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4())
    function rc_component(; name, R = 1, C = 1)
        local sys
        @parameters R=R C=C
        @named p = Pin()
        @named n = Pin()
        @named resistor = Resistor(R = R) # test parent scope default of @named
        @named capacitor = Capacitor(C = C)
        eqs = [connect(p, resistor.p);
               connect(resistor.n, capacitor.p);
               connect(capacitor.n, n)]
        @named sys = System(eqs, t, [], [R, C])
        compose(sys, [p, n, resistor, capacitor]; name = name)
    end

    @named ground = Ground()
    @named shape = Constant(k = 1)
    @named source = Voltage()
    @named rc_comp = rc_component()
    eqs = [connect(shape.output, source.V)
           connect(source.p, rc_comp.p)
           connect(source.n, rc_comp.n)
           connect(source.n, ground.g)]
    @named sys′ = System(eqs, t)
    @named sys_inner_outer = compose(sys′, [ground, shape, source, rc_comp])
    @test_nowarn show(IOBuffer(), MIME"text/plain"(), sys_inner_outer)
    expand_connections(sys_inner_outer)
    sys_inner_outer = mtkcompile(sys_inner_outer)
    @test !isempty(ModelingToolkit.defaults(sys_inner_outer))
    u0 = [rc_comp.capacitor.v => 0.0]
    prob = ODEProblem(sys_inner_outer, u0, (0, 10.0), sparse = true)
    sol_inner_outer = solve(prob, Rodas4())
    @test sol[sys.capacitor.v] ≈ sol_inner_outer[rc_comp.capacitor.v]

    prob = ODEProblem(sys, [sys.capacitor.v => 0.0], (0, 10.0))
    sol = solve(prob, Tsit5())

    @test sol[sys.resistor.p.i] == sol[sys.capacitor.p.i]
    @test sol[sys.resistor.n.i] == -sol[sys.capacitor.p.i]
    @test sol[sys.capacitor.n.i] == -sol[sys.capacitor.p.i]
    @test iszero(sol[sys.ground.g.i])
    @test iszero(sol[sys.ground.g.v])
    @test sol[sys.resistor.v] == sol[sys.source.p.v] - sol[sys.capacitor.p.v]
end
#using Plots
#plot(sol)

include("common/serial_inductor.jl")
@testset "Serial inductor" begin
    sys = mtkcompile(ll_model)
    @test length(equations(sys)) == 2
    u0 = unknowns(sys) .=> 0
    @test_nowarn ODEProblem(
        sys, [], (0, 10.0), guesses = u0, warn_initialize_determined = false)
    prob = DAEProblem(sys, D.(unknowns(sys)) .=> 0, (0, 0.5), guesses = u0)
    sol = solve(prob, DFBDF())
    @test sol.retcode == SciMLBase.ReturnCode.Success

    sys2 = mtkcompile(ll2_model)
    @test length(equations(sys2)) == 3
    u0 = [sys.inductor2.i => 0]
    prob = ODEProblem(sys, u0, (0, 10.0))
    sol = solve(prob, FBDF())
    @test SciMLBase.successful_retcode(sol)
end

@testset "Compose/extend" begin
    @variables x1(t) x2(t) x3(t) x4(t)
    @named sys1_inner = System([D(x1) ~ x1], t)
    @named sys1_partial = compose(System([D(x2) ~ x2], t; name = :foo), sys1_inner)
    @named sys1 = extend(System([D(x3) ~ x3], t; name = :foo), sys1_partial)
    @named sys2 = compose(System([D(x4) ~ x4], t; name = :foo), sys1)
    @test_nowarn sys2.sys1.sys1_inner.x1 # test the correct nesting

    # compose tests
    function record_fun(; name)
        pars = @parameters a=10 b=100
        System(Equation[], t, [], pars; name)
    end

    function first_model(; name)
        @named foo = record_fun()

        defs = Dict()
        defs[foo.a] = 3
        defs[foo.b] = 300
        pars = @parameters x=2 y=20
        compose(System(Equation[], t, [], pars; name, defaults = defs), foo)
    end
    @named goo = first_model()
    @unpack foo = goo
    @test ModelingToolkit.defaults(goo)[foo.a] == 3
    @test ModelingToolkit.defaults(goo)[foo.b] == 300
end

function Load(; name)
    R = 1
    @named p = Pin()
    @named n = Pin()
    @named resistor = Resistor(R = R)
    eqs = [connect(p, resistor.p);
           connect(resistor.n, n)]
    @named sys = System(eqs, t)
    compose(sys, [p, n, resistor]; name = name)
end

function Circuit(; name)
    R = 1
    @named ground = Ground()
    @named load = Load()
    @named resistor = Resistor(R = R)
    eqs = [connect(load.p, ground.g);
           connect(resistor.p, ground.g)]
    @named sys = System(eqs, t)
    compose(sys, [ground, resistor, load]; name = name)
end

@named foo = Circuit()
@test mtkcompile(foo) isa ModelingToolkit.AbstractSystem

# BLT tests
@testset "BLT ordering" begin
    function parallel_rc_model(i; name, shape, source, ground, R, C)
        resistor = Resistor(name = Symbol(:resistor, i), R = R, T_dep = true)
        capacitor = Capacitor(name = Symbol(:capacitor, i), C = C)
        heat_capacitor = HeatCapacitor(name = Symbol(:heat_capacitor, i))

        rc_eqs = [connect(shape.output, source.V)
                  connect(source.p, resistor.p)
                  connect(resistor.n, capacitor.p)
                  connect(capacitor.n, source.n, ground.g)
                  connect(resistor.heat_port, heat_capacitor.port)]

        compose(System(rc_eqs, t, name = Symbol(name, i)),
            [resistor, capacitor, source, ground, shape, heat_capacitor])
    end
    V = 2.0
    @named shape = Constant(k = V)
    @named source = Voltage()
    @named ground = Ground()
    N = 50
    Rs = 10 .^ range(0, stop = -4, length = N)
    Cs = 10 .^ range(-3, stop = 0, length = N)
    rc_systems = map(1:N) do i
        parallel_rc_model(i; name = :rc, source, ground, shape, R = Rs[i], C = Cs[i])
    end
    @variables E(t) = 0.0
    eqs = [
        D(E) ~ sum(((i, sys),) -> getproperty(sys, Symbol(:resistor, i)).heat_port.Q_flow,
        enumerate(rc_systems))
    ]
    @named _big_rc = System(eqs, t, [E], [])
    @named big_rc = compose(_big_rc, rc_systems)
    ts = TearingState(expand_connections(big_rc))
    # this is block upper triangular, so `istriu` needs a little leeway
    @test istriu(but_ordered_incidence(ts)[1], -2)
end

# Test using constants inside subsystems
@testset "Constants inside subsystems" begin
    function FixedResistor(; name, R = 1.0)
        @named oneport = OnePort()
        @unpack v, i = oneport
        @constants R = R
        eqs = [
            v ~ i * R
        ]
        extend(System(eqs, t, [], [R]; name = name), oneport)
    end
    capacitor = Capacitor(; name = :c1, C = 1.0)
    resistor = FixedResistor(; name = :r1)
    ground = Ground(; name = :ground)
    rc_eqs = [connect(capacitor.n, resistor.p)
              connect(resistor.n, capacitor.p)
              connect(capacitor.n, ground.g)]

    @named _rc_model = System(rc_eqs, t)
    @named rc_model = compose(_rc_model,
        [resistor, capacitor, ground])
    sys = mtkcompile(rc_model)
    prob = ODEProblem(sys, [sys.c1.v => 0.0], (0, 10.0))
    sol = solve(prob, Tsit5())
end

@testset "docstrings (#1155)" begin
    """
    Hey there, Pin1!
    """
    @connector function Pin1(; name)
        @independent_variables t
        sts = @variables v(t)=1.0 i(t)=1.0
        System(Equation[], t, sts, []; name = name)
    end
    @test string(Base.doc(Pin1)) == "Hey there, Pin1!\n"

    """
    Hey there, Pin2!
    """
    @component function Pin2(; name)
        @independent_variables t
        sts = @variables v(t)=1.0 i(t)=1.0
        System(Equation[], t, sts, []; name = name)
    end
    @test string(Base.doc(Pin2)) == "Hey there, Pin2!\n"
end

@testset "Issue#3016 Hierarchical indexing" begin
    @mtkmodel Inner begin
        @parameters begin
            p
        end
    end
    @mtkmodel Outer begin
        @components begin
            inner = Inner()
        end
        @variables begin
            x(t)
        end
        @equations begin
            x ~ inner.p
        end
    end

    @named outer = Outer()
    simp = mtkcompile(outer)

    @test sort(propertynames(outer)) == [:inner, :t, :x]
    @test propertynames(simp) == propertynames(outer)
    @test sort(propertynames(outer.inner)) == [:p, :t]
    @test propertynames(simp.inner) == propertynames(outer.inner)

    for sym in (:t, :x)
        @test_nowarn getproperty(simp, sym)
        @test_nowarn getproperty(outer, sym)
    end
    @test_nowarn simp.inner.p
    @test_nowarn outer.inner.p
    @test_throws ArgumentError simp.inner₊p
    @test_throws ArgumentError outer.inner₊p
end

@testset "`getproperty` on `mtkcompile(complete(sys))`" begin
    @mtkmodel Foo begin
        @variables begin
            x(t)
        end
    end
    @mtkmodel Bar begin
        @components begin
            foo = Foo()
        end
        @equations begin
            D(foo.x) ~ foo.x
        end
    end
    @named bar = Bar()
    cbar = complete(bar)
    ss = mtkcompile(cbar)
    @test isequal(cbar.foo.x, ss.foo.x)
end

@testset "Issue#3275: Metadata retained on `complete`" begin
    @variables x(t) y(t)
    @named inner = System(D(x) ~ x, t)
    @named outer = System(D(y) ~ y, t; systems = [inner], metadata = [Int => "test"])
    @test getmetadata(outer, Int, nothing) == "test"
    sys = complete(outer)
    @test getmetadata(sys, Int, nothing) == "test"
end

@testset "Causal connections generate causal equations" begin
    # test interpretation of `Equality` cset as causal connection
    @named input = RealInput()
    @named comp1 = System(Equation[], t; systems = [input])
    @named output = RealOutput()
    @named comp2 = System(Equation[], t; systems = [output])
    @named sys = System([connect(comp2.output, comp1.input)], t; systems = [comp1, comp2])
    eq = only(equations(expand_connections(sys)))
    # as opposed to `output.u ~ input.u`
    @test isequal(eq, comp1.input.u ~ comp2.output.u)

    # test causal ordering of true causal cset
    @named input = RealInput()
    @named comp1 = System(Equation[], t; systems = [input])
    @named output = RealOutput()
    @named comp2 = System(Equation[], t; systems = [output])
    @named sys = System([connect(comp2.output.u, comp1.input.u)], t; systems = [
        comp1, comp2])
    eq = only(equations(expand_connections(sys)))
    # as opposed to `output.u ~ input.u`
    @test isequal(eq, comp1.input.u ~ comp2.output.u)
end
