using Test
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit.BipartiteGraphs
using ModelingToolkit.StructuralTransformations
include("../examples/rc_model.jl")

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

include("../examples/rc_model.jl")

completed_rc_model = complete(rc_model)
@test isequal(completed_rc_model.resistor.n.i, resistor.n.i)
@test ModelingToolkit.n_extra_equations(capacitor) == 2
@test length(equations(structural_simplify(rc_model, allow_parameter = false))) == 2
sys = structural_simplify(rc_model)
@test length(equations(sys)) == 1
check_contract(sys)
@test !isempty(ModelingToolkit.defaults(sys))
u0 = [capacitor.v => 0.0
      capacitor.p.i => 0.0
      resistor.v => 0.0]
prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())
check_rc_sol(sol)
prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())
check_rc_sol(sol)

# test uniform default
defs = [capacitor.v => 0.0
        capacitor.p.i => 0.0
        resistor.v => 0.0
        resistor.R => 1234.0]
prob2 = ODEProblem(sys, defs, (0, 10.0))
@test 1234.0 in prob2.p
prob3 = ODAEProblem(sys, defs, (0, 10.0))
@test 1234.0 in prob3.p

# https://discourse.julialang.org/t/using-optimization-parameters-in-modelingtoolkit/82099
let
    @parameters param_r1 param_c1
    @named resistor = Resistor(R = param_r1)
    @named capacitor = Capacitor(C = param_c1)
    @named source = ConstantVoltage(V = 1.0)
    @named ground = Ground()

    rc_eqs = [connect(source.p, resistor.p)
              connect(resistor.n, capacitor.p)
              connect(capacitor.n, source.n)
              connect(capacitor.n, ground.g)]

    @named _rc_model = ODESystem(rc_eqs, t)
    @named rc_model = compose(_rc_model,
                              [resistor, capacitor, source, ground])
    sys = structural_simplify(rc_model)
    u0 = [
        capacitor.v => 0.0,
    ]

    params = [param_r1 => 1.0, param_c1 => 1.0]
    tspan = (0.0, 10.0)

    prob = ODAEProblem(sys, u0, tspan, params)
    @test solve(prob, Tsit5()).retcode == :Success
end

let
    # 1478
    @named resistor2 = Resistor(R = R)
    rc_eqs2 = [connect(source.p, resistor.p)
               connect(resistor.n, resistor2.p)
               connect(resistor2.n, capacitor.p)
               connect(capacitor.n, source.n)
               connect(capacitor.n, ground.g)]

    @named _rc_model2 = ODESystem(rc_eqs2, t)
    @named rc_model2 = compose(_rc_model2,
                               [resistor, resistor2, capacitor, source, ground])
    sys2 = structural_simplify(rc_model2)
    prob2 = ODAEProblem(sys2, u0, (0, 10.0))
    sol2 = solve(prob2, Tsit5())
    @test sol2[source.p.i] == sol2[rc_model2.source.p.i] == -sol2[capacitor.i]
end

# Outer/inner connections
function rc_component(; name, R = 1, C = 1)
    @parameters R=R C=C
    @named p = Pin()
    @named n = Pin()
    @named resistor = Resistor(R = R) # test parent scope default of @named
    @named capacitor = Capacitor(C = ParentScope(C))
    eqs = [connect(p, resistor.p);
           connect(resistor.n, capacitor.p);
           connect(capacitor.n, n)]
    @named sys = ODESystem(eqs, t, [], [R, C])
    compose(sys, [p, n, resistor, capacitor]; name = name)
end

@named ground = Ground()
@named source = ConstantVoltage(V = 1)
@named rc_comp = rc_component()
eqs = [connect(source.p, rc_comp.p)
       connect(source.n, rc_comp.n)
       connect(source.n, ground.g)]
@named sys′ = ODESystem(eqs, t)
@named sys_inner_outer = compose(sys′, [ground, source, rc_comp])
@test_nowarn show(IOBuffer(), MIME"text/plain"(), sys_inner_outer)
expand_connections(sys_inner_outer, debug = true)
sys_inner_outer = structural_simplify(sys_inner_outer)
@test !isempty(ModelingToolkit.defaults(sys_inner_outer))
u0 = [rc_comp.capacitor.v => 0.0
      rc_comp.capacitor.p.i => 0.0
      rc_comp.resistor.v => 0.0]
prob = ODEProblem(sys_inner_outer, u0, (0, 10.0), sparse = true)
sol_inner_outer = solve(prob, Rodas4())
@test sol[capacitor.v] ≈ sol_inner_outer[rc_comp.capacitor.v]

u0 = [
    capacitor.v => 0.0,
]
prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Tsit5())

@test sol[resistor.p.i] == sol[capacitor.p.i]
@test sol[resistor.n.i] == -sol[capacitor.p.i]
@test sol[capacitor.n.i] == -sol[capacitor.p.i]
@test iszero(sol[ground.g.i])
@test iszero(sol[ground.g.v])
@test sol[resistor.v] == sol[source.p.v] - sol[capacitor.p.v]
#using Plots
#plot(sol)

include("../examples/serial_inductor.jl")
sys = structural_simplify(ll_model)
@test length(equations(sys)) == 2
check_contract(sys)
u0 = states(sys) .=> 0
@test_nowarn ODEProblem(sys, u0, (0, 10.0))
@test_nowarn ODAEProblem(sys, u0, (0, 10.0))
prob = DAEProblem(sys, Differential(t).(states(sys)) .=> 0, u0, (0, 0.5))
@test_nowarn sol = solve(prob, DFBDF())

@variables t x1(t) x2(t) x3(t) x4(t)
D = Differential(t)
@named sys1_inner = ODESystem([D(x1) ~ x1], t)
@named sys1_partial = compose(ODESystem([D(x2) ~ x2], t; name = :foo), sys1_inner)
@named sys1 = extend(ODESystem([D(x3) ~ x3], t; name = :foo), sys1_partial)
@named sys2 = compose(ODESystem([D(x4) ~ x4], t; name = :foo), sys1)
@test_nowarn sys2.sys1.sys1_inner.x1 # test the correct nesting

# compose tests
@parameters t

function record_fun(; name)
    pars = @parameters a=10 b=100
    ODESystem(Equation[], t, [], pars; name)
end

function first_model(; name)
    @named foo = record_fun()

    defs = Dict()
    defs[foo.a] = 3
    defs[foo.b] = 300
    pars = @parameters x=2 y=20
    compose(ODESystem(Equation[], t, [], pars; name, defaults = defs), foo)
end
@named goo = first_model()
@unpack foo = goo
@test ModelingToolkit.defaults(goo)[foo.a] == 3
@test ModelingToolkit.defaults(goo)[foo.b] == 300

#=
model Circuit
  Ground ground;
  Load load;
  Resistor resistor;
equation
  connect(load.p , ground.p);
  connect(resistor.p, ground.p);
end Circuit;
model Load
  extends TwoPin;
  Resistor resistor;
equation
  connect(p, resistor.p);
  connect(resistor.n, n);
end Load;
=#

function Load(; name)
    R = 1
    @named p = Pin()
    @named n = Pin()
    @named resistor = Resistor(R = R)
    eqs = [connect(p, resistor.p);
           connect(resistor.n, n)]
    @named sys = ODESystem(eqs, t)
    compose(sys, [p, n, resistor]; name = name)
end

function Circuit(; name)
    R = 1
    @named ground = Ground()
    @named load = Load()
    @named resistor = Resistor(R = R)
    eqs = [connect(load.p, ground.g);
           connect(resistor.p, ground.g)]
    @named sys = ODESystem(eqs, t)
    compose(sys, [ground, resistor, load]; name = name)
end

@named foo = Circuit()
@test structural_simplify(foo) isa ModelingToolkit.AbstractSystem

# BLT tests
using LinearAlgebra
function parallel_rc_model(i; name, source, ground, R, C)
    resistor = HeatingResistor(name = Symbol(:resistor, i), R = R)
    capacitor = Capacitor(name = Symbol(:capacitor, i), C = C)
    heat_capacitor = HeatCapacitor(name = Symbol(:heat_capacitor, i))

    rc_eqs = [connect(source.p, resistor.p)
              connect(resistor.n, capacitor.p)
              connect(capacitor.n, source.n, ground.g)
              connect(resistor.h, heat_capacitor.h)]

    compose(ODESystem(rc_eqs, t, name = Symbol(name, i)),
            [resistor, capacitor, source, ground, heat_capacitor])
end
V = 2.0
@named source = ConstantVoltage(V = V)
@named ground = Ground()
N = 50
Rs = 10 .^ range(0, stop = -4, length = N)
Cs = 10 .^ range(-3, stop = 0, length = N)
rc_systems = map(1:N) do i
    parallel_rc_model(i; name = :rc, source = source, ground = ground, R = Rs[i], C = Cs[i])
end;
@variables E(t) = 0.0
eqs = [
    D(E) ~ sum(((i, sys),) -> getproperty(sys, Symbol(:resistor, i)).h.Q_flow,
               enumerate(rc_systems)),
]
@named _big_rc = ODESystem(eqs, t, [E], [])
@named big_rc = compose(_big_rc, rc_systems)
ts = TearingState(expand_connections(big_rc))
@test istriu(but_ordered_incidence(ts)[1])
