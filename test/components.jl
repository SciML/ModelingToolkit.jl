using Test
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit.BipartiteGraphs

function check_contract(sys)
    s = structure(sys)
    @unpack fullvars, graph = s
    eqs = equations(sys)
    var2idx = Dict(enumerate(fullvars))
    for (i, eq) in enumerate(eqs)
        actual = union(ModelingToolkit.vars(eq.lhs), ModelingToolkit.vars(eq.rhs))
        actual = filter(!ModelingToolkit.isparameter, collect(actual))
        current = Set(fullvars[𝑠neighbors(graph, i)])
        @test isempty(setdiff(actual, current))
    end
end

include("../examples/rc_model.jl")

sys = structural_simplify(rc_model)
check_contract(sys)
@test !isempty(ModelingToolkit.defaults(sys))
u0 = [
      capacitor.v => 0.0
      capacitor.p.i => 0.0
      resistor.v => 0.0
     ]
prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())

@test sol[resistor.p.i] == sol[capacitor.p.i]
@test sol[resistor.n.i] == -sol[capacitor.p.i]
@test sol[capacitor.n.i] == -sol[capacitor.p.i]
@test iszero(sol[ground.g.i])
@test iszero(sol[ground.g.v])
@test sol[resistor.v] == sol[source.p.v] - sol[capacitor.p.v]

# Outer/inner connections
function rc_component(;name)
    R = 1
    C = 1
    @named p = Pin()
    @named n = Pin()
    @named resistor = Resistor(R=R)
    @named capacitor = Capacitor(C=C)
    eqs = [
           connect(p, resistor.p);
           connect(resistor.n, capacitor.p);
           connect(capacitor.n, n);
          ]
    @named sys = ODESystem(eqs, t)
    compose(sys, [p, n, resistor, capacitor]; name=name)
end

@named ground = Ground()
@named source = ConstantVoltage(V=1)
@named rc_comp = rc_component()
eqs = [
       connect(source.p, rc_comp.p)
       connect(source.n, rc_comp.n)
       connect(source.n, ground.g)
      ]
@named sys′ = ODESystem(eqs, t)
@named sys_inner_outer = compose(sys′, [ground, source, rc_comp])
expand_connections(sys_inner_outer, debug=true)
sys_inner_outer = structural_simplify(sys_inner_outer)
@test !isempty(ModelingToolkit.defaults(sys_inner_outer))
u0 = [
      rc_comp.capacitor.v => 0.0
      rc_comp.capacitor.p.i => 0.0
      rc_comp.resistor.v => 0.0
     ]
prob = ODEProblem(sys_inner_outer, u0, (0, 10.0))
sol_inner_outer = solve(prob, Rodas4())
@test sol[capacitor.v] ≈ sol_inner_outer[rc_comp.capacitor.v]

u0 = [
      capacitor.v => 0.0
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
check_contract(sys)
u0 = [
      inductor1.i => 0.0
      inductor2.i => 0.0
      inductor2.v => 0.0
     ]
prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())

prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Tsit5())

@variables t x1(t) x2(t) x3(t) x4(t)
D = Differential(t)
@named sys1_inner = ODESystem([D(x1) ~ x1], t)
@named sys1_partial = compose(ODESystem([D(x2) ~ x2], t; name=:foo), sys1_inner)
@named sys1 = extend(ODESystem([D(x3) ~ x3], t; name=:foo), sys1_partial)
@named sys2 = compose(ODESystem([D(x4) ~ x4], t; name=:foo), sys1)
@test_nowarn sys2.sys1.sys1_inner.x1 # test the correct nesting


# compose tests
@parameters t

function record_fun(;name)
    pars = @parameters a=10 b=100
    ODESystem(Equation[], t, [], pars; name)
end

function first_model(;name)
    @named foo=record_fun()

    defs = Dict()
    defs[foo.a] = 3
    defs[foo.b] = 300
    pars = @parameters x=2 y=20
    compose(ODESystem(Equation[], t, [], pars; name, defaults=defs), foo)
end
@named goo = first_model()
@unpack foo = goo
@test ModelingToolkit.defaults(goo)[foo.a] == 3
@test ModelingToolkit.defaults(goo)[foo.b] == 300
