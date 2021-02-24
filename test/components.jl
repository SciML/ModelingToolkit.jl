using Test
using ModelingToolkit, StructuralTransformations, OrdinaryDiffEq

# Basic electric components
const t = Sym{ModelingToolkit.Parameter{Real}}(:t)
function Pin(name)
    @variables v(t) i(t)
    ODESystem(Equation[], t, [v, i], [], name=name, default_u0=[v=>1.0, i=>1.0])
end

function Ground(name)
    g = Pin(:g)
    eqs = [g.v ~ 0]
    ODESystem(eqs, t, [], [], systems=[g], name=name)
end

function ConstantVoltage(name; V = 1.0)
    val = V
    p = Pin(:p)
    n = Pin(:n)
    @parameters V
    eqs = [
           V ~ p.v - n.v
           0 ~ p.i + n.i
          ]
    ODESystem(eqs, t, [], [V], systems=[p, n], default_p=Dict(V => val), name=name)
end

function Resistor(name; R = 1.0)
    val = R
    p = Pin(:p)
    n = Pin(:n)
    @variables v(t)
    @parameters R
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           v ~ p.i * R
          ]
    ODESystem(eqs, t, [v], [R], systems=[p, n], default_p=Dict(R => val), name=name)
end

function Capacitor(name; C = 1.0)
    val = C
    p = Pin(:p)
    n = Pin(:n)
    @variables v(t)
    @parameters C
    D = Differential(t)
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           D(v) ~ p.i / C
          ]
    ODESystem(eqs, t, [v], [C], systems=[p, n], default_p=Dict(C => val), name=name)
end

R = 1.0
C = 1.0
V = 1.0
resistor = Resistor(:resistor, R=R)
capacitor = Capacitor(:capacitor, C=C)
source = ConstantVoltage(:source, V=V)
ground = Ground(:ground)

function connect(ps...)
    eqs = [
           0 ~ sum(p->p.i, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps)-1
        push!(eqs, ps[i].v ~ ps[i+1].v)
    end

    return eqs
end
rc_eqs = [
          connect(source.p, resistor.p)
          connect(resistor.n, capacitor.p)
          connect(capacitor.n, source.n, ground.g)
         ]

rc_model = ODESystem(rc_eqs, t, systems=[resistor, capacitor, source, ground], name=:rc)
sys = alias_elimination(rc_model)
@test ModelingToolkit.default_p(sys) == Dict(
                                             capacitor.C => 1.0,
                                             source.V => 1.0,
                                             resistor.R => 1.0,
                                            )
u0 = [
      capacitor.v => 0.0
      capacitor.p.i => 0.0
      resistor.v => 0.0
     ]
prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas5())
#using Plot
#plot(sol)
