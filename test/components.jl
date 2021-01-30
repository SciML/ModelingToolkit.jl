using ModelingToolkit
using Test

# Basic electric components
@parameters t
function Pin(name)
    @variables v(t) i(t)
    ODESystem(Equation[], t, [v, i], [], name=name, default_u0=[v=>1.0, i=>1.0])
end

function Ground(name)
    p = Pin(name)
    eqs = [p.v ~ 0]
    ODESystem(eqs, t, [], [], systems=[p], name=name)
end

function OnePort(name)
    p = Pin(:p)
    n = Pin(:n)
    @variables v(t) i(t)
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           i ~ p.i
          ]
    v, i, ODESystem(eqs, t, [v, i], [], systems=[p, n], name=name)
end

function ConstantVoltage(name; V = 1.0)
    val = V
    v, i, op = OnePort(name)
    @parameters V
    push!(op.eqs, v ~ V)
    V = ModelingToolkit.value(V)
    push!(op.ps, V)
    op.default_p[V] = val
    op
end

function Resistor(name; R = 1.0)
    val = R
    v, i, op = OnePort(name)
    @parameters R
    push!(op.eqs, v ~ R * i)
    R = ModelingToolkit.value(R)
    push!(op.ps, R)
    op.default_p[R] = val
    op
end

function Capacitor(name; C = 1.0)
    val = C
    v, i, op = OnePort(name)
    @parameters C
    D = Differential(t)
    push!(op.eqs, D(v) ~ i / C)
    C = ModelingToolkit.value(C)
    push!(op.ps, C)
    op.default_p[C] = val
    op
end

function Inductor(name; L = 1.0)
    val = L
    v, i, op = OnePort(name)
    @parameters L
    D = Differential(t)
    push!(op.eqs, D(i) ~ v / L)
    L = ModelingToolkit.value(L)
    push!(op.ps, L)
    op.default_p[L] = val
    op
end

R1 = Resistor(:R1)
C1 = Capacitor(:C1)
I1 = Inductor(:I1)
@test isequal(R1.p.i, (@variables R1₊p₊i(t))[1])
@parameters R
@test ModelingToolkit.get_default_p(flatten(R1)) == Dict(R => 1.0)
@variables n₊v(t) p₊v(t) n₊i(t) p₊i(t)
@test ModelingToolkit.get_default_u0(flatten(R1)) == Dict(
                                                          n₊v => 1.0,
                                                          p₊v => 1.0,
                                                          n₊i => 1.0,
                                                          p₊i => 1.0,
                                                         )
