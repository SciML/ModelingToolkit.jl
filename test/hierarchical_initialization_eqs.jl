using ModelingToolkit, OrdinaryDiffEq

t = only(@parameters(t))
D = Differential(t)
"""
A simple linear resistor model

![Resistor](https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcTpJkiEyqh-BRx27pvVH0GLZ4MP_D1oriBwJhnZdgIq7m17z9VKUWaW9MeNQAz1rTML2ho&usqp=CAU)
"""
@component function Resistor(; name, R = 1.0)
    systems = @named begin
        p = Pin()
        n = Pin()
    end
    vars = @variables begin
        v(t), [guess = 0.0]
        i(t), [guess = 0.0]
    end
    params = @parameters begin
        R = R, [description = "Resistance of this Resistor"]
    end
    eqs = [v ~ p.v - n.v
           i ~ p.i
           p.i + n.i ~ 0
           # Ohm's Law
           v ~ i * R]
    return ODESystem(eqs, t, vars, params; systems, name)
end
@connector Pin begin
    v(t)
    i(t), [connect = Flow]
end
@component function ConstantVoltage(; name, V = 1.0)
    systems = @named begin
        p = Pin()
        n = Pin()
    end
    vars = @variables begin
        v(t), [guess = 0.0]
        i(t), [guess = 0.0]
    end
    params = @parameters begin
        V = 10
    end
    eqs = [v ~ p.v - n.v
           i ~ p.i
           p.i + n.i ~ 0
           v ~ V]
    return ODESystem(eqs, t, vars, params; systems, name)
end

@component function Capacitor(; name, C = 1.0)
    systems = @named begin
        p = Pin()
        n = Pin()
    end
    vars = @variables begin
        v(t), [guess = 0.0]
        i(t), [guess = 0.0]
    end
    params = @parameters begin
        C = C
    end
    initialization_eqs = [
        v ~ 0
    ]
    eqs = [v ~ p.v - n.v
           i ~ p.i
           p.i + n.i ~ 0
           C * D(v) ~ i]
    return ODESystem(eqs, t, vars, params; systems, name, initialization_eqs)
end

@component function Ground(; name)
    systems = @named begin
        g = Pin()
    end
    eqs = [
        g.v ~ 0
    ]
    return ODESystem(eqs, t, [], []; systems, name)
end

@component function Inductor(; name, L = 1.0)
    systems = @named begin
        p = Pin()
        n = Pin()
    end
    vars = @variables begin
        v(t), [guess = 0.0]
        i(t), [guess = 0.0]
    end
    params = @parameters begin
        (L = L)
    end
    eqs = [v ~ p.v - n.v
           i ~ p.i
           p.i + n.i ~ 0
           L * D(i) ~ v]
    return ODESystem(eqs, t, vars, params; systems, name)
end

"""
This is an RLC model.  This should support markdown.  That includes
HTML as well.
"""
@component function RLCModel(; name)
    systems = @named begin
        resistor = Resistor(R = 100)
        capacitor = Capacitor(C = 0.001)
        inductor = Inductor(L = 1)
        source = ConstantVoltage(V = 30)
        ground = Ground()
    end
    initialization_eqs = [
        inductor.i ~ 0
    ]
    eqs = [connect(source.p, inductor.n)
           connect(inductor.p, resistor.p, capacitor.p)
           connect(resistor.n, ground.g, capacitor.n, source.n)]
    return ODESystem(eqs, t, [], []; systems, name, initialization_eqs)
end
"""Run model RLCModel from 0 to 10"""
function simple()
    @mtkbuild model = RLCModel()
    u0 = []
    prob = ODEProblem(model, u0, (0.0, 10.0))
    sol = solve(prob)
end
@test SciMLBase.successful_retcode(simple())

@named model = RLCModel()
@test length(ModelingToolkit.get_initialization_eqs(model)) == 1
syslist = ModelingToolkit.get_systems(model)
@test length(ModelingToolkit.get_initialization_eqs(syslist[1])) == 0
@test length(ModelingToolkit.get_initialization_eqs(syslist[2])) == 1
@test length(ModelingToolkit.get_initialization_eqs(syslist[3])) == 0
@test length(ModelingToolkit.get_initialization_eqs(syslist[4])) == 0
@test length(ModelingToolkit.get_initialization_eqs(syslist[5])) == 0
@test length(ModelingToolkit.initialization_equations(model)) == 2

u0 = []
prob = ODEProblem(structural_simplify(model), u0, (0.0, 10.0))
sol = solve(prob, Rodas5P())
@test length(sol.u[end]) == 2
@test length(equations(prob.f.initializeprob.f.sys)) == 0
@test length(unknowns(prob.f.initializeprob.f.sys)) == 0
