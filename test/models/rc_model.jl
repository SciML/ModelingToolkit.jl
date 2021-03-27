include("electrical_components.jl")

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
