using ModelingToolkit, Test, DifferentialEquations

@parameters t
@variables u(t)
D = Differential(t)

eqs = [ D(u) ~ -u ]

affect1!(integ, u, p, ctx) = integ.u[u.u] += 10

@named sys = ODESystem(eqs, t, [u], [], discrete_events=[[4.0, ]=>(affect1!, [u], [], nothing)])
prob = ODEProblem(sys, [u=> 10.0], (0, 10.0))
sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4+1][1] > 10.0

# context
function affect2!(integ, u, p, ctx)
    integ.u[u.u] += ctx[1]
    ctx[1] *= 2
end
ctx1 = [10.0, ]
@named sys = ODESystem(eqs, t, [u], [], discrete_events=[[4.0, 8.0]=>(affect2!, [u], [], ctx1)])
prob = ODEProblem(sys, [u=> 10.0], (0, 10.0))
sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4+1][1] > 10.0
i8 = findfirst(==(8.0), sol[:t])
@test sol.u[i8+1][1] > 20.0
@test ctx1[1] == 40.0

# parameter
function affect3!(integ, u, p, ctx)
    integ.u[u.u] += integ.p[p.a]
    integ.p[p.a] *= 2
end

@parameters a = 10.0
@named sys = ODESystem(eqs, t, [u], [a], discrete_events=[[4.0, 8.0]=>(affect3!, [u], [a], nothing)])
prob = ODEProblem(sys, [u=> 10.0], (0, 10.0))

sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4+1][1] > 10.0
i8 = findfirst(==(8.0), sol[:t])
@test sol.u[i8+1][1] > 20.0

# rename parameter
function affect3!(integ, u, p, ctx)
    integ.u[u.u] += integ.p[p.b]
    integ.p[p.b] *= 2
end

@named sys = ODESystem(eqs, t, [u], [a], discrete_events=[[4.0, 8.0]=>(affect3!, [u], [a=> :b], nothing)])
prob = ODEProblem(sys, [u=> 10.0], (0, 10.0))

sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4+1][1] > 10.0
i8 = findfirst(==(8.0), sol[:t])
@test sol.u[i8+1][1] > 20.0

# same name
@variables v(t)
@test_throws ErrorException ODESystem(eqs, t, [u], [a], discrete_events=[[4.0, 8.0]=>(affect3!, [u, v => :u], [a], nothing)]; name=:sys)

@test_nowarn ODESystem(eqs, t, [u], [a], discrete_events=[[4.0, 8.0]=>(affect3!, [u], [a => :u], nothing)]; name=:sys)

@named resistor = ODESystem(D(v) ~ v, t, [v], [])

# nested namespace
ctx = [0]
function affect4!(integ, u, p, ctx)
    ctx[1] += 1
    @test u.resistor₊v == 1
end
s1 = compose(ODESystem(Equation[], t, [], [], name=:s1, discrete_events=1.0=>(affect4!, [resistor.v], [], ctx)), resistor)
s2 = structural_simplify(s1)
prob = ODEProblem(s2, [resistor.v=> 10.0], (0, 2.01))
sol = solve(prob, Tsit5())
@test ctx[1] == 2

include("../examples/rc_model.jl")

function affect5!(integ, u,p,ctx)
    @test integ.u[u.capacitor₊v] ≈ 0.3
    integ.p[p.C] *= 200
end

@named rc_model = ODESystem(rc_eqs, t, continuous_events=[[capacitor.v ~ 0.3]=>(affect5!, [capacitor.v], [capacitor.C => :C], nothing)])
rc_model = compose(rc_model, [resistor, capacitor, source, ground])

sys = structural_simplify(rc_model)
u0 = [capacitor.v => 0.0
      capacitor.p.i => 0.0
      resistor.v => 0.0]

prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())
@test all(sol[rc_model.capacitor.v] .< 0.4)

function affect6!(integ, u,p,ctx)
    @test integ.u[u.v] ≈ 0.3
    integ.p[p.C] *= 200
end

function Capacitor(; name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C = C
    D = Differential(t)
    eqs = [
        D(v) ~ i / C,
    ]
    extend(ODESystem(eqs, t, [], ps; name = name, continuous_events=[[v ~ 0.3]=>(affect6!, [v], [C], nothing)]), oneport)
end

# hyerarchical - result should be identical

@named capacitor2 = Capacitor(C = C)

rc_eqs2 = [connect(source.p, resistor.p)
          connect(resistor.n, capacitor2.p)
          connect(capacitor2.n, source.n)
          connect(capacitor2.n, ground.g)]

@named rc_model2 = ODESystem(rc_eqs2, t)
rc_model2 = compose(rc_model2, [resistor, capacitor2, source, ground])

sys2 = structural_simplify(rc_model2)
u0 = [capacitor2.v => 0.0
      capacitor2.p.i => 0.0
      resistor.v => 0.0]

prob2 = ODEProblem(sys2, u0, (0, 10.0))
sol2 = solve(prob2, Rodas4())
@test all(sol2[rc_model2.capacitor2.v] .== sol[rc_model.capacitor.v])
