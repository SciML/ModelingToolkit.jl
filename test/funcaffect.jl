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


