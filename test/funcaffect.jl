using ModelingToolkit, Test, DifferentialEquations

@parameters t
@variables u(t)
D = Differential(t)

eqs = [ D(u) ~ -u ]

affect1!(integ, ctx; u) = integ.u[u] += 10

@named sys = ODESystem(eqs, t, [u], [], discrete_events=[[4.0, ]=>(affect1!, [u], [], nothing)])
prob = ODEProblem(sys, [u=> 10.0], (0, 10.0))
sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4+1][1] > 10.0

# context
function affect2!(integ, ctx; u)
    integ.u[u] += ctx[1]
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
function affect3!(integ, ctx; u, a)
    integ.u[u] += integ.p[a]
    integ.p[a] *= 2
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
function affect3!(integ, ctx; u, b)
    integ.u[u] += integ.p[b]
    integ.p[b] *= 2
end

@named sys = ODESystem(eqs, t, [u], [a], discrete_events=[[4.0, 8.0]=>(affect3!, [u], [a=> :b], nothing)])
prob = ODEProblem(sys, [u=> 10.0], (0, 10.0))

sol = solve(prob, Tsit5())
i4 = findfirst(==(4.0), sol[:t])
@test sol.u[i4+1][1] > 10.0
i8 = findfirst(==(8.0), sol[:t])
@test sol.u[i8+1][1] > 20.0

# same name
@test_throws ErrorException ODESystem(eqs, t, [u], [a], discrete_events=[[4.0, 8.0]=>(affect3!, [u], [a=> :u], nothing)]; name=:sys)




