using OrdinaryDiffEq, ModelingToolkit, Test
@parameters t
@variables y[1:3](t)
@parameters k[1:3]
D = Differential(t)

eqs = [D(y[1]) ~ -k[1]*y[1] + k[3]*y[2]*y[3],
       D(y[2]) ~  k[1]*y[1] - k[3]*y[2]*y[3] - k[2]*y[2]^2,
       0 ~  y[1] + y[2] + y[3] - 1]

sys = ODESystem(eqs,t,y,k)
@test_throws ArgumentError ODESystem(eqs,y[1])
M = calculate_massmatrix(sys)
@test M == [1 0 0
            0 1 0
            0 0 0]

f = ODEFunction(sys)
prob_mm = ODEProblem(f,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))
sol = solve(prob_mm,Rodas5(),reltol=1e-8,abstol=1e-8)

function rober(du,u,p,t)
  y₁,y₂,y₃ = u
  k₁,k₂,k₃ = p
  du[1] = -k₁*y₁ + k₃*y₂*y₃
  du[2] =  k₁*y₁ - k₃*y₂*y₃ - k₂*y₂^2
  du[3] =  y₁ + y₂ + y₃ - 1
  nothing
end
f = ODEFunction(rober,mass_matrix=M)
prob_mm2 = ODEProblem(f,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))
sol2 = solve(prob_mm2,Rodas5(),reltol=1e-8,abstol=1e-8,tstops=sol.t,adaptive=false)

# MTK expression are canonicalized, so the floating point numbers are slightly
# different
@test Array(sol) ≈ Array(sol2)
