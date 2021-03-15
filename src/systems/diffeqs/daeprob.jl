using ModelingToolkit, DiffEqBase, OrdinaryDiffEq, DifferentialEquations

function eq_to_num(eq::Equation)
    lhs(eq) - rhs(eq)
end

lhs(eq::Equation) = getfield(eq, :lhs)
rhs(eq::Equation) = getfield(eq, :rhs)

"""
    implicitize(eq::Equation)::Equation

transform an equation of the form y = f(x) 
into f(x) - y = 0
"""
function implicitize(eq::Equation)
    eq_to_num(eq) ~ 0 
end

@parameters t r
@variables u[1:3](t)
D = Differential(t)
# symbolic representation of rober
eqs = [
    D(u[1]) ~ -u[1] + r * u[2] * u[3],
    D(u[2]) ~ u[1] - u[2]^2 - u[2] * u[3],
    1.0 ~ u[1] + u[2] + u[3]
]
nums = eq_to_num.(eqs)
@variables du[1:3] # temporaries instead of D.(u)
p = [r]
f_expr = build_function(nums, du, u, p, t) # du,u,p,t
fiip = f_expr[2]

myf = eval(fiip)

out = zeros(3)
u₀ = [1.0, 0, 0]
du₀ = [-0.04, 0.04, 0.0]
tspan = (0.0, 100.0)
ps = rand(1)
myf(out,du0,u0,ps, tspan)

Base.remove_linenums!(fiip)

:(function (var"##out#264", var"##arg#261", var"##arg#262", var"##arg#263", t)
      let du₁ = var"##arg#261"[1], du₂ = var"##arg#261"[2], du₃ = var"##arg#261"[3], var"u₁(t)" = var"##arg#262"[1], var"u₂(t)" = var"##arg#262"[2], var"u₃(t)" = var"##arg#262"[3], r = var"##arg#263"[1]
          #= C:\Users\Anand\.julia\packages\SymbolicUtils\V6S2E\src\code.jl:310 =# @inbounds begin
                  var"##out#264"[1] = (+)((*)(-1, r, var"u₂(t)", var"u₃(t)"), var"u₁(t)", (Differential(t))(var"u₁(t)"))
                  var"##out#264"[2] = (+)((Differential(t))(var"u₂(t)"), (*)(-1, var"u₁(t)"), (^)(var"u₂(t)", 2), (*)(var"u₂(t)", var"u₃(t)"))
                  var"##out#264"[3] = (+)(1.0, (*)(-1, var"u₁(t)"), (*)(-1, var"u₂(t)"), (*)(-1, var"u₃(t)"))
                  nothing
              end
      end
  end)