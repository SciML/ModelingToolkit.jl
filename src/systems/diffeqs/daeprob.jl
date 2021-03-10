# WIP 

# "The goal is to make a DAEProblem and DAEFunction constructor"
#    - constructors that take sys::AbstractODESystem though,
#  since there are already constructors that take an f(du, u, p, t) or f!(out, du, u, p, t)

# scratchpad
# using ModelingToolkit, DiffEqBase, OrdinaryDiffEq, DifferentialEquations

# function f(out,du,u,p,t)
#   out[1] = - 0.04u[1]              + 1e4*u[2]*u[3] - du[1]
#   out[2] = + 0.04u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3] - du[2]
#   out[3] = u[1] + u[2] + u[3] - 1.0
# end

# u₀ = [1.0, 0, 0]
# du₀ = [-0.04, 0.04, 0.0]
# tspan = (0.0,100.0)

# differential_vars = [true,true,false]
# prob = DAEProblem(f,du₀,u₀,tspan,differential_vars=differential_vars)
# solve(prob, Tsit5()) # ERROR: You cannot use an ODE Algorithm with a DAEProblem
# sys = modelingtoolkitize(prob) # fails
# sol = solve(prob) # figures out it needs a DAE alg


# @parameters t
# @variables u[1:3](t)
# D = Differential(t)
# # symbolic representation of rober
# e1 = D(u[1]) ~ -0.04u[1] + 1e4 * u[2] * u[3]
# e2 = D(u[2]) ~ 0.04u[1] - 3e7 * u[2]^2 - 1e4 * u[2] * u[3]
# e3 = 1.0 ~ u[1] + u[2] + u[3]

# eqs = [e1, e2, e3]


function eq_to_num(eq::Equation)::Num
    getfield(eq, :rhs) - getfield(eq, :lhs)
end

"""
    implicitize(eq::Equation)::Equation

transform an equation of the form y = f(x) 
into f(x) - y = 0
"""
function implicitize(eq::Equation)::Equation
    eq_to_num(eq) ~ 0 
end

# ieqs = implicitize.(eqs)
