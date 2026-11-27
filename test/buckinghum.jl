include("../src/systems/buckinghum.jl")
using DynamicQuantities
using LinearAlgebra
using Test

@variables f, d, v, ρ, μ
vars1_arr = [ρ, d, v, μ, f]

vars1_quants = [DynamicQuantities.Quantity(0, mass=1, length=-1, time=-1), DynamicQuantities.Quantity(0, length=1) , DynamicQuantities.Quantity(0, length=1, time=-1),  DynamicQuantities.Quantity(0, mass=1, length=1, time=-2), DynamicQuantities.Quantity(0, mass=1,  length=-1, time=-1)]

@test isequal(buckinghumFun(vars1_quants, vars1_arr)[1], μ/(d*v*ρ))
@test isequal(buckinghumFun(vars1_quants, vars1_arr)[2], f/(ρ))


@parameters t, α
@variables  a(t), b(t), c(t), d(t)
@variables π1(t) π2(t) π3(t) π4(t) π5(t)
D = ModelingToolkit.Differential(t)
vars2_arr = [a, b, c, d]
vars2_quants =[DynamicQuantities.Quantity(0, mass=1, length=-3), DynamicQuantities.Quantity(0, mass=1, length=-1, time=-1), DynamicQuantities.Quantity(0, length=1, time=-1), DynamicQuantities.Quantity(0, length=1)]

original_equations_map = Dict(D(a) => α*t + b, D(b) => c*d/a)
original_equations = [k ~ v for (k,v) in original_equations_map]

@named sys = ODESystem(original_equations, t, vars2_arr, [α];defaults=Dict(α => 0.5))
pi_terms = buckinghumFun(vars2_quants, vars2_arr)
new_sys = transform_sys(pi_terms, sys, [π1, π2, π3, π4, π5])

@test isequal(pi_terms[1], (a*c*d)/b)
@test isequal(equations(new_sys)[1].rhs ,((b + t*α)*c*d) / b + (-(c^2)*(d^2)) / (b^2))

