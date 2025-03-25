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
D = ModelingToolkit.Differential(t)
vars2_arr = [a, b, c, d]

vars2_quants =[DynamicQuantities.Quantity(0, mass=1, length=-3), DynamicQuantities.Quantity(0, mass=1, length=-1, time=-1), DynamicQuantities.Quantity(0, length=1, time=-1), DynamicQuantities.Quantity(0, length=1)]

pi_terms = buckinghumFun(vars2_quants, vars2_arr)
original_equations = Dict(D(a) => α*t + b, D(b) => c*d/a)

@test isequal(pi_terms[1], (a*c*d)/b)
@test isequal(transform_eqs_from_pi_terms(pi_terms, original_equations, t)[1], ((b + t*α)*c*d) / b + (-(c^2)*(d^2)) / (b^2))

