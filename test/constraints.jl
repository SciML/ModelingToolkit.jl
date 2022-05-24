using ModelingToolkit, DiffEqBase, LinearAlgebra

# Define some variables
@parameters t x y
@variables u(..)

ConstrainedEquation([x ~ 0, y < 1 / 2], u(t, x, y) ~ x + y^2)
