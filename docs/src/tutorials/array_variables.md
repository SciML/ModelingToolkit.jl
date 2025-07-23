# Working with Array Variables

This tutorial demonstrates how to create and use array (vector) variables in ModelingToolkit.jl. Array variables are useful when modeling systems with multiple similar components, such as multi-body systems, discretized PDEs, or any system with repeated structure.

## Basic Array Variable Syntax

The basic syntax for creating array variables is:

```@example array_vars
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

# Create a vector of 3 variables
@variables x[1:3](t)
```

This creates three variables: `x[1](t)`, `x[2](t)`, and `x[3](t)`. You can access individual elements:

```@example array_vars
x[1]  # First element
x[2]  # Second element
x[3]  # Third element
```

## Array Parameters

Similarly, you can create array parameters:

```@example array_vars
@parameters k[1:3] m[1:3]
```

## Building Systems with Array Variables

Here's an example of a mass-spring system with three masses connected in series:

```@example array_vars
# Define the system
n = 3  # number of masses

@variables x[1:n](t) v[1:n](t)
@parameters k[1:n] m[1:n] c[1:n]

# Create equations
eqs = Equation[]

# First mass (connected to wall and second mass)
push!(eqs, D(x[1]) ~ v[1])
push!(eqs, m[1] * D(v[1]) ~ -k[1] * x[1] + k[2] * (x[2] - x[1]) - c[1] * v[1])

# Middle masses
for i in 2:(n-1)
    push!(eqs, D(x[i]) ~ v[i])
    push!(eqs, m[i] * D(v[i]) ~ k[i] * (x[i-1] - x[i]) + k[i+1] * (x[i+1] - x[i]) - c[i] * v[i])
end

# Last mass
push!(eqs, D(x[n]) ~ v[n])
push!(eqs, m[n] * D(v[n]) ~ k[n] * (x[n-1] - x[n]) - c[n] * v[n])

@named mass_spring_system = System(eqs, t)
```

## Initial Conditions and Parameter Values for Arrays

When solving the system, you need to provide initial conditions and parameter values for each array element:

```@example array_vars
using OrdinaryDiffEq

# Initial conditions
u0 = [
    x[1] => 1.0,
    x[2] => 0.5,
    x[3] => 0.0,
    v[1] => 0.0,
    v[2] => 0.0,
    v[3] => 0.0
]

# Parameter values
p = [
    k[1] => 10.0,
    k[2] => 10.0,
    k[3] => 10.0,
    m[1] => 1.0,
    m[2] => 1.0,
    m[3] => 1.0,
    c[1] => 0.1,
    c[2] => 0.1,
    c[3] => 0.1
]

# Create and solve the problem
sys_simplified = mtkcompile(mass_spring_system)
prob = ODEProblem(sys_simplified, u0, (0.0, 10.0), p)
sol = solve(prob, Tsit5())

# Plot the solution
using Plots
plot(sol, idxs=[x[1], x[2], x[3]], labels=["x₁" "x₂" "x₃"])
```

## Multi-dimensional Arrays

You can also create multi-dimensional arrays:

```@example array_vars_2d
using ModelingToolkit
using ModelingToolkit: t_nounits as t

# 2D array of variables
@variables u[1:3, 1:3](t)

# Access elements
u[1, 1]  # Element at row 1, column 1
u[2, 3]  # Element at row 2, column 3
```

## Common Patterns and Tips

### Using Array Comprehensions

You can use Julia's array comprehension syntax to create equations more concisely:

```@example array_vars
n = 4
@variables y[1:n](t)
@parameters a[1:n] b[1:n]

# Create equations using comprehension
eqs = [D(y[i]) ~ a[i] * y[i] + b[i] for i in 1:n]

@named linear_system = System(eqs, t)
```

### Setting Default Values for Arrays

You can set default values when creating array variables:

```@example array_vars
@variables z[1:3](t) = [1.0, 2.0, 3.0]
@parameters p[1:3] = [0.1, 0.2, 0.3]
```

### Working with Variable-sized Arrays

When the array size needs to be a parameter, you can use symbolic arrays in `@mtkmodel`:

```@example array_vars_model
using ModelingToolkit

@mtkmodel VariableSizeModel begin
    @structural_parameters begin
        N = 5  # Array size
    end
    @parameters begin
        a[1:N]
    end
    @variables begin
        x(t)[1:N] = zeros(N)
    end
    @equations begin
        [D(x[i]) ~ -a[i] * x[i] for i in 1:N]...
    end
end
```

## Troubleshooting

### Common Error: Scalar Variable with Vector Initial Condition

If you accidentally provide a vector initial condition to a scalar variable:

```julia
@variables x(t)  # Scalar variable
prob = ODEProblem(sys, [x => [1.0, 2.0]], tspan, params)  # ERROR!
```

You'll get an error about `no method matching zero(::Type{Vector{Float64}})`. The solution is to either:

1. Use array variables: `@variables x[1:2](t)`
2. Provide scalar initial conditions: `[x => 1.0]`

### Performance Considerations

For large arrays, consider using:
- Static arrays (`@SVector`, `@SMatrix`) for small, fixed-size arrays
- Sparse arrays for systems with sparse connectivity
- Symbolic simplification with `mtkcompile` to eliminate redundant calculations

## Summary

Array variables in ModelingToolkit provide a powerful way to model systems with repeated structure. Key points:

- Use `@variables x[1:n](t)` syntax to create array variables
- Access elements with standard Julia indexing: `x[i]`
- Provide initial conditions and parameters for each array element
- Use array comprehensions for concise equation definition
- Multi-dimensional arrays are supported
- For variable-sized arrays, use `@mtkmodel` with structural parameters