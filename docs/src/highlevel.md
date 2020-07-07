# High Level API

The high-level API allows modelers to interactively build models in a symbolic
manner. It is designed as a semi-DSL for easily building large complex models
and manipulating the models to generate optimal forms to be used in numerical
methods.

## High-Level API Documentation

```@docs
@parameters
@variables
@derivatives
Base.:~(::Expression, ::Expression)
modelingtoolkitize
```

## Differentiation Functions

```@docs
ModelingToolkit.gradient
ModelingToolkit.jacobian
ModelingToolkit.hessian
```

## Additional High-Level Explanations and Tips

### The Auto-Detecting System Constructors

For the high-level interface, the system constructors, such as `ODESystem`, have
high-level constructors, which just take in the required equations and automatically
parse the expressions to figure out the states and parameters of the system.
The following high-level constructors exist:

```julia
ODESystem(eqs)
NonlinearSystem(eqs)
```

### Direct Tracing

Because ModelingToolkit expressions respect Julia semantics, one way
to generate symbolic expressions is to simply place ModelingToolkit
variables as inputs into existing Julia code. For example, the following
uses the standard Julia function for the Lorenz equations to generate
the symbolic expression for the Lorenz equations:

```julia
function lorenz(du,u,p,t)
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
@variables t u[1:3](t) du[1:3](t)
@parameters p[1:3]
lorenz(du,u,p,t)
du
```

```julia
3-element Array{Operation,1}:
                 10.0 * (u₂(t) - u₁(t))
         u₁(t) * (28.0 - u₃(t)) - u₂(t)
u₁(t) * u₂(t) - 2.6666666666666665 * u₃(t)
```

Or similarly:

```julia
@variables t x(t) y(t) z(t) dx(t) dy(t) dz(t)
@parameters σ ρ β
du = [dx,dy,dz]
u = [x,y,z]
p = [σ,ρ,β]
lorenz(du,u,p,t)
du
```

```julia
3-element Array{Operation,1}:
                10.0 * (y(t) - x(t))
         x(t) * (28.0 - z(t)) - y(t)
x(t) * y(t) - 2.6666666666666665 * z(t)
```

### Intermediate Calculations

The system building functions can handle intermediate calculations by simply
defining and using an `Operation` of `Variable`s. For example:

```julia
@variables x y z
@parameters σ ρ β
a = y - x
eqs = [0 ~ σ*a,
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs, [x,y,z], [σ,ρ,β])
nlsys_func = generate_function(ns)[2] # second is the inplace version
```

expands to:

```julia
:((var"##MTIIPVar#368", var"##MTKArg#365", var"##MTKArg#366")->begin
          @inbounds begin
                  let (x, y, z, σ, ρ, β) = (var"##MTKArg#365"[1], var"##MTKArg#365"[2], var"##MTKArg#365"[3], var"##MTKArg#366"[1], var"##MTKArg#366"[2], var"##MTKArg#366"[3])
                      var"##MTIIPVar#368"[1] = (*)(σ, (-)(y, x))
                      var"##MTIIPVar#368"[2] = (-)((*)(x, (-)(ρ, z)), y)
                      var"##MTIIPVar#368"[3] = (-)((*)(x, y), (*)(β, z))
                  end
              end
          nothing
      end)
```

In addition, the Jacobian calculations take into account intermediate variables
to appropriately handle them.
