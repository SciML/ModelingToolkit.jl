# High Level API

The high level API allows modelers to interactively build models in a symbolic
manner. It is designed as a semi-DSL for easily building large complex models
and manipulating the models to generate optimal forms for utilizing in numerical
methods.

## Examples

### Example 1: Symbolically Building an ODEProblem for DifferentialEquations.jl

Let's build an ODE. First we define some variables. In a differential equation
system, we need to differentiate between our (dependent) variables
and parameters. Therefore we label them as follows:

```julia
using ModelingToolkit

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t
```

Then we build the system:

```julia
eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
```

Each operation builds an `Operation` type, and thus `eqs` is an array of
`Operation` and `Variable`s. This holds a tree of the full system that can be
analyzed by other programs. We can turn this into a `ODESystem` via:

```julia
de = ODESystem(eqs)
```

This `ODESystem` can then be used to generate an `ODEProblem` by supplying the
constructor with a map from the states of the system to their initial condition
values and from the parameters of the system to their values. For example:

```julia
u0 = [x => 1.0
      y => 0.0
      z => 0.0]

p  = [σ => 10.0
      ρ => 28.0
      β => 8/3]
tspan = (0.0,100.0)
prob = ODEProblem(sys,u0,tspan,p;jac=true)
```

Note that the additional `jac=true` tells the system to symbolically generate
an optimized Jacobian function to enhance the differential equation solvers.

### Example 2: Building a Component-Based ODEProblem with Sparse Jacobians

### Example 3: Building Nonlinear Systems to Solve with NLsolve.jl

We can also build nonlinear systems. Let's say we wanted to solve for the steady
state of the previous ODE. This is the nonlinear system defined by where the
derivatives are zero. We use (unknown) variables for our nonlinear system.

```julia
using ModelingToolkit

@variables x y z
@parameters σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs, [x,y,z], [σ,ρ,β])
nlsys_func = generate_function(ns)[2] # second is the inplace version
```

which generates:

```julia
(var"##MTIIPVar#405", u, p)->begin
        @inbounds begin
                @inbounds begin
                        let (x, y, z, σ, ρ, β) = (u[1], u[2], u[3], p[1], p[2], p[3])
                            var"##MTIIPVar#405"[1] = (*)(σ, (-)(y, x))
                            var"##MTIIPVar#405"[2] = (-)((*)(x, (-)(ρ, z)), y)
                            var"##MTIIPVar#405"[3] = (-)((*)(x, y), (*)(β, z))
                        end
                    end
            end
        nothing
    end
```

We can use this to build a nonlinear function for use with NLsolve.jl:

```julia
f = eval(nlsys_func)
du = zeros(3); u = ones(3)
f(du,u,(10.0,26.0,2.33))
du

#=
3-element Array{Float64,1}:
  0.0
 24.0
 -1.33
 =#
```

We can similarly ask to generate the in-place Jacobian function:

```julia
j_func = generate_jacobian(ns)[2] # second is in-place
j! = eval(j_func)
```

which gives:

```julia
:((var"##MTIIPVar#582", u, p)->begin
          #= C:\Users\accou\.julia\dev\ModelingToolkit\src\utils.jl:70 =#
          #= C:\Users\accou\.julia\dev\ModelingToolkit\src\utils.jl:71 =#
          #= C:\Users\accou\.julia\dev\ModelingToolkit\src\utils.jl:71 =# @inbounds begin
                  #= C:\Users\accou\.julia\dev\ModelingToolkit\src\utils.jl:72 =#
                  #= C:\Users\accou\.julia\dev\ModelingToolkit\src\utils.jl:53 =# @inbounds begin
                          #= C:\Users\accou\.julia\dev\ModelingToolkit\src\utils.jl:53 =#
                          let (x, y, z, σ, ρ, β) = (u[1], u[2], u[3], p[1], p[2], p[3])
                              var"##MTIIPVar#582"[1] = (*)(σ, -1)
                              var"##MTIIPVar#582"[2] = (-)(ρ, z)
                              var"##MTIIPVar#582"[3] = y
                              var"##MTIIPVar#582"[4] = σ
                              var"##MTIIPVar#582"[5] = -1
                              var"##MTIIPVar#582"[6] = x
                              var"##MTIIPVar#582"[7] = 0
                              var"##MTIIPVar#582"[8] = (*)(x, -1)
                              var"##MTIIPVar#582"[9] = (*)(-1, β)
                          end
                      end
              end
          #= C:\Users\accou\.julia\dev\ModelingToolkit\src\utils.jl:74 =#
          nothing
      end)
```

Now we can call `nlsolve` by enclosing our parameters into the functions:

```julia
nlsolve((out, x) -> f(out, x, params), (out, x) -> j!(out, x, params), ones(3))
```

If one would like the generated function to be a Julia function instead of an expression, and allow this
function to be used from within the same world-age, one simply needs to pass `Val{false}` to tell it to
generate the function, i.e.:

```julia
nlsys_func = generate_function(ns, [x,y,z], [σ,ρ,β], Val{false})[2]
```

which uses GeneralizedGenerated.jl to build a same world-age function on the fly without eval.

## High Level API Documentation

```@docs
@parameters
@variables
@derivatives
Base.:~(::Expression, ::Expression)
```

## Additional High Level Explanations and Tips

## The Auto-Detecting System Constructors

For the high level interface, the system constructors such as `ODESystem` have
high level constructors which just take in the required equations and automatically
parse the expressions to figure out the states and parameters of the system.
The following high level constructors exist:

```julia
ODESystem(eqs)
NonlinearSystem(eqs)
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
ns = NonlinearSystem(eqs, [x,y,z])
nlsys_func = generate_function(ns, [x,y,z], [σ,ρ,β])
```

expands to:

```julia
:((##365, u, p)->begin
          let (x, y, z, σ, ρ, β) = (u[1], u[2], u[3], p[1], p[2], p[3])
              ##365[1] = σ * (y - x)
              ##365[2] = x * (ρ - z) - y
              ##365[3] = x * y - β * z
          end
      end)
```

In addition, the Jacobian calculations take into account intermediate variables
to appropriately handle them.
