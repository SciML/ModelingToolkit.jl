# High Level API

The high-level API allows modelers to interactively build models in a symbolic
manner. It is designed as a semi-DSL for easily building large complex models
and manipulating the models to generate optimal forms to be used in numerical
methods.

## Examples

### Example 1: Symbolically Building an ODEProblem for DifferentialEquations.jl

Let's build an ODE. First, we define some variables. In a differential equation
system, we need to differentiate between our (dependent) variables
and parameters. Therefore, we label them as follows:

```julia
using ModelingToolkit

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
sys = ODESystem(eqs)
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
prob = ODEProblem(sys,u0,tspan,p;jac=true,sparse=true)
```

Note that the additional `jac=true` tells the system to symbolically generate
an optimized Jacobian function to enhance the differential equation solvers,
and `sparse` tells it to build the ODEProblem with all of the enhancements
setup for sparse Jacobians.

### Example 2: Building a Component-Based ODEProblem

In addition, we can then use ModelingToolkit to compose multiple ODE subsystems.
Let's define two interacting Lorenz equations:

```julia
lorenz1 = ODESystem(eqs,name=:lorenz1)
lorenz2 = ODESystem(eqs,name=:lorenz2)

@variables α
@parameters γ
connections = [0 ~ lorenz1.x + lorenz2.y + sin(α*γ)]
connected = ODESystem(connections,[α],[γ],systems=[lorenz1,lorenz2])
```

which is now a differential-algebraic equation (DAE) of 7 variables, which has
two independent Lorenz systems and an algebraic equation that determines `α`
such that an implicit constraint holds. We can then define the resulting
`ODEProblem` and send it over to DifferentialEquations.jl.

### Example 3: Building Nonlinear Systems to Solve with NLsolve.jl

In this example we will go one step deeper and showcase the direct function
generation capabilities in ModelingToolkit.jl to build nonlinear systems.
Let's say we wanted to solve for the steady state of the previous ODE. This is
the nonlinear system defined by where the derivatives are zero. We use (unknown)
variables for our nonlinear system.

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
params = (10.0,26.0,2.33)
f(du,u,params)
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

Now, we can call `nlsolve` by enclosing our parameters into the functions:

```julia
using NLsolve
nlsolve((out, x) -> f(out, x, params), (out, x) -> j!(out, x, params), ones(3))
```

If one would like the generated function to be a Julia function instead of an expression, and allow this
function to be used from within the same world-age, one simply needs to pass `Val{false}` to tell it to
generate the function, i.e.:

```julia
nlsys_func = generate_function(ns, [x,y,z], [σ,ρ,β], expression=Val{false})[2]
```

which uses GeneralizedGenerated.jl to build the same world-age function on the fly without eval.

## High-Level API Documentation

```@docs
@parameters
@variables
@derivatives
Base.:~(::Expression, ::Expression)
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

Because the ModelingToolkit `Expression` types obey Julia semantics, one can
directly transform existing Julia functions into ModelingToolkit symbolic
representations of the function by simply inputting the symbolic values into
the function and using what is returned. For example, let's take the following
numerical PDE discretization:

```julia
using ModelingToolkit, LinearAlgebra, SparseArrays

# Define the constants for the PDE
const α₂ = 1.0
const α₃ = 1.0
const β₁ = 1.0
const β₂ = 1.0
const β₃ = 1.0
const r₁ = 1.0
const r₂ = 1.0
const _DD = 100.0
const γ₁ = 0.1
const γ₂ = 0.1
const γ₃ = 0.1
const N = 8
const X = reshape([i for i in 1:N for j in 1:N],N,N)
const Y = reshape([j for i in 1:N for j in 1:N],N,N)
const α₁ = 1.0.*(X.>=4*N/5)

const Mx = Tridiagonal([1.0 for i in 1:N-1],[-2.0 for i in 1:N],[1.0 for i in 1:N-1])
const My = copy(Mx)
Mx[2,1] = 2.0
Mx[end-1,end] = 2.0
My[1,2] = 2.0
My[end,end-1] = 2.0

# Define the discretized PDE as an ODE function
function f!(du,u,p,t)
   A = @view  u[:,:,1]
   B = @view  u[:,:,2]
   C = @view  u[:,:,3]
  dA = @view du[:,:,1]
  dB = @view du[:,:,2]
  dC = @view du[:,:,3]
  mul!(MyA,My,A)
  mul!(AMx,A,Mx)
  @. DA = _DD*(MyA + AMx)
  @. dA = DA + α₁ - β₁*A - r₁*A*B + r₂*C
  @. dB = α₂ - β₂*B - r₁*A*B + r₂*C
  @. dC = α₃ - β₃*C + r₁*A*B - r₂*C
end
```

We can then define the corresponding arrays as ModelingToolkit variables:

```julia
# Define the initial condition as normal arrays
@variables du[1:N,1:N,1:3] u[1:N,1:N,1:3] MyA[1:N,1:N] AMx[1:N,1:N] DA[1:N,1:N]
f!(du,u,nothing,0.0)
```

The output, here the in-place modified `du`, is a symbolic representation of
each output of the function. We can then utilize this in the ModelingToolkit
functionality. For example, let's compute the sparse Jacobian function and
compile a fast multithreaded version:

```julia
jac = sparse(ModelingToolkit.jacobian(vec(du),vec(u),simplify=false))
multithreadedjac = eval(ModelingToolkit.build_function(vec(jac),u,multithread=true)[2])
```

### modelingtoolkitize

For some `DEProblem` types, automatic tracing functionality is already included
via the `modelingtoolkitize` function. Take, for example, the Robertson ODE
defined as an `ODEProblem` for DifferentialEquations.jl:

```julia
using DifferentialEquations
function rober(du,u,p,t)
  y₁,y₂,y₃ = u
  k₁,k₂,k₃ = p
  du[1] = -k₁*y₁+k₃*y₂*y₃
  du[2] =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  du[3] =  k₂*y₂^2
  nothing
end
prob = ODEProblem(rober,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))
```

If we want to get a symbolic representation, we can simply call `modelingtoolkitize`
on the `prob`, which will return an `ODESystem`:

```julia
sys = modelingtoolkitize(prob)
```

Using this, we can symbolically build the Jacobian and then rebuild the ODEProblem:

```julia
jac = eval(ModelingToolkit.generate_jacobian(sys)[2])

f = ODEFunction(rober, jac=jac)
prob_jac = ODEProblem(f,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))
```

```@docs
modelingtoolkitize
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
