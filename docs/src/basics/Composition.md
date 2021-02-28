# [Composing Models and Building Reusable Components](@ref components)

The symbolic models of ModelingToolkit can be composed together to
easily build large models. The composition is lazy and only instantiated
at the time of conversion to numerical models, allowing a more performant
way in terms of computation time and memory.

## Simple Model Composition Example

The following is an example of building a model in a library with
an optional forcing function, and allowing the user to specify the
forcing later. Here, the library author defines a component named
`decay`. The user then builds two `decay` components and connects them,
saying the forcing term of `decay1` is a constant while the forcing term
of `decay2` is the value of the state variable `x`.

```julia
using ModelingToolkit

function decay(;name)
  @parameters t a
  @variables x(t) f(t)
  D = Differential(t)
  ODESystem([
      D(x) ~ -a*x + f
    ];
    name=name)
end

@named decay1 = decay()
@named decay2 = decay()

@parameters t
D = Differential(t)
@named connected = ODESystem([
                        decay2.f ~ decay1.x
                        D(decay1.f) ~ 0
                      ], t, systems=[decay1, decay2])

equations(connected)

#4-element Vector{Equation}:
# Differential(t)(decay1₊f(t)) ~ 0
# decay2₊f(t) ~ decay1₊x(t)
# Differential(t)(decay1₊x(t)) ~ decay1₊f(t) - (decay1₊a*(decay1₊x(t)))
# Differential(t)(decay2₊x(t)) ~ decay2₊f(t) - (decay2₊a*(decay2₊x(t)))

simplified_sys = structural_simplify(connected)

equations(simplified_sys)

#3-element Vector{Equation}:
# Differential(t)(decay1₊f(t)) ~ 0
# Differential(t)(decay1₊x(t)) ~ decay1₊f(t) - (decay1₊a*(decay1₊x(t)))
# Differential(t)(decay2₊x(t)) ~ decay1₊x(t) - (decay2₊a*(decay2₊x(t)))
```

Now we can solve the system:

```julia
x0 = [
  decay1.x => 1.0
  decay1.f => 0.0
  decay2.x => 1.0
]
p = [
  decay1.a => 0.1
  decay2.a => 0.2
]

using DifferentialEquations
prob = ODEProblem(simplified_sys, x0, (0.0, 100.0), p)
sol = solve(prob, Tsit5())
sol[decay2.f]
```

## Basics of Model Composition

Every `AbstractSystem` has a `system` keyword argument for specifying
subsystems. A model is the composition of itself and its subsystems.
For example, if we have:

```julia
@named sys = ODESystem(eqs,indepvar,states,ps,system=[subsys])
```

the `equations` of `sys` is the concatenation of `get_eqs(sys)` and
`equations(subsys)`, the states are the concatenation of their states,
etc. When the `ODEProblem` or `ODEFunction` is generated from this
system, it will build and compile the functions associated with this
composition.

The new equations within the higher level system can access the variables
in the lower level system by namespacing via the `nameof(subsys)`. For
example, let's say there is a variable `x` in `states` and a variable
`x` in `subsys`. We can declare that these two variables are the same
by specifying their equality: `x ~ subsys.x` in the `eqs` for `sys`.
This algebraic relationship can then be simplified by transformations
like `structural_simplify` which will be described later.

### Numerics with Composed Models

These composed models can then be directly transformed into their
associated `SciMLProblem` type using the standard constructors. When
this is done, the initial conditions and parameters must be specified
in their namespaced form. For example:

```julia
u0 = [
  x => 2.0
  subsys.x => 2.0
]
```

Note that any default values within the given subcomponent will be
used if no override is provided at construction time. If any values for
initial conditions or parameters are unspecified an error will be thrown.

When the model is numerically solved, the solution can be accessed via
its symbolic values. For example, if `sol` is the `ODESolution`, one
can use `sol[x]` and `sol[subsys.x]` to access the respective timeseries
in the solution. All other indexing rules stay the same, so `sol[x,1:5]`
accesses the first through fifth values of `x`. Note that this can be
done even if the variable `x` is eliminated from the system from
transformations like `alias_elimination` or `tearing`: the variable
will be lazily reconstructed on demand.

## Structural Simplify

In many cases, the nicest way to build a model may leave a lot of
unnecessary variables. Thus one may want to remove these equations
before numerically solving. The `structural_simplify` function removes
these trivial equality relationships and trivial singularity equations,
i.e. equations which result in `0~0` expressions, in over-specified systems.

## Inheritance and Combine (TODO)

Model inheritance can be done in two ways: explicitly or implicitly.
The explicit way is to shadow variables with equality expressions.
For example, let's assume we have three separate systems which we
want to compose to a single one. This is how one could explicitly
forward all states and parameters to the higher level system:

```julia
using ModelingToolkit, OrdinaryDiffEq, Plots

## Library code

@parameters t
D = Differential(t)

@variables S(t), I(t), R(t)
N = S + I + R
@parameters β,γ

@named seqn = ODESystem([D(S) ~ -β*S*I/N])
@named ieqn = ODESystem([D(I) ~ β*S*I/N-γ*I])
@named reqn = ODESystem([D(R) ~ γ*I])

@named sir = ODESystem([   
                    S ~ ieqn.S,
                    I ~ seqn.I,
                    R ~ ieqn.R,
                    ieqn.S ~ seqn.S,
                    seqn.I ~ ieqn.I,
                    seqn.R ~ reqn.R,
                    ieqn.R ~ reqn.R,
                    reqn.I ~ ieqn.I], t, [S,I,R], [β,γ],
                    systems=[seqn,ieqn,reqn],
                    default_p = [
                        seqn.β => β
                        ieqn.β => β
                        ieqn.γ => γ
                        reqn.γ => γ
                    ])
```

Note that the states are forwarded by an equality relationship, while
the parameters are forwarded through a relationship in their default
values. The user of this model can then solve this model simply by
specifying the values at the highest level:

```julia
sireqn_simple = structural_simplify(sir)

equations(sireqn_simple)

# 3-element Vector{Equation}:
#Differential(t)(seqn₊S(t)) ~ -seqn₊β*ieqn₊I(t)*seqn₊S(t)*(((ieqn₊I(t)) + (reqn₊R(t)) + (seqn₊S(t)))^-1)
#Differential(t)(ieqn₊I(t)) ~ ieqn₊β*ieqn₊I(t)*seqn₊S(t)*(((ieqn₊I(t)) + (reqn₊R(t)) + (seqn₊S(t)))^-1) - (ieqn₊γ*(ieqn₊I(t)))
#Differential(t)(reqn₊R(t)) ~ reqn₊γ*ieqn₊I(t)

## User Code

u0 = [seqn.S => 990.0,
      ieqn.I => 10.0,
      reqn.R => 0.0]

p = [
    β => 0.5
    γ => 0.25
]

tspan = (0.0,40.0)
prob = ODEProblem(sireqn_simple,u0,tspan,p,jac=true)
sol = solve(prob,Tsit5())
sol[reqn.R]
```

However, one can similarly simplify this process of inheritance by
using `combine` which concatenates all of the vectors within the
systems. For example, we could equivalently have done:

```julia
@named sir = combine([seqn,ieqn,reqn])
```

## Tearing Problem Construction

Some system types, specifically `ODESystem` and `NonlinearSystem`, can be further
reduced if `structural_simplify` has already been applied to them. This is done
by using the alternative problem constructors, `ODAEProblem` and `BlockNonlinearProblem`
respectively. In these cases, the constructor uses the knowledge of the
strongly connected components calculated during the process of simplification
as the basis for building pre-simplified nonlinear systems in the implicit
solving. In summary: these problems are structurally modified, but could be
more efficient and more stable.

## Automatic Model Promotion (TODO)

In many cases one might want to compose models of different types. For example,
one may want to include a `NonlinearSystem` as a set of algebraic equations
within an `ODESystem`, or one may want to use an `ODESystem` as a subsystem of
an `SDESystem`. In these cases, the compostion works automatically by promoting
the model via `promote_system`. System promotions exist in the cases where a
mathematically-trivial definition of the promotion exists.
