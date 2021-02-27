# [Composing Models and Building Reusable Components](@ref components)

The symbolic models of ModelingToolkit can be composed together to
easily build large models. The composition is lazy and only instantiated
at the time of conversion to numerical models, allowing a more performant
way in terms of computation time and memory.

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
like `alias_elimination` or `tearing` which will be described later.

### Simple Model Composition Example

The following is an example of building a model in a library with
an optional forcing function, and allowing the user to specify the
forcing later. Here, the library author defines a component named
`decay`. The user then builds two `decay` components and connects them,
saying the forcing term of `decay1` is a constant while the forcing term
of `decay2` is the value of the state variable `x`.

```julia
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
connected = ODESystem([
                        decay2.f ~ decay1.x
                        D(decay1.f) ~ 0
                      ], t, systems=[decay1, decay2])
```

## Combine

## Alias Elimination

## The Tearing Transformation

## Automatic Model Promotion
