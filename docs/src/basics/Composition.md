# [Composing Models and Building Reusable Components](@id components)

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
of `decay2` is the value of the unknown variable `x`.

```@example composition
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

function decay(; name)
    @parameters a
    @variables x(t) f(t)
    ODESystem([
            D(x) ~ -a * x + f
        ], t;
        name = name)
end

@named decay1 = decay()
@named decay2 = decay()

connected = compose(
    ODESystem([decay2.f ~ decay1.x
               D(decay1.f) ~ 0], t; name = :connected), decay1, decay2)

equations(connected)

#4-element Vector{Equation}:
# Differential(t)(decay1₊f(t)) ~ 0
# decay2₊f(t) ~ decay1₊x(t)
# Differential(t)(decay1₊x(t)) ~ decay1₊f(t) - (decay1₊a*(decay1₊x(t)))
# Differential(t)(decay2₊x(t)) ~ decay2₊f(t) - (decay2₊a*(decay2₊x(t)))

simplified_sys = structural_simplify(connected)

equations(simplified_sys)
```

Now we can solve the system:

```@example composition
x0 = [decay1.x => 1.0
      decay1.f => 0.0
      decay2.x => 1.0]
p = [decay1.a => 0.1
     decay2.a => 0.2]

using OrdinaryDiffEq
prob = ODEProblem(simplified_sys, x0, (0.0, 100.0), p)
sol = solve(prob, Tsit5())
sol[decay2.f]
```

## Basics of Model Composition

Every `AbstractSystem` has a `system` keyword argument for specifying
subsystems. A model is the composition of itself and its subsystems.
For example, if we have:

```julia
@named sys = compose(ODESystem(eqs, indepvar, unknowns, ps), subsys)
```

the `equations` of `sys` is the concatenation of `get_eqs(sys)` and
`equations(subsys)`, the unknowns are the concatenation of their unknowns,
etc. When the `ODEProblem` or `ODEFunction` is generated from this
system, it will build and compile the functions associated with this
composition.

The new equations within the higher level system can access the variables
in the lower level system by namespacing via the `nameof(subsys)`. For
example, let's say there is a variable `x` in `unknowns` and a variable
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
u0 = [x => 2.0
      subsys.x => 2.0]
```

Note that any default values within the given subcomponent will be
used if no override is provided at construction time. If any values for
initial conditions or parameters are unspecified, an error will be thrown.

When the model is numerically solved, the solution can be accessed via
its symbolic values. For example, if `sol` is the `ODESolution`, one
can use `sol[x]` and `sol[subsys.x]` to access the respective timeseries
in the solution. All other indexing rules stay the same, so `sol[x,1:5]`
accesses the first through fifth values of `x`. Note that this can be
done even if the variable `x` is eliminated from the system from
transformations like `alias_elimination` or `tearing`: the variable
will be lazily reconstructed on demand.

### Variable scope and parameter expressions

In some scenarios, it could be useful for model parameters to be expressed
in terms of other parameters, or shared between common subsystems.
To facilitate this, ModelingToolkit supports symbolic expressions
in default values, and scoped variables.

With symbolic parameters, it is possible to set the default value of a parameter or initial condition to an expression of other variables.

```julia
# ...
sys = ODESystem(
# ...
# directly in the defaults argument
    defaults = Pair{Num, Any}[x => u,
    y => σ,
    z => u - 0.1])
# by assigning to the parameter
sys.y = u * 1.1
```

In a hierarchical system, variables of the subsystem get namespaced by the name of the system they are in. This prevents naming clashes, but also enforces that every unknown and parameter is local to the subsystem it is used in. In some cases it might be desirable to have variables and parameters that are shared between subsystems, or even global. This can be accomplished as follows.

```julia
@parameters a b c d e f

# a is a local variable
b = ParentScope(b) # b is a variable that belongs to one level up in the hierarchy
c = ParentScope(ParentScope(c)) # ParentScope can be nested
d = DelayParentScope(d) # skips one level before applying ParentScope
e = DelayParentScope(e, 2) # second argument allows skipping N levels
f = GlobalScope(f)

p = [a, b, c, d, e, f]

level0 = ODESystem(Equation[], t, [], p; name = :level0)
level1 = ODESystem(Equation[], t, [], []; name = :level1) ∘ level0
parameters(level1)
#level0₊a
#b
#c
#level0₊d
#level0₊e
#f
level2 = ODESystem(Equation[], t, [], []; name = :level2) ∘ level1
parameters(level2)
#level1₊level0₊a
#level1₊b
#c
#level0₊d
#level1₊level0₊e
#f
level3 = ODESystem(Equation[], t, [], []; name = :level3) ∘ level2
parameters(level3)
#level2₊level1₊level0₊a
#level2₊level1₊b
#level2₊c
#level2₊level0₊d
#level1₊level0₊e
#f
```

## Structural Simplify

In many cases, the nicest way to build a model may leave a lot of
unnecessary variables. Thus one may want to remove these equations
before numerically solving. The `structural_simplify` function removes
these trivial equality relationships and trivial singularity equations,
i.e. equations which result in `0~0` expressions, in over-specified systems.

## Inheritance and Combine

Model inheritance can be done in two ways: implicitly or explicitly. First, one
can use the `extend` function to extend a base model with another set of
equations, unknowns, and parameters. An example can be found in the
[acausal components tutorial](@ref acausal).

The explicit way is to shadow variables with equality expressions. For example,
let's assume we have three separate systems which we want to compose to a single
one. This is how one could explicitly forward all unknowns and parameters to the
higher level system:

```@example compose
using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D

## Library code
@variables S(t), I(t), R(t)
N = S + I + R
@parameters β, γ

@named seqn = ODESystem([D(S) ~ -β * S * I / N], t)
@named ieqn = ODESystem([D(I) ~ β * S * I / N - γ * I], t)
@named reqn = ODESystem([D(R) ~ γ * I], t)

sir = compose(
    ODESystem(
        [
            S ~ ieqn.S,
            I ~ seqn.I,
            R ~ ieqn.R,
            ieqn.S ~ seqn.S,
            seqn.I ~ ieqn.I,
            seqn.R ~ reqn.R,
            ieqn.R ~ reqn.R,
            reqn.I ~ ieqn.I],
        t,
        [S, I, R],
        [β, γ],
        defaults = [seqn.β => β
                    ieqn.β => β
                    ieqn.γ => γ
                    reqn.γ => γ], name = :sir),
    seqn,
    ieqn,
    reqn)
```

Note that the unknowns are forwarded by an equality relationship, while
the parameters are forwarded through a relationship in their default
values. The user of this model can then solve this model simply by
specifying the values at the highest level:

```@example compose
sireqn_simple = structural_simplify(sir)

equations(sireqn_simple)
```

```@example compose
## User Code

u0 = [seqn.S => 990.0,
    ieqn.I => 10.0,
    reqn.R => 0.0]

p = [β => 0.5
     γ => 0.25]

tspan = (0.0, 40.0)
prob = ODEProblem(sireqn_simple, u0, tspan, p, jac = true)
sol = solve(prob, Tsit5())
sol[reqn.R]
```

## Tearing Problem Construction

Some system types (specifically `NonlinearSystem`) can be further
reduced if `structural_simplify` has already been applied to them. This is done
by using the alternative problem constructors (`BlockNonlinearProblem`).
In these cases, the constructor uses the knowledge of the
strongly connected components calculated during the process of simplification
as the basis for building pre-simplified nonlinear systems in the implicit
solving. In summary: these problems are structurally modified, but could be
more efficient and more stable.

## Components with discontinuous dynamics

When modeling, e.g., impacts, saturations or Coulomb friction, the dynamic
equations are discontinuous in either the unknown or one of its derivatives. This
causes the solver to take very small steps around the discontinuity, and
sometimes leads to early stopping due to `dt <= dt_min`. The correct way to
handle such dynamics is to tell the solver about the discontinuity by a
root-finding equation, which can be modeling using [`ODESystem`](@ref)'s event
support. Please see the tutorial on [Callbacks and Events](@ref events) for
details and examples.
