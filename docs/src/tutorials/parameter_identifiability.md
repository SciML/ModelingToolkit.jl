# Parameter Identifiability in ODE Models

Using ordinary differential equations for modeling real-world processes is commonplace and the problem of parameter identifiability is one of the key design challenges. We say that a parameter is identifiable if we can recover its value from experimental data. When we describe the parameter without actual data at hand is _structurally_ identifiable. In this tutorial, we will show how to use `StructuralIdentifiability.jl` with `ModelingToolkit.jl` to assess parameter identifiability.

We will start with determining local identifiability, where a parameter is known up to finitely many values, and then proceed to determining global identifiability properties, that is, which parameters can be identified uniquely.

To install `StructuralIdentifiability.jl`, simply run
```julia
using Pkg
Pkg.add("StructuralIdentifiability")
```

The package has a standalone data structure for ordinary differential equations but is also compatible with `ODESystem` type from `ModelingToolkit.jl`.

Let's start with local identifiability!
## Local Identifiability
### Input System

We will consider the following model:

$$\begin{cases}
\frac{d\,x_4}{d\,t} = - \frac{k_5 x_4}{k_6 + x_4},\\
\frac{d\,x_5}{d\,t} = \frac{k_5 x_4}{k_6 + x_4} - \frac{k_7 x_5}{(k_8 + x_5 + x_6)},\\
\frac{d\,x_6}{d\,t} = \frac{k_7 x_5}{(k_8 + x_5 + x_6)} - \frac{k_9  x_6  (k_{10} - x_6) }{k_{10}},\\
\frac{d\,x_7}{d\,t} = \frac{k_9  x_6  (k_{10} - x_6)}{ k_{10}},\\
y_1 = x_4,\\
y_2 = x_5\end{cases}$$

This model describes the biohydrogenation[^1] process[^2] with unknown initial conditions.

### Using the `ODESystem` object
To define the ode system in Julia, we use `ModelingToolkit.jl`.

We first define the parameters, variables, differential equations and the output equations.
```julia
using StructuralIdentifiability, ModelingToolkit

# define parameters and variables
@variables t x4(t) x5(t) x6(t) x7(t) y1(t) [output=true] y2(t) [output=true]
@parameters k5 k6 k7 k8 k9 k10
D = Differential(t)

# define equations
eqs = [
    D(x4) ~ - k5 * x4 / (k6 + x4),
    D(x5) ~ k5 * x4 / (k6 + x4) - k7 * x5/(k8 + x5 + x6),
    D(x6) ~ k7 * x5 / (k8 + x5 + x6) - k9 * x6 * (k10 - x6) / k10,
    D(x7) ~ k9 * x6 * (k10 - x6) / k10,
    y1 ~ x4,
    y2 ~ x5
]

# define the system
de = ODESystem(eqs, t, name=:Biohydrogenation)

```

After that we are ready to check the system for local identifiability:
```julia
# query local identifiability
# we pass the ode-system
local_id_all = assess_local_identifiability(de, 0.99)
                # [ Info: Preproccessing `ModelingToolkit.ODESystem` object
                # 6-element Vector{Bool}:
                #  1
                #  1
                #  1
                #  1
                #  1
                #  1
```
We can see that all states (except $x_7$) and all parameters are locally identifiable with probability 0.99. 

Let's try to check specific parameters and their combinations
```julia
to_check = [k5, k7, k10/k9, k5+k6]
local_id_some = assess_local_identifiability(de, to_check, 0.99)
                # 4-element Vector{Bool}:
                #  1
                #  1
                #  1
                #  1
```

Notice that in this case, everything (except the state variable $x_7$) is locally identifiable, including combinations such as $k_{10}/k_9, k_5+k_6$

## Global Identifiability

In this tutorial, let us cover an example problem of querying the ODE for globally identifiable parameters.

### Input System

Let us consider the following four-dimensional model with two outputs:

$$\begin{cases}
    x_1'(t) = -b  x_1(t) + \frac{1 }{ c + x_4(t)},\\
    x_2'(t) = \alpha  x_1(t) - \beta  x_2(t),\\
    x_3'(t) = \gamma  x_2(t) - \delta  x_3(t),\\
    x_4'(t) = \sigma  x_4(t)  \frac{(\gamma x_2(t) - \delta x_3(t))}{ x_3(t)},\\
    y(t) = x_1(t)
\end{cases}$$

This model describes enzyme dynamics[^3]. Let us run a global identifiability check on this model. We will use the default settings: the probability of correctness will be `p=0.99` and we are interested in identifiability of all possible parameters

Global identifiability needs information about local identifiability first, but the function we chose here will take care of that extra step for us.

__Note__: as of writing this tutorial, UTF-symbols such as Greek characters are not supported by one of the project's dependencies, see (this issue)[https://github.com/SciML/StructuralIdentifiability.jl/issues/43].

```julia
using StructuralIdentifiability, ModelingToolkit
@parameters b c a beta g delta sigma
@variables t x1(t) x2(t) x3(t) x4(t) y(t) [output=true]
D = Differential(t)

eqs = [
    D(x1) ~ -b * x1 + 1/(c + x4),
    D(x2) ~ a * x1 - beta * x2,
    D(x3) ~ g * x2 - delta * x3,
    D(x4) ~ sigma * x4 * (g * x2 - delta * x3)/x3,
    y~x1
]


ode = ODESystem(eqs, t, name=:GoodwinOsc)

@time global_id = assess_identifiability(ode)
                # 28.961573 seconds (88.92 M allocations: 5.541 GiB, 4.01% gc time)
                # Dict{Nemo.fmpq_mpoly, Symbol} with 7 entries:
                #   c     => :globally
                #   a     => :nonidentifiable
                #   g     => :nonidentifiable
                #   delta => :locally
                #   sigma => :globally
                #   beta  => :locally
                #   b     => :globally
```
We can see that 

Let us consider the same system but with two inputs and we will try to find out identifiability with probability `0.9` for parameters `c` and `b`:

```julia
using StructuralIdentifiability, ModelingToolkit
@parameters b c a beta g delta sigma
@variables t x1(t) x2(t) x3(t) x4(t) y(t) [output=true] u1(t) [input=true] u2(t) [input=true]
D = Differential(t)

eqs = [
    D(x1) ~ -b * x1 + 1/(c + x4),
    D(x2) ~ a * x1 - beta * x2 - u1,
    D(x3) ~ g * x2 - delta * x3 + u2,
    D(x4) ~ sigma * x4 * (g * x2 - delta * x3)/x3,
    y~x1
]

# check only 2 parameters
to_check = [b, c]

ode = ODESystem(eqs, t, name=:GoodwinOsc)

global_id = assess_identifiability(ode, to_check, 0.9)
            # Dict{Num, Symbol} with 2 entries:
            #   b => :globally
            #   c => :globally
```

Both parameters $b, c$ are globally identifiable with probability 0.9.

[^1]:
    > R. Munoz-Tamayo, L. Puillet, J.B. Daniel, D. Sauvant, O. Martin, M. Taghipoor, P. Blavy [*Review: To be or not to be an identifiable model. Is this a relevant question in animal science modelling?*](https://doi.org/10.1017/S1751731117002774), Animal, Vol 12 (4), 701-712, 2018. The model is the ODE system (3) in Supplementary Material 2, initial conditions are assumed to be unknown.

[^2]:
    > Moate P.J., Boston R.C., Jenkins T.C. and Lean I.J., [*Kinetics of Ruminal Lipolysis of Triacylglycerol and Biohydrogenationof Long-Chain Fatty Acids: New Insights from Old Data*](doi:10.3168/jds.2007-0398), Journal of Dairy Science 91, 731–742, 2008

[^3]:
    > Goodwin, B.C. [*Oscillatory behavior in enzymatic control processes*](https://doi.org/10.1016/0065-2571(65)90067-1), Advances in Enzyme Regulation, Vol 3 (C), 425-437, 1965