# Parameter Identifiability in ODE Models

Using ordinary differential equations in modeling processes is commonplace and the challenge of parameter identifiability is one of the key design challenges. In this tutorial, we will show how to use `StructuralIdentifiability.jl` with `ModelingToolkit.jl` to assess parameter identifiability.

We will start with determining local identifiability, where a parameter is known up to finitely many values, and then proceed to determining global identifiability properties, that is, which parameters can be identified uniquely.

## Local Identifiability
### Input System

We will consider a simple two-species competition model

$$\begin{cases}
\frac{d\,x_4}{d\,t} = - \frac{k_5 x_4}{k_6 + x_4},\\
\frac{d\,x_5}{d\,t} = \frac{k_5 x_4}{k_6 + x_4} - \frac{k_7 x_5}{(k_8 + x_5 + x_6)},\\
\frac{d\,x_6}{d\,t} = \frac{k_7 x_5}{(k_8 + x_5 + x_6)} - \frac{k_9  x_6  (k_{10} - x_6) }{k_{10}},\\
\frac{d\,x_7}{d\,t} = \frac{k_9  x_6  (k_{10} - x_6)}{ k_{10}},\\
y_1 = x_4,\\
y_2 = x_5\end{cases}$$

This model describes the biohydrogenation[^1] process[^2] with unknown initial conditions.

### Using the `ODESystem` object
To define the system in Julia, we use `ModelingToolkit.jl`.

We first define the parameters, variables, differential equations and the output equations. Notice that the system does not have any input functions, so inputs will be an empty array.

```@example
using StructuralIdentifiability, ModelingToolkit

# define parameters and variables
@variables t x4(t) x5(t) x6(t) x7(t) y1(t) y2(t)
@parameters k5 k6 k7 k8 k9 k10
D = Differential(t)

# define equations
eqs = [
    D(x4) ~ - k5 * x4 / (k6 + x4),
    D(x5) ~ k5 * x4 / (k6 + x4) - k7 * x5/(k8 + x5 + x6),
    D(x6) ~ k7 * x5 / (k8 + x5 + x6) - k9 * x6 * (k10 - x6) / k10,
    D(x7) ~ k9 * x6 * (k10 - x6) / k10
]

# define observed functions
observed = [
    y1 ~ x4,
    y2 ~ x5
]

# define the system
de = ODESystem(eqs, t, [x4, x5, x6, x7], [k5, k6, k7, k8, k9, k10], observed=observed, name=:Biohydrogenation)

# no input functions:
inputs = []

# we want to check everything
to_check = []

# query local identifiability
# we pass the ode-system
local_id_all = assess_local_identifiability(de, inputs, to_check, 0.99)

# let's try to check specific parameters and their combinations
to_check = [k5, k7, k10/k9, k5+k6]
local_id_some = assess_local_identifiability(de, inputs, to_check, 0.99)

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

Global identifiability needs information about local identifiability first, hence the function we chose here will take care of that extra step for us.

```@repl
using StructuralIdentifiability, ModelingToolkit
@parameters b c α β γ δ σ
@variables t x1(t) x2(t) x3(t) x4(t) y(t)
D = Differential(t)

eqs = [
    D(x1) ~ -b * x1 + 1/(c + x4),
    D(x2) ~ α * x1 - β * x2,
    D(x3) ~ γ * x2 - δ * x3,
    D(x4) ~ σ * x4 * (γ * x2 - δ * x3)/x3
]

observed = [
    y~x1
]

# no inputs
inputs = []

# check all parameters
to_check = []

ode = ODESystem(eqs, t, [x1, x2, x3, x4], [b, c, α, β, γ, δ, σ], observed=observed, name=:GoodwinOsc)

@time global_id = assess_identifiability(ode, inputs, to_check, 0.99)
```

Let us consider the same system but with two inputs and we will try to find out identifiability with probability `0.9` for parameters `c` and `b`:

```@repl
using StructuralIdentifiability, ModelingToolkit
@parameters b c α β γ δ σ
@variables t x1(t) x2(t) x3(t) x4(t) y(t) u1 u2
D = Differential(t)

eqs = [
    D(x1) ~ -b * x1 + 1/(c + x4),
    D(x2) ~ α * x1 - β * x2 - u1,
    D(x3) ~ γ * x2 - δ * x3 + u2,
    D(x4) ~ σ * x4 * (γ * x2 - δ * x3)/x3
]

observed = [
    y~x1
]

# no inputs
inputs = [u1, u2]

# check all parameters
to_check = [b, c]

ode = ODESystem(eqs, t, [x1, x2, x3, x4], [b, c, α, β, γ, δ, σ], observed=observed, name=:GoodwinOsc)

@time global_id = assess_identifiability(ode, inputs, to_check, 0.9)
```


[^1]:
    > R. Munoz-Tamayo, L. Puillet, J.B. Daniel, D. Sauvant, O. Martin, M. Taghipoor, P. Blavy [*Review: To be or not to be an identifiable model. Is this a relevant question in animal science modelling?*](https://doi.org/10.1017/S1751731117002774), Animal, Vol 12 (4), 701-712, 2018. The model is the ODE system (3) in Supplementary Material 2, initial conditions are assumed to be unknown.

[^2]:
    > Moate P.J., Boston R.C., Jenkins T.C. and Lean I.J., [*Kinetics of Ruminal Lipolysis of Triacylglycerol and Biohydrogenationof Long-Chain Fatty Acids: New Insights from Old Data*](doi:10.3168/jds.2007-0398), Journal of Dairy Science 91, 731–742, 2008

[^3]:
    > Goodwin, B.C. [*Oscillatory behavior in enzymatic control processes*](https://doi.org/10.1016/0065-2571(65)90067-1), Advances in Enzyme Regulation, Vol 3 (C), 425-437, 1965