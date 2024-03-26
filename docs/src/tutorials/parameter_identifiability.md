# Parameter Identifiability in ODE Models

Ordinary differential equations are commonly used for modeling real-world processes. The problem of parameter identifiability is one of the key design challenges for mathematical models. A parameter is said to be _identifiable_ if one can recover its value from experimental data. _Structural_ identifiability is a theoretical property of a model that answers this question. In this tutorial, we will show how to use `StructuralIdentifiability.jl` with `ModelingToolkit.jl` to assess identifiability of parameters in ODE models. The theory behind `StructuralIdentifiability.jl` is presented in paper [^4].

We will start by illustrating **local identifiability** in which a parameter is known up to _finitely many values_, and then proceed to determining **global identifiability**, that is, which parameters can be identified _uniquely_.

The package has a standalone data structure for ordinary differential equations, but is also compatible with `ODESystem` type from `ModelingToolkit.jl`.

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
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel Biohydrogenation begin
    @variables begin
        x4(t)
        x5(t)
        x6(t)
        x7(t)
        y1(t), [output = true]
        y2(t), [output = true]
    end
    @parameters begin
        k5
        k6
        k7
        k8
        k9
        k10
    end
    # define equations
    @equations begin
        D(x4) ~ -k5 * x4 / (k6 + x4)
        D(x5) ~ k5 * x4 / (k6 + x4) - k7 * x5 / (k8 + x5 + x6)
        D(x6) ~ k7 * x5 / (k8 + x5 + x6) - k9 * x6 * (k10 - x6) / k10
        D(x7) ~ k9 * x6 * (k10 - x6) / k10
        y1 ~ x4
        y2 ~ x5
    end
end

# define the system
@mtkbuild de = Biohydrogenation()
```

After that, we are ready to check the system for local identifiability:

```julia
# query local identifiability
# we pass the ode-system
local_id_all = assess_local_identifiability(de, prob_threshold = 0.99)
```

We can see that all unknowns (except $x_7$) and all parameters are locally identifiable with probability 0.99.

Let's try to check specific parameters and their combinations

```julia
to_check = [de.k5, de.k7, de.k10 / de.k9, de.k5 + de.k6]
local_id_some = assess_local_identifiability(
    de, funcs_to_check = to_check, prob_threshold = 0.99)
```

Notice that in this case, everything (except the unknown variable $x_7$) is locally identifiable, including combinations such as $k_{10}/k_9, k_5+k_6$

## Global Identifiability

In this part tutorial, let us cover an example problem of querying the ODE for globally identifiable parameters.

### Input System

Let us consider the following four-dimensional model with two outputs:

$$\begin{cases}
x_1'(t) = -b  x_1(t) + \frac{1 }{ c + x_4(t)},\\
x_2'(t) = \alpha  x_1(t) - \beta  x_2(t),\\
x_3'(t) = \gamma  x_2(t) - \delta  x_3(t),\\
x_4'(t) = \sigma  x_4(t)  \frac{(\gamma x_2(t) - \delta x_3(t))}{ x_3(t)},\\
y(t) = x_1(t)
\end{cases}$$

We will run a global identifiability check on this enzyme dynamics[^3] model. We will use the default settings: the probability of correctness will be `p=0.99` and we are interested in identifiability of all possible parameters.

Global identifiability needs information about local identifiability first, but the function we chose here will take care of that extra step for us.

```julia
using StructuralIdentifiability, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel GoodwinOsc begin
    @parameters begin
        b
        c
        α
        β
        γ
        δ
        σ
    end
    @variables begin
        x1(t)
        x2(t)
        x3(t)
        x4(t)
        y(t), [output = true]
        y2(t), [output = true]
    end
    @equations begin
        D(x1) ~ -b * x1 + 1 / (c + x4)
        D(x2) ~ α * x1 - β * x2
        D(x3) ~ γ * x2 - δ * x3
        D(x4) ~ σ * x4 * (γ * x2 - δ * x3) / x3
        y ~ x1 + x2
        y2 ~ x2
    end
end

@named ode = GoodwinOsc()

global_id = assess_identifiability(ode)
```

We can see that only parameters `a, g` are unidentifiable, and everything else can be uniquely recovered.

Let us consider the same system but with two inputs, and we will find out identifiability with probability `0.9` for parameters `c` and `b`:

```julia
using StructuralIdentifiability, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel GoodwinOscillator begin
    @parameters begin
        b
        c
        α
        β
        γ
        δ
        σ
    end
    @variables begin
        x1(t)
        x2(t)
        x3(t)
        x4(t)
        y(t), [output = true]
        y2(t), [output = true]
        u1(t), [input = true]
        u2(t), [input = true]
    end
    @equations begin
        D(x1) ~ -b * x1 + 1 / (c + x4)
        D(x2) ~ α * x1 - β * x2 - u1
        D(x3) ~ γ * x2 - δ * x3 + u2
        D(x4) ~ σ * x4 * (γ * x2 - δ * x3) / x3
        y ~ x1 + x2
        y2 ~ x2
    end
end

@mtkbuild ode = GoodwinOscillator()

# check only 2 parameters
to_check = [ode.b, ode.c]

global_id = assess_identifiability(ode, funcs_to_check = to_check, prob_threshold = 0.9)
```

Both parameters `b, c` are globally identifiable with probability `0.9` in this case.

[^1]: > R. Munoz-Tamayo, L. Puillet, J.B. Daniel, D. Sauvant, O. Martin, M. Taghipoor, P. Blavy [*Review: To be or not to be an identifiable model. Is this a relevant question in animal science modelling?*](https://doi.org/10.1017/S1751731117002774), Animal, Vol 12 (4), 701-712, 2018. The model is the ODE system (3) in Supplementary Material 2, initial conditions are assumed to be unknown.
[^2]: > Moate P.J., Boston R.C., Jenkins T.C. and Lean I.J., [*Kinetics of Ruminal Lipolysis of Triacylglycerol and Biohydrogenationof Long-Chain Fatty Acids: New Insights from Old Data*](https://doi.org/10.3168/jds.2007-0398), Journal of Dairy Science 91, 731–742, 2008
[^3]: > Goodwin, B.C. [*Oscillatory behavior in enzymatic control processes*](https://doi.org/10.1016/0065-2571(65)90067-1), Advances in Enzyme Regulation, Vol 3 (C), 425-437, 1965
[^4]: > Dong, R., Goodbrake, C., Harrington, H. A., & Pogudin, G. [*Computing input-output projections of dynamical models with applications to structural identifiability*](https://arxiv.org/pdf/2111.00991). arXiv preprint arXiv:2111.00991.
