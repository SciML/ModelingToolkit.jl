# Automated Sparse Analytical Jacobians

In many cases where you have large stiff differential equations, getting a
sparse Jacobian can be essential for performance. In this tutorial, we will show
how to use `modelingtoolkitize` to regenerate an `ODEProblem` code with
the analytical solution to the sparse Jacobian, along with the sparsity
pattern required by DifferentialEquations.jl's solvers to specialize the solving
process.

First, let's start out with an implementation of the 2-dimensional Brusselator
partial differential equation discretized using finite differences:

```@example sparsejac
using OrdinaryDiffEq, ModelingToolkit

const N = 32
const xyd_brusselator = range(0, stop = 1, length = N)
brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0
limit(a, N) = a == N + 1 ? 1 : a == 0 ? N : a
function brusselator_2d_loop(du, u, p, t)
    A, B, alpha, dx = p
    alpha = alpha / dx^2
    @inbounds for I in CartesianIndices((N, N))
        i, j = Tuple(I)
        x, y = xyd_brusselator[I[1]], xyd_brusselator[I[2]]
        ip1, im1, jp1, jm1 = limit(i + 1, N), limit(i - 1, N), limit(j + 1, N),
        limit(j - 1, N)
        du[i, j, 1] = alpha * (u[im1, j, 1] + u[ip1, j, 1] + u[i, jp1, 1] + u[i, jm1, 1] -
                       4u[i, j, 1]) +
                      B + u[i, j, 1]^2 * u[i, j, 2] - (A + 1) * u[i, j, 1] +
                      brusselator_f(x, y, t)
        du[i, j, 2] = alpha * (u[im1, j, 2] + u[ip1, j, 2] + u[i, jp1, 2] + u[i, jm1, 2] -
                       4u[i, j, 2]) +
                      A * u[i, j, 1] - u[i, j, 1]^2 * u[i, j, 2]
    end
end
p = (3.4, 1.0, 10.0, step(xyd_brusselator))

function init_brusselator_2d(xyd)
    N = length(xyd)
    u = zeros(N, N, 2)
    for I in CartesianIndices((N, N))
        x = xyd[I[1]]
        y = xyd[I[2]]
        u[I, 1] = 22 * (y * (1 - y))^(3 / 2)
        u[I, 2] = 27 * (x * (1 - x))^(3 / 2)
    end
    u
end
u0 = init_brusselator_2d(xyd_brusselator)
prob = ODEProblem(brusselator_2d_loop, u0, (0.0, 11.5), p)
```

Now let's use `modelingtoolkitize` to generate the symbolic version:

```@example sparsejac
@mtkbuild sys = modelingtoolkitize(prob);
nothing # hide
```

Now we regenerate the problem using `jac=true` for the analytical Jacobian
and `sparse=true` to make it sparse:

```@example sparsejac
sparseprob = ODEProblem(sys, Pair[], (0.0, 11.5), jac = true, sparse = true)
```

Hard? No! How much did that help?

```@example sparsejac
using BenchmarkTools
@btime solve(prob, save_everystep = false);
return nothing # hide
```

```@example sparsejac
@btime solve(sparseprob, save_everystep = false);
return nothing # hide
```

Notice though that the analytical solution to the Jacobian can be quite expensive.
Thus in some cases we may only want to get the sparsity pattern. In this case,
we can simply do:

```@example sparsejac
sparsepatternprob = ODEProblem(sys, Pair[], (0.0, 11.5), sparse = true)
@btime solve(sparsepatternprob, save_everystep = false);
return nothing # hide
```
