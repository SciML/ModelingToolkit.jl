# Automated Sparse Parallelism of ODEs via Tracing

Because the ModelingToolkit `Expression` types obey Julia semantics, one can
directly transform existing Julia functions into ModelingToolkit symbolic
representations of the function by simply inputting the symbolic values into
the function and using what is returned. For example, let's take [the following
numerical PDE discretization](https://www.stochasticlifestyle.com/solving-systems-stochastic-pdes-using-gpus-julia/):

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
const N = 32
const X = reshape([i for i in 1:N for j in 1:N],N,N)
const Y = reshape([j for i in 1:N for j in 1:N],N,N)
const α₁ = 1.0.*(X.>=4*N/5)

const Mx = Array(Tridiagonal([1.0 for i in 1:N-1],[-2.0 for i in 1:N],[1.0 for i in 1:N-1]))
const My = copy(Mx)
Mx[2,1] = 2.0
Mx[end-1,end] = 2.0
My[1,2] = 2.0
My[end,end-1] = 2.0

# Define the discretized PDE as an ODE function
function f(u,p,t)
    A = u[:,:,1]
    B = u[:,:,2]
    C = u[:,:,3]
    MyA = My*A
    AMx = A*Mx
    DA = @. _DD*(MyA + AMx)
    dA = @. DA + α₁ - β₁*A - r₁*A*B + r₂*C
    dB = @. α₂ - β₂*B - r₁*A*B + r₂*C
    dC = @. α₃ - β₃*C + r₁*A*B - r₂*C
    cat(dA,dB,dC,dims=3)
end
```

We can build the ModelingToolkit version of this model by tracing the
model function:

```julia
# Define the initial condition as normal arrays
@variables u[1:N,1:N,1:3]
du = simplify.(f(u,nothing,0.0))
```

The output, here the in-place modified `du`, is a symbolic representation of
each output of the function. We can then utilize this in the ModelingToolkit
functionality. For example, let's build a parallel version of `f` first:

```julia
fastf = eval(ModelingToolkit.build_function(du,u,
            parallel=ModelingToolkit.MultithreadedForm())[2])
```

Now let's compute the sparse Jacobian function and compile a fast multithreaded version:

```julia
jac = ModelingToolkit.sparsejacobian(vec(du),vec(u))
fjac = eval(ModelingToolkit.build_function(jac,u,
            parallel=ModelingToolkit.MultithreadedForm())[2])
```

It takes awhile for this to generate, but the results will be worth it!
Now let's setup the parabolic PDE to be solved by DifferentialEquations.jl.
We will setup the vanilla version and the sparse multithreaded
version:

```julia
using OrdinaryDiffEq
u0 = zeros(N,N,3)
MyA = zeros(N,N);
AMx = zeros(N,N);
DA = zeros(N,N);
prob = ODEProblem(f!,u0,(0.0,10.0))
fastprob = ODEProblem(ODEFunction((du,u,p,t)->fastf(du,u),
                                   jac = (du,u,p,t) -> fjac(du,u),
                                   jac_prototype = similar(jac,Float64)),
                                   u0,(0.0,10.0))
```

Let's see the timing difference:

```julia
using BenchmarkTools
@btime solve(prob, TRBDF2()) # 33.073 s (895404 allocations: 23.87 GiB)
@btime solve(fastprob, TRBDF2()) # 209.670 ms (8208 allocations: 109.25 MiB)
```

Boom, an automatic 157x acceleration that grows as the size of the problem
increases!
