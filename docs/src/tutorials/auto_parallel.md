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

We can build the ModelingToolkit version of this model by tracing the
model function:

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
jac = sparse(ModelingToolkit.jacobian(vec(du),vec(u)))
fjac = eval(ModelingToolkit.build_function(jac,u,parallel=ModelingToolkit.MultithreadedForm())[2])
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
prob_jac = ODEProblem(ODEFunction(f!,jac = (du,u,p,t) -> fjac(du,u), jac_prototype = similar(jac,Float64)),u0,(0.0,10.0))
```

Let's see the timing difference:

```julia
using BenchmarkTools
@btime solve(prob, TRBDF2(autodiff=false)) # 5.995 s (85786 allocations: 149.26 MiB)
@btime solve(prob_jac, TRBDF2()) # 250.409 ms (7411 allocations: 97.73 MiB)
```

Boom, and automatic 24x acceleration that grows as the size of the problem
increases!
