# Solving Nonlinear Systems with NLsolve

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
