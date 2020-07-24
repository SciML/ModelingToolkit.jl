using ModelingToolkit

@parameters t σ ρ β
@variables x(t) y(t) z(t) a(t) p(t) q(t)
@derivatives D'~t

# Define a differential equation
eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(q) ~ tan(x)*sin((ρ-z)-y),
       D(p) ~ exp(x)*(ρ-cos(z))-y,
       D(a) ~ x*(ρ-z)-sin(y),
       D(z) ~ x*y - β*z]

ModelingToolkit.simplified_expr.(eqs)[1]
:(derivative(x(t), t) = σ * (y(t) - x(t))).args
de = ODESystem(eqs)

using BenchmarkTools



display(@benchmark calculate_jacobian($de))
display(@benchmark calculate_factorized_W($de))

using Profile

de = ODESystem(eqs)
@profile calculate_factorized_W(de)

# master
#
#BenchmarkTools.Trial: 
#  memory estimate:  192 bytes
#  allocs estimate:  2
#  --------------
#  minimum time:     584.940 ns (0.00% GC)
#  median time:      726.746 ns (0.00% GC)
#  mean time:        806.295 ns (0.77% GC)
#  maximum time:     14.499 μs (94.09% GC)
#  --------------
#  samples:          10000
#  evals/sample:     183
#BenchmarkTools.Trial: 
#  memory estimate:  320 bytes
#  allocs estimate:  4
#  --------------
#  minimum time:     981.600 ns (0.00% GC)
#  median time:      1.306 μs (0.00% GC)
#  mean time:        1.474 μs (0.00% GC)
#  maximum time:     8.318 μs (0.00% GC)
#  --------------
#  samples:          10000
#  evals/sample:     10
#
# symutils
#
# BenchmarkTools.Trial: 
#   memory estimate:  192 bytes
#   allocs estimate:  2
#   --------------
#   minimum time:     571.560 ns (0.00% GC)
#   median time:      708.547 ns (0.00% GC)
#   mean time:        782.288 ns (1.41% GC)
#   maximum time:     20.184 μs (95.86% GC)
#   --------------
#   samples:          10000
#   evals/sample:     182
# BenchmarkTools.Trial: 
#   memory estimate:  320 bytes
#   allocs estimate:  4
#   --------------
#   minimum time:     1.003 μs (0.00% GC)
#   median time:      1.193 μs (0.00% GC)
#   mean time:        1.380 μs (0.00% GC)
#   maximum time:     6.631 μs (0.00% GC)
#   --------------
#   samples:          10000
#   evals/sample:     10
#
#
# - Master
#
# julia> @btime calculate_factorized_W(de) setup = (de = ODESystem(eqs);)
#   4.698 s (19477345 allocations: 615.74 MiB)
# (Expression[-1 + -1 * __MTKWgamma * σ __MTKWgamma * σ … ModelingToolkit.Constant(0) ModelingToolkit.Constant(0); __MTKWgamma * inv(-1 + -1 * __MTKWgamma * σ) * one(x(t) * (-1 * z(t) + ρ)) * (-1 * z(t) + ρ) -1 + -1 * (__MTKWgamma * inv(-1 + -1 * __MTKWgamma * σ) * one(x(t) * (-1 * z(t) + ρ)) * (-1 * z(t) + ρ) * __MTKWgamma * σ) + -(one(y(t))) * __MTKWgamma … ModelingToolkit.Constant(0) -1 * __MTKWgamma * one(x(t) * (-1 * z(t) + ρ)) * x(t); … ; __MTKWgamma * inv(-1 + -1 * __MTKWgamma * σ) * one(x(t) * (-1 * z(t) + ρ)) * (-1 * z(t) + ρ) inv(-1 + -1 * (__MTKWgamma * inv(-1 + -1 * __MTKWgamma * σ) * one(x(t) * (-1 * z(t) + ρ)) * (-1 * z(t) + ρ) * __MTKWgamma * σ) + -(one(y(t))) * __MTKWgamma) * (-1 * (__MTKWgamma * inv(-1 + -1 * __MTKWgamma * σ) * one(x(t) * (-1 * z(t) + ρ)) * (-1 * z(t) + ρ) * __MTKWgamma * σ) + -(one(sin(y(t)))) * __MTKWgamma * cos(y(t))) … ModelingToolkit.Constant(-1) -1 * (-1 * inv(-1 + -1 * (__MTKWgamma * inv(-1 + -1 * __MTKWgamma * σ) * one(x(t) * (-1 * z(t) + ρ)) * (-1 * z(t) + ρ) * __MTKWgamma * σ) + -(one(y(t))) * __MTKWgamma) * (-1 * (__MTKWgamma * inv(-1 + -1 * __MTKWgamma * σ) * one(x(t) * (-1 * z(t) + ρ)) * (-1 * z(t) + ρ) * __MTKWgamma * σ) + -(one(sin(y(t)))) * __MTKWgamma * cos(y(t))) * __MTKWgamma * one(x(t) * (-1 * z(t) + ρ)) * x(t)) + -1 * __MTKWgamma * one(x(t) * (-1 * z(t) + ρ)) * x(t); __MTKWgamma * inv(-1 + -1 * __MTKWgamma * σ) * one(x(t) * y(t)) * y(t) inv(-1 + -1 * (__MTKWgamma * inv(-1 + -1 * __MTKWgamma * σ) * one(x(t) * (-1 * z(t) + ρ)) * (-1 * z(t) + ρ) * __MTKWgamma * σ) + -(one(y(t))) * __MTKWgamma) * (-1 * (inv(-1 + -1 * __MTKWgamma * σ) * one(x(t) * y(t)) * y(t) * __MTKWgamma ^ 2 * σ) + __MTKWgamma * one(x(t) * y(t)) * x(t)) … ModelingToolkit.Constant(0) -1 + -1 * (-1 * inv(-1 + -1 * (__MTKWgamma * inv(-1 + -1 * __MTKWgamma * σ) * one(x(t) * (-1 * z(t) + ρ)) * (-1 * z(t) + ρ) * __MTKWgamma * σ) + -(one(y(t))) * __MTKWgamma) * (-1 * (inv(-1 + -1 * __MTKWgamma * σ) * one(x(t) * y(t)) * y(t) * __MTKWgamma ^ 2 * σ) + __MTKWgamma * one(x(t) * y(t)) * x(t)) * __MTKWgamma * one(x(t) * (-1 * z(t) + ρ)) * x(t)) + -(one(z(t) * β)) * __MTKWgamma * β], Expression[-1 * __MTKWgamma ^ -1 + -1σ σ … ModelingToolkit.Constant(0) ModelingToolkit.Constant(0); one(x(t) * (-1 * z(t) + ρ)) * inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) -(one(y(t))) + -1 * (one(x(t) * (-1 * z(t) + ρ)) * inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -1 * __MTKWgamma ^ -1 … ModelingToolkit.Constant(0) -1 * one(x(t) * (-1 * z(t) + ρ)) * x(t); … ; one(x(t) * (-1 * z(t) + ρ)) * inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) inv(-(one(y(t))) + -1 * (one(x(t) * (-1 * z(t) + ρ)) * inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -1 * __MTKWgamma ^ -1) * (-1 * (one(x(t) * (-1 * z(t) + ρ)) * inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -(one(sin(y(t)))) * cos(y(t))) … -1 * __MTKWgamma ^ -1 -1 * (-1 * one(x(t) * (-1 * z(t) + ρ)) * inv(-(one(y(t))) + -1 * (one(x(t) * (-1 * z(t) + ρ)) * inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -1 * __MTKWgamma ^ -1) * x(t) * (-1 * (one(x(t) * (-1 * z(t) + ρ)) * inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -(one(sin(y(t)))) * cos(y(t)))) + -1 * one(x(t) * (-1 * z(t) + ρ)) * x(t); one(x(t) * y(t)) * inv(-1 * __MTKWgamma ^ -1 + -1σ) * y(t) inv(-(one(y(t))) + -1 * (one(x(t) * (-1 * z(t) + ρ)) * inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -1 * __MTKWgamma ^ -1) * (-1 * (one(x(t) * y(t)) * inv(-1 * __MTKWgamma ^ -1 + -1σ) * y(t) * σ) + one(x(t) * y(t)) * x(t)) … ModelingToolkit.Constant(0) -1 * (-1 * one(x(t) * (-1 * z(t) + ρ)) * inv(-(one(y(t))) + -1 * (one(x(t) * (-1 * z(t) + ρ)) * inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -1 * __MTKWgamma ^ -1) * x(t) * (-1 * (one(x(t) * y(t)) * inv(-1 * __MTKWgamma ^ -1 + -1σ) * y(t) * σ) + one(x(t) * y(t)) * x(t))) + -(one(z(t) * β)) * β + -1 * __MTKWgamma ^ -1])
#
# symutils
#
# julia> @btime calculate_factorized_W(de) setup = (de = ODESystem(eqs);)
#   3.712 s (15254171 allocations: 480.70 MiB)
# (Expression[-1 + -1 * __MTKWgamma * σ __MTKWgamma * σ … ModelingToolkit.Constant(0) ModelingToolkit.Constant(0); (-1 * z(t) + ρ) * __MTKWgamma * inv(-1 + -1 * __MTKWgamma * σ) -1 + -1 * (inv(-1 + -1 * __MTKWgamma * σ) * (-1 * z(t) + ρ) * __MTKWgamma ^ 2 * σ) + -1__MTKWgamma … ModelingToolkit.Constant(0) -1 * __MTKWgamma * x(t); … ; (-1 * z(t) + ρ) * __MTKWgamma * inv(-1 + -1 * __MTKWgamma * σ) inv(-1 + -1 * (inv(-1 + -1 * __MTKWgamma * σ) * (-1 * z(t) + ρ) * __MTKWgamma ^ 2 * σ) + -1__MTKWgamma) * (-1 * (inv(-1 + -1 * __MTKWgamma * σ) * (-1 * z(t) + ρ) * __MTKWgamma ^ 2 * σ) + -1 * __MTKWgamma * cos(y(t))) … ModelingToolkit.Constant(-1) -1 * (-1 * inv(-1 + -1 * (inv(-1 + -1 * __MTKWgamma * σ) * (-1 * z(t) + ρ) * __MTKWgamma ^ 2 * σ) + -1__MTKWgamma) * (-1 * (inv(-1 + -1 * __MTKWgamma * σ) * (-1 * z(t) + ρ) * __MTKWgamma ^ 2 * σ) + -1 * __MTKWgamma * cos(y(t))) * __MTKWgamma * x(t)) + -1 * __MTKWgamma * x(t); __MTKWgamma * inv(-1 + -1 * __MTKWgamma * σ) * y(t) inv(-1 + -1 * (inv(-1 + -1 * __MTKWgamma * σ) * (-1 * z(t) + ρ) * __MTKWgamma ^ 2 * σ) + -1__MTKWgamma) * (-1 * (inv(-1 + -1 * __MTKWgamma * σ) * y(t) * __MTKWgamma ^ 2 * σ) + __MTKWgamma * x(t)) … ModelingToolkit.Constant(0) -1 + -1 * (-1 * inv(-1 + -1 * (inv(-1 + -1 * __MTKWgamma * σ) * (-1 * z(t) + ρ) * __MTKWgamma ^ 2 * σ) + -1__MTKWgamma) * (-1 * (inv(-1 + -1 * __MTKWgamma * σ) * y(t) * __MTKWgamma ^ 2 * σ) + __MTKWgamma * x(t)) * __MTKWgamma * x(t)) + -1 * __MTKWgamma * β], Expression[-1 * __MTKWgamma ^ -1 + -1σ σ … ModelingToolkit.Constant(0) ModelingToolkit.Constant(0); inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) -1 + -1 * (inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -1 * __MTKWgamma ^ -1 … ModelingToolkit.Constant(0) -1 * x(t); … ; inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) inv(-1 + -1 * (inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -1 * __MTKWgamma ^ -1) * (-1 * (inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -1 * cos(y(t))) … -1 * __MTKWgamma ^ -1 -1 * (-1 * inv(-1 + -1 * (inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -1 * __MTKWgamma ^ -1) * x(t) * (-1 * (inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -1 * cos(y(t)))) + -1 * x(t); inv(-1 * __MTKWgamma ^ -1 + -1σ) * y(t) inv(-1 + -1 * (inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -1 * __MTKWgamma ^ -1) * (x(t) + -1 * (inv(-1 * __MTKWgamma ^ -1 + -1σ) * y(t) * σ)) … ModelingToolkit.Constant(0) -1 * (-1 * inv(-1 + -1 * (inv(-1 * __MTKWgamma ^ -1 + -1σ) * (-1 * z(t) + ρ) * σ) + -1 * __MTKWgamma ^ -1) * x(t) * (x(t) + -1 * (inv(-1 * __MTKWgamma ^ -1 + -1σ) * y(t) * σ))) + -1 * __MTKWgamma ^ -1 + -1β])
