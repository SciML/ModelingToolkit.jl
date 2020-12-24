using Revise # DELETE ME

using Test
using ModelingToolkit
using LinearAlgebra.BLAS

@parameters t
@variables vars[1:10]
@derivatives D'~t

const sleep_time = 0.1
function untrace_dupe_f(a::Float64, t::Float64)
    # This is a function that is both expensive and untraceable (symbolically),
    # and whose result is used several times in our ODE system of equations.
    sleep(sleep_time)  # Expensive function!
    BLAS.dot(10, fill(a, 10), 1, fill(t, 20), 2)  # Untraceable result!
end
@test untrace_dupe_f(5., 12.) == 600.0

# Confirm that our function cannot be traced symbolically...
@test_throws MethodError untrace_dupe_f(vars[1], t)

# Register the function so we can now use it without having to trace it.
@register untrace_dupe_f(a, t)

# Build an ODE system of equations that uses several duplicated, untraceable terms.
dupe_terms = (untrace_dupe_f(x, t) for x in vars)
eqs = (D.(vars) .~ prod(dupe_terms))
@test length(eqs) == length(vars)

# Build the ODE functions.
ode_system = ODESystem(eqs, t, vars, [])
ode_function_naive = ODEFunction(ode_system; jac=false, tgrad=false, eval_expression=false)
ode_function_deduplicated = ODEFunction(ode_system; jac=false, tgrad=false, eval_expression=false, deduplicate_terms=dupe_terms)

# Run both functions and compare...
u0 = rand(Float64, length(x))
t0 = 10.0
result_naive = ode_function_naive(u0, [], t0)
result_deduplicated = ode_function_deduplicated(u0, [], t0)
@test result_naive == result_deduplicated

# The deduplicated version should run much faster, since it does not need to recompute the terms multiple times.
fuzzy_factor = 0.9
expected_speedup = fuzzy_factor * length(eqs)
@test @elapsed(ode_function_naive(u0, [], t0)) > (expected_speedup * @elapsed(ode_function_deduplicated(u0, [], t0)))