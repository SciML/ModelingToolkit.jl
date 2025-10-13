module ODEPrecompileTest
using ModelingToolkit

function system(; kwargs...)
    # Define some variables
    @independent_variables t
    @parameters σ ρ β
    @variables x(t) y(t) z(t)
    D = Differential(t)

    # Define a differential equation
    eqs = [D(x) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z]

    @named de = System(eqs, t)
    de = complete(de)
    return ODEFunction(de; kwargs...)
end

# Build an ODEFunction as part of the module's precompilation. These cases
# will not work, because the generated RGFs are put into the ModelingToolkit cache.
const f_bad = system()
const f_noeval_bad = system(; eval_expression = false)

# Setting eval_expression=false and eval_module=[this module] will ensure
# the RGFs are put into our own cache, initialised below.
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)
const f_noeval_good = system(; eval_expression = false, eval_module = @__MODULE__)

# Eval the expression but into MTK's module, which means it won't be properly cached by
# the package image
const f_eval_bad = system(; eval_expression = true, eval_module = @__MODULE__)

# Change the module the eval'd function is eval'd into to be the containing module,
# which should make it be in the package image
const f_eval_good = system(; eval_expression = true, eval_module = @__MODULE__)

function problem(; kwargs...)
    # Define some variables
    @independent_variables t
    @parameters σ ρ β
    @variables x(t) y(t) z(t)
    D = Differential(t)

    # Define a differential equation
    eqs = [D(x) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z]

    @named de = System(eqs, t)
    de = complete(de)
    return ODEProblem(de, [x => 1, y => 0, z => 0, σ => 10, ρ => 28, β => 8/3], (0.0, 5.0); kwargs...)
end

const prob_eval = problem(; eval_expression = true, eval_module = @__MODULE__)

end
