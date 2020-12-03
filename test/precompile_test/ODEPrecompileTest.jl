module ODEPrecompileTest
    using ModelingToolkit

    function system(; kwargs...)
        # Define some variables
        @parameters t σ ρ β
        @variables x(t) y(t) z(t)
        @derivatives D'~t

        # Define a differential equation
        eqs = [D(x) ~ σ*(y-x),
            D(y) ~ x*(ρ-z)-y,
            D(z) ~ x*y - β*z]

        de = ODESystem(eqs)
        return ODEFunction(de, [x,y,z], [σ,ρ,β]; kwargs...)
    end
    
    # Build an ODEFunction as part of the module's precompilation. This case
    # will not work, because the generated RGFs will be put into
    # ModelingToolkit's RGF cache.
    const f_bad = system()

    # This case will work, because it will be put into our own module's cache.
    using RuntimeGeneratedFunctions
    RuntimeGeneratedFunctions.init(@__MODULE__)
    const f_good = system(; eval_module=@__MODULE__)
end