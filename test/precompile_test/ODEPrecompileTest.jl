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
    
    # Build a simple ODEFunction as part of the module's precompilation.
    const f_bad = system()
    # const f_good = system(; eval_module=@__MODULE__)
end