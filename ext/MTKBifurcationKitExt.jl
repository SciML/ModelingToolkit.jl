module MTKBifurcationKitExt

println("BifurcationKit extension loaded")

### Preparations ###

# Imports
using ModelingToolkit, Setfield
import BifurcationKit

### Creates BifurcationProblem Overloads ###

# When input is a NonlinearSystem.
function BifurcationKit.BifurcationProblem(nsys::NonlinearSystem, u0_guess, p_start, bif_par, args...; plot_var=nothing, record_from_solution=BifurcationKit.record_sol_default, kwargs...)
    # Creates F and J functions.
    ofun = NonlinearFunction(nsys; jac=true)
    F = ofun.f
    J = ofun.jac

    # Computes bifurcation parameter and plot var indexes.
    bif_idx = findfirst(isequal(bif_par), parameters(nsys))
    if !isnothing(plot_var) 
        plot_idx = findfirst(isequal(plot_var), states(nsys))
        record_from_solution = (x, p) -> x[plot_idx]
    end

    # Converts the input state guess.
    u0_guess = ModelingToolkit.varmap_to_vars(u0_guess, states(nsys))

    return BifurcationKit.BifurcationProblem(F, u0_guess, [p_start], (@lens _[1]), args...; record_from_solution = record_from_solution, J = J, kwargs...)
end

# When input is a ODESystem.
function BifurcationKit.BifurcationProblem(osys::ODESystem, u0, p, args...; kwargs...)
    return BifurcationProblem(convert(NonlinearProblem, osys), args...; kwargs...)
end



end # module
