module MTKBifurcationKitExt

### Preparations ###

# Imports
using ModelingToolkit, Setfield
import BifurcationKit

### Creates BifurcationProblem Overloads ###

# When input is a NonlinearSystem.
function BifurcationKit.BifurcationProblem(nsys::NonlinearSystem,
    u0_bif,
    ps,
    bif_par,
    args...;
    plot_var = nothing,
    record_from_solution = BifurcationKit.record_sol_default,
    jac = true,
    kwargs...)
    # Creates F and J functions.
    ofun = NonlinearFunction(nsys; jac = jac)
    F = ofun.f
    J = jac ? ofun.jac : nothing

    # Computes bifurcation parameter and plot var indexes.
    bif_idx = findfirst(isequal(bif_par), parameters(nsys))
    if !isnothing(plot_var)
        plot_idx = findfirst(isequal(plot_var), states(nsys))
        record_from_solution = (x, p) -> x[plot_idx]
    end

    # Converts the input state guess.
    u0_bif = ModelingToolkit.varmap_to_vars(u0_bif, states(nsys))
    ps = ModelingToolkit.varmap_to_vars(ps, parameters(nsys))

    return BifurcationKit.BifurcationProblem(F,
        u0_bif,
        ps,
        (@lens _[bif_idx]),
        args...;
        record_from_solution = record_from_solution,
        J = J,
        kwargs...)
end

# When input is a ODESystem.
function BifurcationKit.BifurcationProblem(osys::ODESystem, args...; kwargs...)
    nsys = NonlinearSystem([0 ~ eq.rhs for eq in equations(osys)],
        states(osys),
        parameters(osys);
        name = osys.name)
    return BifurcationKit.BifurcationProblem(nsys, args...; kwargs...)
end

end # module
