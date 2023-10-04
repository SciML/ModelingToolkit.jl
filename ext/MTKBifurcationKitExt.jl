module MTKBifurcationKitExt

println("BifurcationKit extension loaded")

### Preparations ###

# Imports
using ModelingToolkit, Setfield
import BifurcationKit: BifurcationProblem

### Creates BifurcationProblem Overloads ###

# When input is a NonlinearProblem.
function BifurcationKit.BifurcationProblem(nprob::NonlinearProblem, u0, p, bif_par, args...; plot_variable=nothing, record_from_solution=record_sol_default, kwargs...)
    # check if jac does nto exist and modify.
    F = nprob.f.f
    J = nprob.f.jac.f

    bif_idx = findfirst(isequal(bif_par), parameters(nsys))
    if isnothing(plot_variable) 
        plot_idx = findfirst(isequal(plot_variable), variables(nsys))
        record_from_solution = (x, p) -> x[plot_idx]
    end
    return BifurcationProblem(F, u0, p, (@lens _[bif_idx]), args...; record_from_solution = record_from_solution, J = J, kwargs...)
end

# When input is a NonlinearSystem.
function BifurcationKit.BifurcationProblem(nsys::NonlinearSystem, u0, p, args...; kwargs...)
    return BifurcationProblem(NonlinearProblem(nsys, u0, p; jac=true), args...; kwargs...)
end

# When input is a ODEProblem.
function BifurcationKit.BifurcationProblem(oprob::ODEProblem, u0, p, args...; kwargs...)
    return BifurcationProblem(convert(NonlinearProblem, oprob), args...; kwargs...)
end

# When input is a ODESystem.
function BifurcationKit.BifurcationProblem(osys::ODESystem, u0, p, args...; kwargs...)
    return BifurcationProblem(convert(NonlinearProblem, osys), args...; kwargs...)
end



end # module
