module MTKBifurcationKitExt

### Preparations ###

# Imports
using ModelingToolkit, Setfield
import BifurcationKit

### Observable Plotting Handling ###

# Functor used when the plotting variable is an observable. Keeps track of the required information for computing the observable's value at each point of the bifurcation diagram.
struct ObservableRecordFromSolution{S, T}
    # The equations determining the observables values.
    obs_eqs::S
    # The index of the observable that we wish to plot.
    target_obs_idx::Int64
    # The final index in subs_vals that contains a state.
    state_end_idxs::Int64
    # The final index in subs_vals that contains a param.
    param_end_idxs::Int64
    # The index (in subs_vals) that contain the bifurcation parameter.
    bif_par_idx::Int64
    # A Vector of pairs (Symbolic => value) with the default values of all system variables and parameters.
    subs_vals::T

    function ObservableRecordFromSolution(nsys::NonlinearSystem,
            plot_var,
            bif_idx,
            u0_vals,
            p_vals) where {S, T}
        obs_eqs = observed(nsys)
        target_obs_idx = findfirst(isequal(plot_var, eq.lhs) for eq in observed(nsys))
        state_end_idxs = length(states(nsys))
        param_end_idxs = state_end_idxs + length(parameters(nsys))

        bif_par_idx = state_end_idxs + bif_idx
        # Gets the (base) substitution values for states. 
        subs_vals_states = Pair.(states(nsys), u0_vals)
        # Gets the (base) substitution values for parameters. 
        subs_vals_params = Pair.(parameters(nsys), p_vals)
        # Gets the (base) substitution values for observables. 
        subs_vals_obs = [obs.lhs => substitute(obs.rhs,
            [subs_vals_states; subs_vals_params]) for obs in observed(nsys)]
        # Sometimes observables depend on other observables, hence we make a second update to this vector.
        subs_vals_obs = [obs.lhs => substitute(obs.rhs,
            [subs_vals_states; subs_vals_params; subs_vals_obs]) for obs in observed(nsys)]
        # During the bifurcation process, the value of some states, parameters, and observables may vary (and are calculated in each step). Those that are not are stored in this vector
        subs_vals = [subs_vals_states; subs_vals_params; subs_vals_obs]

        param_end_idxs = state_end_idxs + length(parameters(nsys))
        new{typeof(obs_eqs), typeof(subs_vals)}(obs_eqs,
            target_obs_idx,
            state_end_idxs,
            param_end_idxs,
            bif_par_idx,
            subs_vals)
    end
end
# Functor function that computes the value.
function (orfs::ObservableRecordFromSolution)(x, p)
    # Updates the state values (in subs_vals).
    for state_idx in 1:(orfs.state_end_idxs)
        orfs.subs_vals[state_idx] = orfs.subs_vals[state_idx][1] => x[state_idx]
    end

    # Updates the bifurcation parameters value (in subs_vals).
    orfs.subs_vals[orfs.bif_par_idx] = orfs.subs_vals[orfs.bif_par_idx][1] => p

    # Updates the observable values (in subs_vals).
    for (obs_idx, obs_eq) in enumerate(orfs.obs_eqs)
        orfs.subs_vals[orfs.param_end_idxs + obs_idx] = orfs.subs_vals[orfs.param_end_idxs + obs_idx][1] => substitute(obs_eq.rhs,
            orfs.subs_vals)
    end

    # Substitutes in the value for all states, parameters, and observables into the equation for the designated observable.
    return substitute(orfs.obs_eqs[orfs.target_obs_idx].rhs, orfs.subs_vals)
end

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

    # Converts the input state guess.
    u0_bif_vals = ModelingToolkit.varmap_to_vars(u0_bif,
        states(nsys);
        defaults = nsys.defaults)
    p_vals = ModelingToolkit.varmap_to_vars(ps, parameters(nsys); defaults = nsys.defaults)

    # Computes bifurcation parameter and the plotting function.
    bif_idx = findfirst(isequal(bif_par), parameters(nsys))
    if !isnothing(plot_var)
        # If the plot var is a normal state.
        if any(isequal(plot_var, var) for var in states(nsys))
            plot_idx = findfirst(isequal(plot_var), states(nsys))
            record_from_solution = (x, p) -> x[plot_idx]

            # If the plot var is an observed state.
        elseif any(isequal(plot_var, eq.lhs) for eq in observed(nsys))
            record_from_solution = ObservableRecordFromSolution(nsys,
                plot_var,
                bif_idx,
                u0_bif_vals,
                p_vals)

            # If neither an variable nor observable, throw an error.
        else
            error("The plot variable ($plot_var) was neither recognised as a system state nor observable.")
        end
    end

    return BifurcationKit.BifurcationProblem(F,
        u0_bif_vals,
        p_vals,
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
