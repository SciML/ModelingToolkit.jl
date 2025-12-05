MTKBase.singular_check(ts::TearingState) = StateSelection.singular_check(ts)

function MTKBase.get_initialization_problem_type(sys::System, isys::System;
                                                 warn_initialize_determined = true,
                                                 use_scc = true, kwargs...)
    neqs = length(equations(isys))
    nunknown = length(unknowns(isys))
    ts = get_tearing_state(isys)::TearingState

    if use_scc
        scc_message = """
        `SCCNonlinearProblem` can only be used for initialization of fully determined \
        systems and hence will not be used here.
        """
    else
        scc_message = ""
    end

    if warn_initialize_determined && neqs > nunknown
        @warn overdetermined_initialization_message(neqs, nunknown, scc_message)
    end
    if warn_initialize_determined && neqs < nunknown
        @warn underdetermined_initialization_message(neqs, nunknown, scc_message)
    end

    unassigned_vars = MTKBase.singular_check(ts)
    if neqs == nunknown && isempty(unassigned_vars)
        if use_scc && neqs > 0
            if is_split(isys)
                SCCNonlinearProblem
            else
                @warn """
                `SCCNonlinearProblem` can only be used with `split = true` systems. \
                Simplify your `System` with `split = true` or pass `use_scc = false` to \
                disable this warning
                """
                NonlinearProblem
            end
        else
            NonlinearProblem
        end
    else
        NonlinearLeastSquaresProblem
    end
end

MTKBase.default_missing_guess_value(::Nothing) = MTKBase.MissingGuessValue.Error()
