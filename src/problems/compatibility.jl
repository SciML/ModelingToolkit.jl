"""
    function check_compatible_system(T::Type, sys::System)

Check if `sys` can be used to construct a problem/function of type `T`.
"""
function check_compatible_system end

struct SystemCompatibilityError <: Exception
    msg::String
end

function Base.showerror(io::IO, err::SystemCompatibilityError)
    println(io, err.msg)
    println(io)
    print(io, "To disable this check, pass `check_compatibility = false`.")
end

function check_time_dependent(sys::System, T)
    if !is_time_dependent(sys)
        throw(SystemCompatibilityError("""
        `$T` requires a time-dependent system.
        """))
    end
end

function check_time_independent(sys::System, T)
    if is_time_dependent(sys)
        throw(SystemCompatibilityError("""
        `$T` requires a time-independent system.
        """))
    end
end

function check_is_dde(sys::System)
    altT = get_noise_eqs(sys) === nothing ? ODEProblem : SDEProblem
    if !is_dde(sys)
        throw(SystemCompatibilityError("""
        The system does not have delays. Consider an `$altT` instead.
        """))
    end
end

function check_not_dde(sys::System)
    altT = get_noise_eqs(sys) === nothing ? DDEProblem : SDDEProblem
    if is_dde(sys)
        throw(SystemCompatibilityError("""
        The system has delays. Consider a `$altT` instead.
        """))
    end
end

function check_no_cost(sys::System, T)
    cost = ModelingToolkit.cost(sys)
    if !_iszero(cost)
        throw(SystemCompatibilityError("""
        `$T` will not optimize solutions of systems that have associated cost \
        functions. Solvers for optimal control problems are forthcoming.
        """))
    end
end

function check_has_cost(sys::System, T)
    cost = ModelingToolkit.cost(sys)
    if _iszero(cost)
        throw(SystemCompatibilityError("""
        A system without cost cannot be used to construct a `$T`.
        """))
    end
end

function check_no_constraints(sys::System, T)
    if !isempty(constraints(sys))
        throw(SystemCompatibilityError("""
        A system with constraints cannot be used to construct a `$T`.
        """))
    end
end

function check_has_constraints(sys::System, T)
    if isempty(constraints(sys))
        throw(SystemCompatibilityError("""
        A system without constraints cannot be used to construct a `$T`. Consider an \
        `ODEProblem` instead.
        """))
    end
end

function check_no_jumps(sys::System, T)
    if !isempty(jumps(sys))
        throw(SystemCompatibilityError("""
        A system with jumps cannot be used to construct a `$T`. Consider a \
        `JumpProblem` instead.
        """))
    end
end

function check_has_jumps(sys::System, T)
    if isempty(jumps(sys))
        throw(SystemCompatibilityError("`$T` requires a system with jumps."))
    end
end

function check_no_noise(sys::System, T)
    altT = is_dde(sys) ? SDDEProblem : SDEProblem
    if get_noise_eqs(sys) !== nothing
        throw(SystemCompatibilityError("""
        A system with noise cannot be used to construct a `$T`. Consider an \
        `$altT` instead.
        """))
    end
end

function check_has_noise(sys::System, T)
    altT = is_dde(sys) ? DDEProblem : ODEProblem
    if get_noise_eqs(sys) === nothing
        msg = """
        A system without noise cannot be used to construct a `$T`. Consider an \
        `$altT` instead.
        """
        if !isempty(brownians(sys))
            msg = """
            Systems constructed by defining Brownian variables with `@brownians` must be \
            simplified by calling `mtkcompile` before a `$T` can be constructed.
            """
        end
        throw(SystemCompatibilityError(msg))
    end
end

function check_is_discrete(sys::System, T)
    if !is_discrete_system(sys)
        throw(SystemCompatibilityError("""
        `$T` expects a discrete system. Consider an `ODEProblem` instead. If your system \
        is discrete, ensure `mtkcompile` has been run on it.
        """))
    end
end

function check_is_continuous(sys::System, T)
    altT = has_alg_equations(sys) ? ImplicitDiscreteProblem : DiscreteProblem
    if is_discrete_system(sys)
        throw(SystemCompatibilityError("""
        A discrete system cannot be used to construct a `$T`. Consider a `$altT` instead.
        """))
    end
end

function check_is_explicit(sys::System, T, altT)
    if has_alg_equations(sys)
        throw(SystemCompatibilityError("""
        `$T` expects an explicit system. Consider a `$altT` instead.
        """))
    end
end

function check_is_implicit(sys::System, T, altT)
    if !has_alg_equations(sys)
        throw(SystemCompatibilityError("""
        `$T` expects an implicit system. Consider a `$altT` instead.
        """))
    end
end

function check_no_equations(sys::System, T)
    if !isempty(equations(sys))
        throw(SystemCompatibilityError("""
        A system with equations cannot be used to construct a `$T`. Consider turning the
        equations into constraints instead.
        """))
    end
end

function check_affine(sys::System, T)
    if !isaffine(sys)
        throw(SystemCompatibilityError("""
        A non-affine system cannot be used to construct a `$T`. Consider a
        `NonlinearProblem` instead.
        """))
    end
end
