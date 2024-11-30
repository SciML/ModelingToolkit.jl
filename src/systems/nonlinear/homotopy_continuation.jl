"""
$(TYPEDEF)

A type of Nonlinear problem which specializes on polynomial systems and uses
HomotopyContinuation.jl to solve the system. Requires importing HomotopyContinuation.jl to
create and solve.
"""
struct HomotopyContinuationProblem{uType, H, D, O, SS, U} <:
       SciMLBase.AbstractNonlinearProblem{uType, true}
    """
    The initial values of states in the system. If there are multiple real roots of
    the system, the one closest to this point is returned.
    """
    u0::uType
    """
    A subtype of `HomotopyContinuation.AbstractSystem` to solve. Also contains the
    parameter object.
    """
    homotopy_continuation_system::H
    """
    A function with signature `(u, p) -> resid`. In case of rational functions, this
    is used to rule out roots of the system which would cause the denominator to be
    zero.
    """
    denominator::D
    """
    The `NonlinearSystem` used to create this problem. Used for symbolic indexing.
    """
    sys::NonlinearSystem
    """
    A function which generates and returns observed expressions for the given system.
    """
    obsfn::O
    """
    The HomotopyContinuation.jl solver and start system, obtained through
    `HomotopyContinuation.solver_startsystems`.
    """
    solver_and_starts::SS
    """
    A function which takes a solution of the transformed system, and returns a vector
    of solutions for the original system. This is utilized when converting systems
    to polynomials.
    """
    unpack_solution::U
end

function HomotopyContinuationProblem(::AbstractSystem, _u0, _p; kwargs...)
    error("HomotopyContinuation.jl is required to create and solve `HomotopyContinuationProblem`s. Please run `Pkg.add(\"HomotopyContinuation\")` to continue.")
end

SymbolicIndexingInterface.symbolic_container(p::HomotopyContinuationProblem) = p.sys
SymbolicIndexingInterface.state_values(p::HomotopyContinuationProblem) = p.u0
function SymbolicIndexingInterface.set_state!(p::HomotopyContinuationProblem, args...)
    set_state!(p.u0, args...)
end
function SymbolicIndexingInterface.parameter_values(p::HomotopyContinuationProblem)
    parameter_values(p.homotopy_continuation_system)
end
function SymbolicIndexingInterface.set_parameter!(p::HomotopyContinuationProblem, args...)
    set_parameter!(parameter_values(p), args...)
end
function SymbolicIndexingInterface.observed(p::HomotopyContinuationProblem, sym)
    if p.obsfn !== nothing
        return p.obsfn(sym)
    else
        return SymbolicIndexingInterface.observed(p.sys, sym)
    end
end
