module MTKHomotopyContinuationExt

using ModelingToolkit
using ModelingToolkit.SciMLBase
using ModelingToolkit.Symbolics: unwrap
using ModelingToolkit.SymbolicIndexingInterface
using ModelingToolkit.DocStringExtensions
using HomotopyContinuation
using ModelingToolkit: iscomplete, parameters, has_index_cache, get_index_cache, get_u0,
                       get_u0_p, check_eqs_u0, CommonSolve

const MTK = ModelingToolkit

function contains_variable(x, wrt)
    any(y -> occursin(y, x), wrt)
end

"""
$(TYPEDSIGNATURES)

Check if `x` is polynomial with respect to the variables in `wrt`.
"""
function is_polynomial(x, wrt)
    x = unwrap(x)
    symbolic_type(x) == NotSymbolic() && return true
    iscall(x) || return true
    contains_variable(x, wrt) || return true
    any(isequal(x), wrt) && return true

    if operation(x) in (*, +, -)
        return all(y -> is_polynomial(y, wrt), arguments(x))
    end
    if operation(x) == (^)
        b, p = arguments(x)
        return is_polynomial(b, wrt) && !contains_variable(p, wrt)
    end
    return false
end

"""
$(TYPEDSIGNATURES)

Convert `expr` from a symbolics expression to one that uses `HomotopyContinuation.ModelKit`.
"""
function symbolics_to_hc(expr)
    if iscall(expr)
        if operation(expr) == getindex
            args = arguments(expr)
            return ModelKit.Variable(getname(args[1]), args[2:end]...)
        else
            return operation(expr)(symbolics_to_hc.(arguments(expr))...)
        end
    elseif symbolic_type(expr) == NotSymbolic()
        return expr
    else
        return ModelKit.Variable(getname(expr))
    end
end

"""
$(TYPEDEF)

A subtype of `HomotopyContinuation.AbstractSystem` used to solve `HomotopyContinuationProblem`s.
"""
struct MTKHomotopySystem{F, P, J, V} <: HomotopyContinuation.AbstractSystem
    """
    The generated function for the residual of the polynomial system. In-place.
    """
    f::F
    """
    The parameter object.
    """
    p::P
    """
    The generated function for the jacobian of the polynomial system. In-place.
    """
    jac::J
    """
    The `HomotopyContinuation.ModelKit.Variable` representation of the unknowns of
    the system.
    """
    vars::V
    """
    The number of polynomials in the system. Must also be equal to `length(vars)`.
    """
    nexprs::Int
end

Base.size(sys::MTKHomotopySystem) = (sys.nexprs, length(sys.vars))
ModelKit.variables(sys::MTKHomotopySystem) = sys.vars

function (sys::MTKHomotopySystem)(x, p = nothing)
    sys.f(x, sys.p)
end

function ModelKit.evaluate!(u, sys::MTKHomotopySystem, x, p = nothing)
    sys.f(u, x, sys.p)
end

function ModelKit.evaluate_and_jacobian!(u, U, sys::MTKHomotopySystem, x, p = nothing)
    sys.f(u, x, sys.p)
    sys.jac(U, x, sys.p)
end

SymbolicIndexingInterface.parameter_values(s::MTKHomotopySystem) = s.p

"""
    $(TYPEDSIGNATURES)

Create a `HomotopyContinuationProblem` from a `NonlinearSystem` with polynomial equations.
The problem will be solved by HomotopyContinuation.jl. The resultant `NonlinearSolution`
will contain the polynomial root closest to the point specified by `u0map` (if real roots
exist for the system).
"""
function MTK.HomotopyContinuationProblem(
        sys::NonlinearSystem, u0map, parammap = nothing; eval_expression = false,
        eval_module = ModelingToolkit, kwargs...)
    if !iscomplete(sys)
        error("A completed `NonlinearSystem` is required. Call `complete` or `structural_simplify` on the system before creating a `HomotopyContinuationProblem`")
    end

    dvs = unknowns(sys)
    eqs = equations(sys)

    for eq in eqs
        if !is_polynomial(eq.lhs, dvs) || !is_polynomial(eq.rhs, dvs)
            error("Equation $eq is not a polynomial in the unknowns")
        end
    end

    nlfn, u0, p = MTK.process_SciMLProblem(NonlinearFunction{true}, sys, u0map, parammap;
        jac = true, eval_expression, eval_module)

    hvars = symbolics_to_hc.(dvs)
    mtkhsys = MTKHomotopySystem(nlfn.f, p, nlfn.jac, hvars, length(eqs))

    obsfn = MTK.ObservedFunctionCache(sys; eval_expression, eval_module)

    return MTK.HomotopyContinuationProblem(u0, mtkhsys, sys, obsfn)
end

"""
$(TYPEDSIGNATURES)

Solve a `HomotopyContinuationProblem`. Ignores the algorithm passed to it, and always
uses `HomotopyContinuation.jl`. All keyword arguments are forwarded to
`HomotopyContinuation.solve`. The original solution as returned by `HomotopyContinuation.jl`
will be available in the `.original` field of the returned `NonlinearSolution`.

All keyword arguments have their default values in HomotopyContinuation.jl, except
`show_progress` which defaults to `false`.
"""
function CommonSolve.solve(prob::MTK.HomotopyContinuationProblem,
        alg = nothing; show_progress = false, kwargs...)
    sol = HomotopyContinuation.solve(
        prob.homotopy_continuation_system; show_progress, kwargs...)
    realsols = HomotopyContinuation.results(sol; only_real = true)
    if isempty(realsols)
        u = state_values(prob)
        resid = prob.homotopy_continuation_system(u)
        retcode = SciMLBase.ReturnCode.ConvergenceFailure
    else
        distance, idx = findmin(realsols) do result
            norm(result.solution - state_values(prob))
        end
        u = real.(realsols[idx].solution)
        resid = prob.homotopy_continuation_system(u)
        retcode = SciMLBase.ReturnCode.Success
    end

    return SciMLBase.build_solution(
        prob, :HomotopyContinuation, u, resid; retcode, original = sol)
end

end
