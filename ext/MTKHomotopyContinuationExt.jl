module MTKHomotopyContinuationExt

using ModelingToolkit
using ModelingToolkit.SciMLBase
using ModelingToolkit.Symbolics: unwrap, symtype, BasicSymbolic, simplify_fractions
using ModelingToolkit.SymbolicIndexingInterface
using ModelingToolkit.DocStringExtensions
using HomotopyContinuation
using ModelingToolkit: iscomplete, parameters, has_index_cache, get_index_cache, get_u0,
                       get_u0_p, check_eqs_u0, CommonSolve

const MTK = ModelingToolkit

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

Keyword arguments:
- `eval_expression`: Whether to `eval` the generated functions or use a `RuntimeGeneratedFunction`.
- `eval_module`: The module to use for `eval`/`@RuntimeGeneratedFunction`
- `warn_parametric_exponent`: Whether to warn if the system contains a parametric
  exponent preventing the homotopy from being cached.

All other keyword arguments are forwarded to `HomotopyContinuation.solver_startsystems`.
"""
function MTK.HomotopyContinuationProblem(
        sys::NonlinearSystem, u0map, parammap = nothing; kwargs...)
    prob = MTK._safe_HomotopyContinuationProblem(sys, u0map, parammap; kwargs...)
    prob isa MTK.HomotopyContinuationProblem || throw(prob)
    return prob
end

function MTK._safe_HomotopyContinuationProblem(sys, u0map, parammap = nothing;
        fraction_cancel_fn = SymbolicUtils.simplify_fractions, kwargs...)
    if !iscomplete(sys)
        error("A completed `NonlinearSystem` is required. Call `complete` or `structural_simplify` on the system before creating a `HomotopyContinuationProblem`")
    end
    transformation = MTK.PolynomialTransformation(sys)
    if transformation isa MTK.NotPolynomialError
        return transformation
    end
    result = MTK.transform_system(sys, transformation; fraction_cancel_fn)
    if result isa MTK.NotPolynomialError
        return result
    end
    MTK.HomotopyContinuationProblem(sys, transformation, result, u0map, parammap; kwargs...)
end

function MTK.HomotopyContinuationProblem(
        sys::MTK.NonlinearSystem, transformation::MTK.PolynomialTransformation,
        result::MTK.PolynomialTransformationResult, u0map,
        parammap = nothing; eval_expression = false,
        eval_module = ModelingToolkit, warn_parametric_exponent = true, kwargs...)
    sys2 = result.sys
    denoms = result.denominators
    polydata = transformation.polydata
    new_dvs = transformation.new_dvs
    all_solutions = transformation.all_solutions

    _, u0, p = MTK.process_SciMLProblem(
        MTK.EmptySciMLFunction, sys, u0map, parammap; eval_expression, eval_module)
    nlfn = NonlinearFunction{true}(sys2; jac = true, eval_expression, eval_module)

    denominator = MTK.build_explicit_observed_function(sys2, denoms)
    unpack_solution = MTK.build_explicit_observed_function(sys2, all_solutions)

    hvars = symbolics_to_hc.(new_dvs)
    mtkhsys = MTKHomotopySystem(nlfn.f, p, nlfn.jac, hvars, length(new_dvs))

    obsfn = MTK.ObservedFunctionCache(sys; eval_expression, eval_module)

    has_parametric_exponents = any(d -> d.has_parametric_exponent, polydata)
    if has_parametric_exponents
        if warn_parametric_exponent
            @warn """
            The system has parametric exponents, preventing caching of the homotopy. \
            This will cause `solve` to be slower. Pass `warn_parametric_exponent \
            = false` to turn off this warning
            """
        end
        solver_and_starts = nothing
    else
        solver_and_starts = HomotopyContinuation.solver_startsolutions(mtkhsys; kwargs...)
    end
    return MTK.HomotopyContinuationProblem(
        u0, mtkhsys, denominator, sys, obsfn, solver_and_starts, unpack_solution)
end

"""
$(TYPEDSIGNATURES)

Solve a `HomotopyContinuationProblem`. Ignores the algorithm passed to it, and always
uses `HomotopyContinuation.jl`. The original solution as returned by
`HomotopyContinuation.jl` will be available in the `.original` field of the returned
`NonlinearSolution`.

All keyword arguments except the ones listed below are forwarded to
`HomotopyContinuation.solve`. Note that the solver and start solutions are precomputed,
and only keyword arguments related to the solve process are valid. All keyword
arguments have their default values in HomotopyContinuation.jl, except `show_progress`
which defaults to `false`.

Extra keyword arguments:
- `denominator_abstol`: In case `prob` is solving a rational function, roots which cause
  the denominator to be below `denominator_abstol` will be discarded.
"""
function CommonSolve.solve(prob::MTK.HomotopyContinuationProblem,
        alg = nothing; show_progress = false, denominator_abstol = 1e-7, kwargs...)
    if prob.solver_and_starts === nothing
        sol = HomotopyContinuation.solve(
            prob.homotopy_continuation_system; show_progress, kwargs...)
    else
        solver, starts = prob.solver_and_starts
        sol = HomotopyContinuation.solve(solver, starts; show_progress, kwargs...)
    end
    realsols = HomotopyContinuation.results(sol; only_real = true)
    if isempty(realsols)
        u = state_values(prob)
        retcode = SciMLBase.ReturnCode.ConvergenceFailure
        resid = prob.homotopy_continuation_system(u)
    else
        T = eltype(state_values(prob))
        distance = T(Inf)
        u = state_values(prob)
        resid = nothing
        for result in realsols
            if any(<=(denominator_abstol),
                prob.denominator(real.(result.solution), parameter_values(prob)))
                continue
            end
            for truesol in prob.unpack_solution(result.solution, parameter_values(prob))
                dist = norm(truesol - state_values(prob))
                if dist < distance
                    distance = dist
                    u = T.(real.(truesol))
                    resid = T.(real.(prob.homotopy_continuation_system(result.solution)))
                end
            end
        end
        # all roots cause denominator to be zero
        if isinf(distance)
            u = state_values(prob)
            resid = prob.homotopy_continuation_system(u)
            retcode = SciMLBase.ReturnCode.Infeasible
        else
            retcode = SciMLBase.ReturnCode.Success
        end
    end

    return SciMLBase.build_solution(
        prob, :HomotopyContinuation, u, resid; retcode, original = sol)
end

end
