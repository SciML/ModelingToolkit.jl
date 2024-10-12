module MTKHomotopyContinuationExt

using ModelingToolkit
using ModelingToolkit.SciMLBase
using ModelingToolkit.Symbolics: unwrap
using ModelingToolkit.SymbolicIndexingInterface
using HomotopyContinuation
using ModelingToolkit: iscomplete, parameters, has_index_cache, get_index_cache, get_u0,
                       get_u0_p, check_eqs_u0, CommonSolve

const MTK = ModelingToolkit

function contains_variable(x, wrt)
    any(isequal(x), wrt) && return true
    iscall(x) || return false
    return any(y -> contains_variable(y, wrt), arguments(x))
end

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

struct MTKHomotopySystem{F, P, J, V} <: HomotopyContinuation.AbstractSystem
    f::F
    p::P
    jac::J
    vars::V
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

function MTK.HomotopyContinuationProblem(
        sys::NonlinearSystem, u0map, parammap; compile = :all, eval_expression = false, eval_module = ModelingToolkit, kwargs...)
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

function CommonSolve.solve(prob::MTK.HomotopyContinuationProblem; kwargs...)
    sol = HomotopyContinuation.solve(prob.homotopy_continuation_system; kwargs...)
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
