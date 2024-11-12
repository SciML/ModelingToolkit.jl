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

function contains_variable(x, wrt)
    any(y -> occursin(y, x), wrt)
end

"""
Possible reasons why a term is not polynomial
"""
MTK.EnumX.@enumx NonPolynomialReason begin
    NonIntegerExponent
    ExponentContainsUnknowns
    BaseNotPolynomial
    UnrecognizedOperation
end

function display_reason(reason::NonPolynomialReason.T, sym)
    if reason == NonPolynomialReason.NonIntegerExponent
        pow = arguments(sym)[2]
        "In $sym: Exponent $pow is not an integer"
    elseif reason == NonPolynomialReason.ExponentContainsUnknowns
        pow = arguments(sym)[2]
        "In $sym: Exponent $pow contains unknowns of the system"
    elseif reason == NonPolynomialReason.BaseNotPolynomial
        base = arguments(sym)[1]
        "In $sym: Base $base is not a polynomial in the unknowns"
    elseif reason == NonPolynomialReason.UnrecognizedOperation
        op = operation(sym)
        """
        In $sym: Operation $op is not recognized. Allowed polynomial operations are \
        `*, /, +, -, ^`.
        """
    else
        error("This should never happen. Please open an issue in ModelingToolkit.jl.")
    end
end

mutable struct PolynomialData
    non_polynomial_terms::Vector{BasicSymbolic}
    reasons::Vector{NonPolynomialReason.T}
    has_parametric_exponent::Bool
end

PolynomialData() = PolynomialData(BasicSymbolic[], NonPolynomialReason.T[], false)

struct NotPolynomialError <: Exception
    eq::Equation
    data::PolynomialData
end

function Base.showerror(io::IO, err::NotPolynomialError)
    println(io,
        "Equation $(err.eq) is not a polynomial in the unknowns for the following reasons:")
    for (term, reason) in zip(err.data.non_polynomial_terms, err.data.reasons)
        println(io, display_reason(reason, term))
    end
end

function is_polynomial!(data, y, wrt)
    process_polynomial!(data, y, wrt)
    isempty(data.reasons)
end

"""
$(TYPEDSIGNATURES)

Return information about the polynmial `x` with respect to variables in `wrt`,
writing said information to `data`.
"""
function process_polynomial!(data::PolynomialData, x, wrt)
    x = unwrap(x)
    symbolic_type(x) == NotSymbolic() && return true
    iscall(x) || return true
    contains_variable(x, wrt) || return true
    any(isequal(x), wrt) && return true

    if operation(x) in (*, +, -, /)
        return all(y -> is_polynomial!(data, y, wrt), arguments(x))
    end
    if operation(x) == (^)
        b, p = arguments(x)
        is_pow_integer = symtype(p) <: Integer
        if !is_pow_integer
            push!(data.non_polynomial_terms, x)
            push!(data.reasons, NonPolynomialReason.NonIntegerExponent)
        end
        if symbolic_type(p) != NotSymbolic()
            data.has_parametric_exponent = true
        end

        exponent_has_unknowns = contains_variable(p, wrt)
        if exponent_has_unknowns
            push!(data.non_polynomial_terms, x)
            push!(data.reasons, NonPolynomialReason.ExponentContainsUnknowns)
        end
        base_polynomial = is_polynomial!(data, b, wrt)
        if !base_polynomial
            push!(data.non_polynomial_terms, x)
            push!(data.reasons, NonPolynomialReason.BaseNotPolynomial)
        end
        return base_polynomial && !exponent_has_unknowns && is_pow_integer
    end
    push!(data.non_polynomial_terms, x)
    push!(data.reasons, NonPolynomialReason.UnrecognizedOperation)
    return false
end

"""
$(TYPEDSIGNATURES)

Given a `x`, a polynomial in variables in `wrt` which may contain rational functions,
express `x` as a single rational function with polynomial `num` and denominator `den`.
Return `(num, den)`.
"""
function handle_rational_polynomials(x, wrt)
    x = unwrap(x)
    symbolic_type(x) == NotSymbolic() && return x, 1
    iscall(x) || return x, 1
    contains_variable(x, wrt) || return x, 1
    any(isequal(x), wrt) && return x, 1

    # simplify_fractions cancels out some common factors
    # and expands (a / b)^c to a^c / b^c, so we only need
    # to handle these cases
    x = simplify_fractions(x)
    op = operation(x)
    args = arguments(x)

    if op == /
        # numerator and denominator are trivial
        num, den = args
        # but also search for rational functions in numerator
        n, d = handle_rational_polynomials(num, wrt)
        num, den = n, den * d
    elseif op == +
        num = 0
        den = 1

        # we don't need to do common denominator
        # because we don't care about cases where denominator
        # is zero. The expression is zero when all the numerators
        # are zero.
        for arg in args
            n, d = handle_rational_polynomials(arg, wrt)
            num += n
            den *= d
        end
    else
        return x, 1
    end
    # if the denominator isn't a polynomial in `wrt`, better to not include it
    # to reduce the size of the gcd polynomial
    if !contains_variable(den, wrt)
        return num / den, 1
    end
    return num, den
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

Keyword arguments:
- `eval_expression`: Whether to `eval` the generated functions or use a `RuntimeGeneratedFunction`.
- `eval_module`: The module to use for `eval`/`@RuntimeGeneratedFunction`
- `warn_parametric_exponent`: Whether to warn if the system contains a parametric
  exponent preventing the homotopy from being cached.

All other keyword arguments are forwarded to `HomotopyContinuation.solver_startsystems`.
"""
function MTK.HomotopyContinuationProblem(
        sys::NonlinearSystem, u0map, parammap = nothing; eval_expression = false,
        eval_module = ModelingToolkit, warn_parametric_exponent = true, kwargs...)
    if !iscomplete(sys)
        error("A completed `NonlinearSystem` is required. Call `complete` or `structural_simplify` on the system before creating a `HomotopyContinuationProblem`")
    end

    dvs = unknowns(sys)
    # we need to consider `full_equations` because observed also should be
    # polynomials (if used in equations) and we don't know if observed is used
    # in denominator.
    # This is not the most efficient, and would be improved significantly with
    # CSE/hashconsing.
    eqs = full_equations(sys)

    denoms = []
    has_parametric_exponents = false
    eqs2 = map(eqs) do eq
        data = PolynomialData()
        process_polynomial!(data, eq.lhs, dvs)
        process_polynomial!(data, eq.rhs, dvs)
        has_parametric_exponents |= data.has_parametric_exponent
        if !isempty(data.non_polynomial_terms)
            throw(NotPolynomialError(eq, data))
        end
        num, den = handle_rational_polynomials(eq.rhs - eq.lhs, dvs)

        # make factors different elements, otherwise the nonzero factors artificially
        # inflate the error of the zero factor.
        if iscall(den) && operation(den) == *
            for arg in arguments(den)
                # ignore constant factors
                symbolic_type(arg) == NotSymbolic() && continue
                push!(denoms, abs(arg))
            end
        elseif symbolic_type(den) != NotSymbolic()
            push!(denoms, abs(den))
        end
        return 0 ~ num
    end

    sys2 = MTK.@set sys.eqs = eqs2
    # remove observed equations to avoid adding them in codegen
    MTK.@set! sys2.observed = Equation[]
    MTK.@set! sys2.substitutions = nothing

    nlfn, u0, p = MTK.process_SciMLProblem(NonlinearFunction{true}, sys2, u0map, parammap;
        jac = true, eval_expression, eval_module)

    denominator = MTK.build_explicit_observed_function(sys, denoms)

    hvars = symbolics_to_hc.(dvs)
    mtkhsys = MTKHomotopySystem(nlfn.f, p, nlfn.jac, hvars, length(eqs))

    obsfn = MTK.ObservedFunctionCache(sys; eval_expression, eval_module)

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
        u0, mtkhsys, denominator, sys, obsfn, solver_and_starts)
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
    else
        T = eltype(state_values(prob))
        distance, idx = findmin(realsols) do result
            if any(<=(denominator_abstol),
                prob.denominator(real.(result.solution), parameter_values(prob)))
                return T(Inf)
            end
            norm(result.solution - state_values(prob))
        end
        # all roots cause denominator to be zero
        if isinf(distance)
            u = state_values(prob)
            retcode = SciMLBase.ReturnCode.Infeasible
        else
            u = real.(realsols[idx].solution)
            retcode = SciMLBase.ReturnCode.Success
        end
    end
    resid = prob.homotopy_continuation_system(u)

    return SciMLBase.build_solution(
        prob, :HomotopyContinuation, u, resid; retcode, original = sol)
end

end
