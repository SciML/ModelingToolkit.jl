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

abstract type PolynomialTransformationError <: Exception end

struct MultivarTerm <: PolynomialTransformationError
    term::Any
    vars::Any
end

function Base.showerror(io::IO, err::MultivarTerm)
    println(io,
        "Cannot convert system to polynomial: Found term $(err.term) which is a function of multiple unknowns $(err.vars).")
end

struct MultipleTermsOfSameVar <: PolynomialTransformationError
    terms::Any
    var::Any
end

function Base.showerror(io::IO, err::MultipleTermsOfSameVar)
    println(io,
        "Cannot convert system to polynomial: Found multiple non-polynomial terms $(err.terms) involving the same unknown $(err.var).")
end

struct SymbolicSolveFailure <: PolynomialTransformationError
    term::Any
    var::Any
end

function Base.showerror(io::IO, err::SymbolicSolveFailure)
    println(io,
        "Cannot convert system to polynomial: Unable to symbolically solve $(err.term) for $(err.var).")
end

struct NemoNotLoaded <: PolynomialTransformationError end

function Base.showerror(io::IO, err::NemoNotLoaded)
    println(io,
        "ModelingToolkit may be able to solve this system as a polynomial system if `Nemo` is loaded. Run `import Nemo` and try again.")
end

struct VariablesAsPolyAndNonPoly <: PolynomialTransformationError
    vars::Any
end

function Base.showerror(io::IO, err::VariablesAsPolyAndNonPoly)
    println(io,
        "Cannot convert convert system to polynomial: Variables $(err.vars) occur in both polynomial and non-polynomial terms in the system.")
end

struct NotPolynomialError <: Exception
    transformation_err::Union{PolynomialTransformationError, Nothing}
    eq::Vector{Equation}
    data::Vector{PolynomialData}
end

function Base.showerror(io::IO, err::NotPolynomialError)
    if err.transformation_err !== nothing
        Base.showerror(io, err.transformation_err)
    end
    for (eq, data) in zip(err.eq, err.data)
        if isempty(data.non_polynomial_terms)
            continue
        end
        println(io,
            "Equation $(eq) is not a polynomial in the unknowns for the following reasons:")
        for (term, reason) in zip(data.non_polynomial_terms, data.reasons)
            println(io, display_reason(reason, term))
        end
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
        # `map` because `all` will early exit, but we want to search
        # through everything to get all the non-polynomial terms
        return all(map(y -> is_polynomial!(data, y, wrt), arguments(x)))
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

struct PolynomialTransformationData
    new_var::BasicSymbolic
    term::BasicSymbolic
    inv_term::Vector
end

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

    polydata = map(eqs) do eq
        data = PolynomialData()
        process_polynomial!(data, eq.lhs, dvs)
        process_polynomial!(data, eq.rhs, dvs)
        data
    end

    has_parametric_exponents = any(d -> d.has_parametric_exponent, polydata)

    all_non_poly_terms = mapreduce(d -> d.non_polynomial_terms, vcat, polydata)
    unique!(all_non_poly_terms)

    var_to_nonpoly = Dict{BasicSymbolic, PolynomialTransformationData}()

    is_poly = true
    transformation_err = nothing
    for t in all_non_poly_terms
        # if the term involves multiple unknowns, we can't invert it
        dvs_in_term = map(x -> occursin(x, t), dvs)
        if count(dvs_in_term) > 1
            transformation_err = MultivarTerm(t, dvs[dvs_in_term])
            is_poly = false
            break
        end
        # we already have a substitution solving for `var`
        var = dvs[findfirst(dvs_in_term)]
        if haskey(var_to_nonpoly, var) && !isequal(var_to_nonpoly[var].term, t)
            transformation_err = MultipleTermsOfSameVar([t, var_to_nonpoly[var].term], var)
            is_poly = false
            break
        end
        # we want to solve `term - new_var` for `var`
        new_var = gensym(Symbol(var))
        new_var = unwrap(only(@variables $new_var))
        invterm = Symbolics.ia_solve(
            t - new_var, var; complex_roots = false, periodic_roots = false, warns = false)
        # if we can't invert it, quit
        if invterm === nothing || isempty(invterm)
            transformation_err = SymbolicSolveFailure(t, var)
            is_poly = false
            break
        end
        # `ia_solve` returns lazy terms i.e. `asin(1.0)` instead of `pi/2`
        # this just evaluates the constant expressions
        invterm = Symbolics.substitute.(invterm, (Dict(),))
        # RootsOf implies Symbolics couldn't solve the inner polynomial because
        # `Nemo` wasn't loaded.
        if any(x -> MTK.iscall(x) && MTK.operation(x) == Symbolics.RootsOf, invterm)
            transformation_err = NemoNotLoaded()
            is_poly = false
            break
        end
        var_to_nonpoly[var] = PolynomialTransformationData(new_var, t, invterm)
    end

    if !is_poly
        throw(NotPolynomialError(transformation_err, eqs, polydata))
    end

    subrules = Dict()
    combinations = Vector[]
    new_dvs = []
    for x in dvs
        if haskey(var_to_nonpoly, x)
            _data = var_to_nonpoly[x]
            subrules[_data.term] = _data.new_var
            push!(combinations, _data.inv_term)
            push!(new_dvs, _data.new_var)
        else
            push!(combinations, [x])
            push!(new_dvs, x)
        end
    end
    all_solutions = collect.(collect(Iterators.product(combinations...)))

    denoms = []
    eqs2 = map(eqs) do eq
        t = eq.rhs - eq.lhs
        t = Symbolics.fixpoint_sub(t, subrules; maxiters = length(dvs))
        # the substituted variable occurs outside the substituted term
        poly_and_nonpoly = map(dvs) do x
            haskey(var_to_nonpoly, x) && occursin(x, t)
        end
        if any(poly_and_nonpoly)
            throw(NotPolynomialError(
                VariablesAsPolyAndNonPoly(dvs[poly_and_nonpoly]), eqs, polydata))
        end

        num, den = handle_rational_polynomials(t, new_dvs)
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
    MTK.@set! sys2.unknowns = new_dvs
    # remove observed equations to avoid adding them in codegen
    MTK.@set! sys2.observed = Equation[]
    MTK.@set! sys2.substitutions = nothing

    _, u0, p = MTK.process_SciMLProblem(
        MTK.EmptySciMLFunction, sys, u0map, parammap; eval_expression, eval_module)
    nlfn = NonlinearFunction{true}(sys2; jac = true, eval_expression, eval_module)

    denominator = MTK.build_explicit_observed_function(sys2, denoms)
    unpack_solution = MTK.build_explicit_observed_function(sys2, all_solutions)

    hvars = symbolics_to_hc.(new_dvs)
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
