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

"""
    $(TYPEDSIGNATURES)

Utility function for `safe_HomotopyContinuationProblem`, implemented in the extension.
"""
function _safe_HomotopyContinuationProblem end

"""
    $(TYPEDSIGNATURES)

Return a `HomotopyContinuationProblem` if the extension is loaded and the system is
polynomial. If the extension is not loaded, return `nothing`. If the system is not
polynomial, return the appropriate `NotPolynomialError`.
"""
function safe_HomotopyContinuationProblem(sys::NonlinearSystem, args...; kwargs...)
    if Base.get_extension(ModelingToolkit, :MTKHomotopyContinuationExt) === nothing
        return nothing
    end
    return _safe_HomotopyContinuationProblem(sys, args...; kwargs...)
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

function contains_variable(x, wrt)
    any(y -> occursin(y, x), wrt)
end

"""
Possible reasons why a term is not polynomial
"""
EnumX.@enumx NonPolynomialReason begin
    """
    Exponent of an expression involving unknowns is not an integer.
    """
    NonIntegerExponent
    """
    Exponent is an expression containing unknowns.
    """
    ExponentContainsUnknowns
    """
    The base of an exponent is not a polynomial in the unknowns.
    """
    BaseNotPolynomial
    """
    An expression involves a non-polynomial operation involving unknowns.
    """
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

"""
    $(TYPEDEF)

Information about an expression about its polynomial nature.
"""
mutable struct PolynomialData
    """
    A list of all non-polynomial terms in the expression.
    """
    non_polynomial_terms::Vector{BasicSymbolic}
    """
    Corresponding to `non_polynomial_terms`, a list of reasons why they are
    not polynomial.
    """
    reasons::Vector{NonPolynomialReason.T}
    """
    Whether the polynomial contains parametric exponents of unknowns.
    """
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
    $(TYPEDEF)

Information about how an unknown in the system is substituted for a non-polynomial
expression to turn the system into a polynomial. Used in `PolynomialTransformation`.
"""
struct PolynomialTransformationData
    """
    The new variable to use as an unknown of the transformed system.
    """
    new_var::BasicSymbolic
    """
    The non-polynomial expression being substituted.
    """
    term::BasicSymbolic
    """
    A vector of expressions corresponding to the solutions of
    the non-polynomial expression `term` in terms of the new unknown `new_var`,
    used to backsolve for the original unknown of the system.
    """
    inv_term::Vector{BasicSymbolic}
end

"""
    $(TYPEDEF)

Information representing how to transform a `NonlinearSystem` into a polynomial
system.
"""
struct PolynomialTransformation
    """
    Substitutions mapping non-polynomial terms to temporary unknowns. The system
    is a polynomial in the new unknowns. Currently, each non-polynomial term is a
    function of a single unknown of the original system.
    """
    substitution_rules::Dict{BasicSymbolic, BasicSymbolic}
    """
    A vector of expressions involving unknowns of the transformed system, mapping
    back to solutions of the original system.
    """
    all_solutions::Vector{Vector{BasicSymbolic}}
    """
    The new unknowns of the transformed system.
    """
    new_dvs::Vector{BasicSymbolic}
    """
    The polynomial data for each equation.
    """
    polydata::Vector{PolynomialData}
end

function PolynomialTransformation(sys::NonlinearSystem)
    # we need to consider `full_equations` because observed also should be
    # polynomials (if used in equations) and we don't know if observed is used
    # in denominator.
    # This is not the most efficient, and would be improved significantly with
    # CSE/hashconsing.
    eqs = full_equations(sys)
    dvs = unknowns(sys)

    # Collect polynomial information about all equations
    polydata = map(eqs) do eq
        data = PolynomialData()
        process_polynomial!(data, eq.lhs, dvs)
        process_polynomial!(data, eq.rhs, dvs)
        data
    end

    # Get all unique non-polynomial terms
    # NOTE:
    # Is there a better way to check for uniqueness? `simplify` is relatively slow
    # (maybe use the threaded version?) and `expand` can blow up expression size.
    # Could metatheory help?
    all_non_poly_terms = mapreduce(
        d -> d.non_polynomial_terms, vcat, polydata; init = BasicSymbolic[])
    unique!(all_non_poly_terms)

    # each variable can only be replaced by one non-polynomial expression involving
    # that variable. Keep track of this mapping.
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
        if any(x -> iscall(x) && operation(x) == Symbolics.RootsOf, invterm)
            transformation_err = NemoNotLoaded()
            is_poly = false
            break
        end
        var_to_nonpoly[var] = PolynomialTransformationData(new_var, t, invterm)
    end

    # return the error instead of throwing it, so the user can choose what to do
    # without having to catch the exception
    if !is_poly
        return NotPolynomialError(transformation_err, eqs, polydata)
    end

    subrules = Dict{BasicSymbolic, BasicSymbolic}()
    # corresponding to each unknown in `dvs`, the list of its possible solutions
    # in terms of the new unknown.
    combinations = Vector{BasicSymbolic}[]
    new_dvs = BasicSymbolic[]
    for x in dvs
        if haskey(var_to_nonpoly, x)
            _data = var_to_nonpoly[x]
            # map term to new unknown
            subrules[_data.term] = _data.new_var
            push!(combinations, _data.inv_term)
            push!(new_dvs, _data.new_var)
        else
            push!(combinations, BasicSymbolic[x])
            push!(new_dvs, x)
        end
    end
    all_solutions = vec(collect.(collect(Iterators.product(combinations...))))
    return PolynomialTransformation(subrules, all_solutions, new_dvs, polydata)
end

"""
    $(TYPEDEF)

A struct containing the result of transforming a system into a polynomial system
using the appropriate `PolynomialTransformation`. Also contains the denominators
in the equations, to rule out invalid roots.
"""
struct PolynomialTransformationResult
    sys::NonlinearSystem
    denominators::Vector{BasicSymbolic}
end

"""
    $(TYPEDSIGNATURES)

Transform the system `sys` with `transformation` and return a
`PolynomialTransformationResult`, or a `NotPolynomialError` if the system cannot
be transformed.
"""
function transform_system(sys::NonlinearSystem, transformation::PolynomialTransformation;
        fraction_cancel_fn = simplify_fractions)
    subrules = transformation.substitution_rules
    dvs = unknowns(sys)
    eqs = full_equations(sys)
    polydata = transformation.polydata
    new_dvs = transformation.new_dvs
    all_solutions = transformation.all_solutions

    eqs2 = Equation[]
    denoms = BasicSymbolic[]
    for eq in eqs
        t = eq.rhs - eq.lhs
        t = Symbolics.fixpoint_sub(t, subrules; maxiters = length(dvs))
        # the substituted variable occurs outside the substituted term
        poly_and_nonpoly = map(dvs) do x
            all(!isequal(x), new_dvs) && occursin(x, t)
        end
        if any(poly_and_nonpoly)
            return NotPolynomialError(
                VariablesAsPolyAndNonPoly(dvs[poly_and_nonpoly]), eqs, polydata)
        end
        num, den = handle_rational_polynomials(t, new_dvs; fraction_cancel_fn)
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
        push!(eqs2, 0 ~ num)
    end

    sys2 = @set sys.eqs = eqs2
    @set! sys2.unknowns = new_dvs
    # remove observed equations to avoid adding them in codegen
    @set! sys2.observed = Equation[]
    @set! sys2.substitutions = nothing
    return PolynomialTransformationResult(sys2, denoms)
end

"""
$(TYPEDSIGNATURES)

Given a `x`, a polynomial in variables in `wrt` which may contain rational functions,
express `x` as a single rational function with polynomial `num` and denominator `den`.
Return `(num, den)`.

Keyword arguments:
- `fraction_cancel_fn`: A function which takes a fraction (`operation(expr) == /`) and returns
  a simplified symbolic quantity with common factors in the numerator and denominator are
  cancelled. Defaults to `SymbolicUtils.simplify_fractions`, but can be changed to
  `nothing` to improve performance on large polynomials at the cost of avoiding non-trivial
  cancellation.
"""
function handle_rational_polynomials(x, wrt; fraction_cancel_fn = simplify_fractions)
    x = unwrap(x)
    symbolic_type(x) == NotSymbolic() && return x, 1
    iscall(x) || return x, 1
    contains_variable(x, wrt) || return x, 1
    any(isequal(x), wrt) && return x, 1

    op = operation(x)
    args = arguments(x)

    if op == /
        # numerator and denominator are trivial
        num, den = args
        n1, d1 = handle_rational_polynomials(num, wrt; fraction_cancel_fn)
        n2, d2 = handle_rational_polynomials(den, wrt; fraction_cancel_fn)
        num, den = n1 * d2, d1 * n2
    elseif (op == +) || (op == -)
        num = 0
        den = 1
        if op == -
            args[2] = -args[2]
        end
        for arg in args
            n, d = handle_rational_polynomials(arg, wrt; fraction_cancel_fn)
            num = num * d + n * den
            den *= d
        end
    elseif op == ^
        base, pow = args
        num, den = handle_rational_polynomials(base, wrt; fraction_cancel_fn)
        num ^= pow
        den ^= pow
    elseif op == *
        num = 1
        den = 1
        for arg in args
            n, d = handle_rational_polynomials(arg, wrt; fraction_cancel_fn)
            num *= n
            den *= d
        end
    else
        error("Unhandled operation in `handle_rational_polynomials`. This should never happen. Please open an issue in ModelingToolkit.jl with an MWE.")
    end

    if fraction_cancel_fn !== nothing
        expr = fraction_cancel_fn(num / den)
        if iscall(expr) && operation(expr) == /
            num, den = arguments(expr)
        else
            num, den = expr, 1
        end
    end

    # if the denominator isn't a polynomial in `wrt`, better to not include it
    # to reduce the size of the gcd polynomial
    if !contains_variable(den, wrt)
        return num / den, 1
    end
    return num, den
end
