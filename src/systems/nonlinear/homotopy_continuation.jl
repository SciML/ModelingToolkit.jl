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

abstract type PolynomialTransformationError <: Exception end

struct UnmatchedUnknowns <: PolynomialTransformationError
    unmatched::Vector{BasicSymbolic}
end

function Base.showerror(io::IO, err::UnmatchedUnknowns)
    println(io,
        "Cannot convert system to polynomial: could not find terms to solve for unknowns $(err.unmatched).")
end

struct UnmatchedTerms <: PolynomialTransformationError
    unmatched::Vector{BasicSymbolic}
end

function Base.showerror(io::IO, err::UnmatchedTerms)
    println(io, "Cannot convert system to polynomial: too many non-polynomial terms in system. Unmatched terms are $(err.unmatched).")
end

function no_nemo_warning()
    @warn "ModelingToolkit may be able to symbolically solve some non-polynomial terms in this system for all roots if `Nemo` is loaded. Run `import Nemo` and try again to enable this functionality and possibly obtain additional roots."
end

struct NotPolynomialError <: Exception
    transformation_err::Union{PolynomialTransformationError, Nothing}
    eq::Vector{Equation}
    data::Vector{Any} # PolynomialData but then I have to put that struct above this
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

"""
    $(TYPEDEF)

Information about an expression about its polynomial nature.
"""
mutable struct PolynomialData
    """
    A list of the variables in the expression being searched that occur outside of
    any non-polynomial terms.
    """
    solo_terms::Set{BasicSymbolic}
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

PolynomialData() = PolynomialData(Set{BasicSymbolic}(), BasicSymbolic[], NonPolynomialReason.T[], false)

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
    contains_variable(x, wrt) || return true
    if any(isequal(x), wrt)
        push!(data.solo_terms, x)
        return true
    end
    iscall(x) || return true

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

Information for how to solve for unknowns involved in non-symbolically-solvable
non-polynomial terms to turn the system into a polynomial. Used in
`PolynomialTransformation`.
"""
struct NonlinearSolveTransformation
    """
    The system which solves for the unknowns of the parent system.
    """
    sys::NonlinearSystem
    """
    The input variables to this system representing solutions of non-polynomial terms.
    """
    inputvars::Vector{BasicSymbolic}
end

"""
    $(TYPEDEF)

Information representing how to transform a `NonlinearSystem` into a polynomial
system.
"""
struct PolynomialTransformation
    """
    The stages in which to recover the solution in terms of the original unknowns, in
    order.
    """
    solve_stages::Vector{Any}
    """
    The (previous) stages each stage depends on.
    """
    stage_dependencies::Vector{Vector{Int}}
    """
    Mapping from terms to new unknowns they are replaced by. The system is a
    polynomial in the new unknowns.
    """
    substitution_rules::Dict{BasicSymbolic, BasicSymbolic}
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

    all_solo_vars = mapreduce(d -> d.solo_terms, union, polydata; init = Set{BasicSymbolic}())

    # Graph matches variables to candidates for unknowns of the polynomial system that
    # they occur in. These unknowns can be solo variables that appear outside of
    # non-polynomial terms in the system, or non-polynomials.
    graph = BipartiteGraph(length(dvs), 0)
    # all solo variables are candidates for unknowns
    graph_srcs = dvs
    graph_dsts = BasicSymbolic[]
    for (i, var) in enumerate(dvs)
        var in all_solo_vars || continue
        push!(graph_dsts, var)
        vert = add_vertex!(graph, DST)
        add_edge!(graph, i, vert)
    end

    # buffer to prevent reallocations
    dvs_in_term = Set()
    # for efficient queries
    dvs_to_src = Dict(graph_srcs .=> eachindex(graph_srcs))
    # build out graph with other non-polynomial terms
    for t in all_non_poly_terms
        empty!(dvs_in_term)
        vars!(dvs_in_term, t)
        intersect!(dvs_in_term, dvs)
        push!(graph_dsts, t)
        vert = add_vertex!(graph, DST)
        for var in dvs_in_term
            add_edge!(graph, dvs_to_src[var], vert)
        end
    end

    # Match variables to the candidate unknown we'll use to solve for them.
    # This is a poor man's version of `structural_simplify`, but if we create
    # and simplify a `NonlinearSystem` it makes doing symbolic solving more
    # annoying.
    matching = BipartiteGraphs.complete(maximal_matching(graph))
    inv_matching = invview(matching)
    # matching is from destination to source vertices
    unassigned_dsts = filter(i -> matching[i] == unassigned, 𝑠vertices(graph))
    unassigned_srcs = filter(i -> inv_matching[i] == unassigned, 𝑑vertices(graph))

    # return the error instead of throwing it, so the user can choose what to do
    # without having to catch the exception
    if !isempty(unassigned_srcs)
        return NotPolynomialError(UnmatchedUnknowns(graph_srcs[unassigned_srcs]), eqs, polydata)
    end
    if !isempty(unassigned_dsts)
        return NotPolynomialError(UnmatchedTerms(graph_dsts[unassigned_dsts]), eqs, polydata)
    end

    # At this point, the matching is perfect. Find the SCCs so we know
    # which terms to solve for which variables.
    digraph = DiCMOBiGraph{false}(graph, matching)
    var_sccs = Graphs.strongly_connected_components(digraph)
    foreach(sort!, var_sccs)
    # construct a condensation graph of the SCCs so we can topologically sort them
    scc_graph = MatchedCondensationGraph(digraph, var_sccs)
    toporder = topological_sort(scc_graph)
    var_sccs = var_sccs[toporder]
    # get the corresponding terms
    term_sccs = map(var_sccs) do scc
        map(scc) do src
            inv_matching[src]
        end
    end

    # keep track of which previous SCCs each SCC depends on
    dependencies = Vector{Int}[]
    # the method to solve each stage
    solve_stages = []
    # mapping from terms to the new unknowns they are replaced by
    subrules = Dict{BasicSymbolic, BasicSymbolic}()
    # if we've already emitted the no nemo warning
    warned_no_nemo = false
    for (i, (vscc, tscc)) in enumerate(zip(var_sccs, term_sccs))
        # dependencies are simply outneighbors
        push!(dependencies, collect(Graphs.outneighbors(scc_graph, i)))

        # whether the SCC is solvable with a single variable
        single_scc_solvable = length(vscc) == 1
        # for single-variable SCCs, we use `ia_solve`
        if single_scc_solvable
            varidx = vscc[]
            termidx = tscc[]
            var = graph_srcs[varidx]
            t = graph_dsts[termidx]
            # Create a new variable and representing the non-polynomial term...
            new_var = unwrap(similar_variable(var, Symbol(var)))
            # ...and solve for `var` in terms of this new variable.
            invterm = Symbolics.ia_solve(t - new_var, var; complex_roots = false, periodic_roots = false, warns = false)
            # `ia_solve` returns lazy terms i.e. `asin(1.0)` instead of `pi/2`
            # this just evaluates the constant expressions
            invterm = Symbolics.substitute.(invterm, (Dict(),))
            # if `ia_solve` returns `nothing`, the broadcast above turns it into `(nothing,)`
            if invterm === (nothing,) || isempty(invterm)
                # if we can't invert it, quit
                single_scc_solvable = false
            elseif any(x -> iscall(x) && operation(x) == Symbolics.RootsOf, invterm)
                # RootsOf implies Symbolics couldn't solve the inner polynomial because
                # `Nemo` wasn't loaded.
                warned_no_nemo || no_nemo_warning()
                warned_no_nemo = true
                single_scc_solvable = false
            else
                subrules[t] = new_var
                push!(solve_stages, PolynomialTransformationData(new_var, t, invterm))
            end
        end

        # the SCC was solved with a single variable
        single_scc_solvable && continue

        # Solve using a `NonlinearSolve`.
        vars = graph_srcs[vscc]
        ts = graph_dsts[tscc]
        # the new variables are inputs to the system, so they're parameters
        new_vars = map(vars) do var
            toparam.(unwrap.(similar_variable(var, Symbol(var))))
        end
        eqs = collect(0 .~ (ts .- new_vars))
        scc_sys = complete(NonlinearSystem(eqs; name = Symbol(:scc_, i)))
        push!(solve_stages, NonlinearSolveTransformation(scc_sys, new_vars))
        for (new_var, t) in zip(new_vars, ts)
            subrules[t] = new_var
        end
    end

    return PolynomialTransformation(solve_stages, dependencies, subrules, polydata)
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
    """
    The stages in which to recover the solution in terms of the original unknowns, in
    order.
    """
    solve_stages::Vector{Any}
    """
    The (previous) stages each stage depends on.
    """
    stage_dependencies::Vector{Vector{Int}}
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
    new_dvs = collect(values(subrules))

    eqs2 = Equation[]
    denoms = BasicSymbolic[]
    for eq in eqs
        t = eq.rhs - eq.lhs
        t = Symbolics.fixpoint_sub(t, subrules; maxiters = length(dvs))
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
    return PolynomialTransformationResult(sys2, denoms, transformation.solve_stages, transformation.stage_dependencies)
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
