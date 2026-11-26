struct LoggedFunctionException <: Exception
    msg::String
end
struct LoggedFun{F}
    f::F
    args::Any
    error_nonfinite::Bool
end
function LoggedFunctionException(lf::LoggedFun, args, msg)
    return LoggedFunctionException(
        "Function $(lf.f)($(join(lf.args, ", "))) " * msg * " with input" *
            join("\n  " .* string.(lf.args .=> args)) # one line for each "var => val" for readability
    )
end
Base.showerror(io::IO, err::LoggedFunctionException) = print(io, err.msg)
Base.nameof(lf::LoggedFun) = nameof(lf.f)
SymbolicUtils.promote_symtype(f::LoggedFun, Ts::SU.TypeT...) = SU.promote_symtype(f.f, Ts...)
SU.promote_shape(f::LoggedFun, @nospecialize(shs::SU.ShapeT...)) = SU.promote_shape(f.f, shs...)
function (lf::LoggedFun)(args...)
    val = try
        lf.f(args...) # try to call with numerical input, as usual
    catch err
        throw(LoggedFunctionException(lf, args, "errors")) # Julia automatically attaches original error message
    end
    if lf.error_nonfinite && !isfinite(val)
        throw(LoggedFunctionException(lf, args, "output non-finite value $val"))
    end
    return val
end

function logged_fun(f, args...; error_nonfinite = true) # remember to update error_nonfinite in debug_system() docstring
    # Currently we don't really support complex numbers
    return maketerm(SymbolicT, LoggedFun(f, args, error_nonfinite), args, nothing)
end

function debug_sub(eq::Equation, funcs; kw...)
    return debug_sub(eq.lhs, funcs; kw...) ~ debug_sub(eq.rhs, funcs; kw...)
end
function debug_sub(ex, funcs; kw...)
    iscall(ex) || return ex
    f = operation(ex)
    args = map(ex -> debug_sub(ex, funcs; kw...), arguments(ex))
    return f in funcs ? logged_fun(f, args...; kw...) :
        maketerm(typeof(ex), f, args, metadata(ex))
end

"""
    $(TYPEDSIGNATURES)

A function which returns `NaN` if `condition` fails, and `0.0` otherwise.
"""
function _nan_condition(condition)
    # Sparsity tracers represent comparisons without a Boolean value. Their `ifelse`
    # overloads can combine the constant branches without introducing dependencies.
    return ifelse(condition, 0.0, NaN)
end

@register_symbolic _nan_condition(condition::Bool)

"""
    $(TYPEDSIGNATURES)

A function which takes a condition `expr` and returns `NaN` if it is false,
and zero if it is true. In case the condition is false and `log == true`,
`message` will be logged as an `@error`.
"""
function _debug_assertion(expr::Bool, message::String, log::Bool)
    value = _nan_condition(expr)
    isnan(value) || return value
    log && @error message
    return value
end

@register_symbolic _debug_assertion(expr::Bool, message::String, log::Bool)

"""
Boolean parameter added to models returned from `debug_system` to control logging of
assertions.
"""
const ASSERTION_LOG_VARIABLE = only(@parameters __log_assertions_ₘₜₖ::Bool = false)

"""
    $(TYPEDSIGNATURES)

Get a symbolic expression for all the assertions in `sys`. The expression returns `NaN`
if any of the assertions fail, and `0.0` otherwise. If `ASSERTION_LOG_VARIABLE` is a
parameter in the system, it will control whether the message associated with each
assertion is logged when it fails.
"""
function get_assertions_expr(sys::AbstractSystem)
    asserts = assertions(sys)
    term = 0
    if is_parameter(sys, ASSERTION_LOG_VARIABLE)
        for (k, v) in asserts
            term += _debug_assertion(k, "Assertion $k failed:\n$v", ASSERTION_LOG_VARIABLE)
        end
    else
        for (k, v) in asserts
            term += _nan_condition(k)
        end
    end
    return term
end

function SciMLBase.diagnose_symbolic_instability(sys::AbstractSystem, u, uprev)
    diagnosis = String[]

    #check for assertion failures
    unks = unknowns(sys)
    curr_substitution_map = Dict{SymbolicT, SymbolicT}(zip(unks, u))
    prev_substitution_map = Dict{SymbolicT, SymbolicT}(zip(unknowns(sys), uprev))

    for (cond, msg) in assertions(sys)
        subclauses = String[]
        find_failing_subterms(cond, prev_substitution_map, curr_substitution_map, subclauses)
        if !isempty(subclauses)
            push!(diagnosis, "\n\nAssertion violated: $cond - \"$msg\"")
            append!(diagnosis, subclauses)
        end
    end

    #find singularity causes in equations
    # The equations and the observed equations are walked separately, not as
    # `full_equations`: substituting the observed equations into the equations is
    # expensive on large systems, and printing the result expands every shared
    # subexpression (exponentially in the depth of the observed chain).
    singularities = String[]
    values = copy(prev_substitution_map)
    add_array_values!(values, unks)
    ctx = SingularityAnalysis(sys, values, singularities)
    obs = observed(sys)
    # value the observed variables that depend only on the unknowns (`obs` is sorted so
    # that each one only depends on those before it)
    for eq in obs
        val = value_at_state(eq.rhs, ctx, eq)
        val === nothing && continue
        values[eq.lhs] = val
        push!(ctx.evaluable_atoms, eq.lhs)
    end
    for eqs in (equations(sys), obs), eq in eqs
        find_singular_subterms(eq, eq.rhs, ctx)
        ctx.budget[] <= 0 && break
    end
    if ctx.unlisted[] > 0
        push!(singularities, "(and $(ctx.unlisted[]) more of these)")
    end
    if ctx.budget[] <= 0
        push!(
            singularities,
            "(analysis stopped after $DIAGNOSIS_MAX_TERMS subexpressions; the remaining equations were not checked)"
        )
    end
    if !isempty(singularities)
        push!(diagnosis, "\nSymbolic Analysis of MTK System:")
        append!(diagnosis, singularities)
    end

    return isempty(diagnosis) ? "" : join(diagnosis, "\n")
end

# The analysis runs whenever an integration fails, so its cost must stay bounded on large
# systems: it visits each distinct subexpression once, and stops after this many.
const DIAGNOSIS_MAX_TERMS = 1_000_000
# Equations and subexpressions longer than this are abbreviated in the messages.
const DIAGNOSIS_MAX_EXPRESSION_LENGTH = 200
# At most this many findings are listed; the rest are only counted.
const DIAGNOSIS_MAX_FINDINGS = 50
# A finding in an array term lists at most this many of its elements.
const DIAGNOSIS_MAX_ELEMENTS = 5

struct SingularityAnalysis{S}
    # substitutes the state the integrator failed from
    subber::S
    diagnosis::Vector{String}
    visited::IdDict{SymbolicT, Nothing}
    # what a subexpression may depend on to evaluate to a number: the analysis has no
    # parameter values, so it can only evaluate subexpressions of the unknowns and time
    evaluable_atoms::Set{SymbolicT}
    evaluable::IdDict{SymbolicT, Bool}
    budget::Base.RefValue{Int}
    # findings not listed because there were more than `DIAGNOSIS_MAX_FINDINGS`
    unlisted::Base.RefValue{Int}
    # equations with an evaluation that threw, which is reported once per equation
    failed::IdDict{Equation, Nothing}
end

function SingularityAnalysis(sys::AbstractSystem, substitution_map, diagnosis)
    subber = SymbolicUtils.IRSubstituter{true}(get_irstructure(sys), substitution_map)
    atoms = Set{SymbolicT}(keys(substitution_map))
    is_time_dependent(sys) && push!(atoms, get_iv(sys)::SymbolicT)
    return SingularityAnalysis(
        subber, diagnosis, IdDict{SymbolicT, Nothing}(), atoms, IdDict{SymbolicT, Bool}(),
        Ref(DIAGNOSIS_MAX_TERMS), Ref(0), IdDict{Equation, Nothing}()
    )
end

"""
    $(TYPEDSIGNATURES)

Whether `expr` depends only on the unknowns and time, so that substituting the state gives a
number. Substituting into anything else (a subexpression with a parameter) only builds a new,
often much larger, symbolic expression, which on large systems runs out of memory.
"""
function is_evaluable(expr, ctx::SingularityAnalysis)
    expr = unwrap(expr)
    expr isa SymbolicT || return true
    expr in ctx.evaluable_atoms && return true
    SymbolicUtils.isconst(expr) && return true
    # a variable that is not an unknown, e.g. a parameter
    SymbolicUtils.iscall(expr) || return false
    op = SymbolicUtils.operation(expr)
    # `p(t)` (a discrete), `Initial(x)`, `Pre(p)`, ...
    (op isa SymbolicT || op isa SU.Operator) && return false
    return get!(ctx.evaluable, expr) do
        all(arg -> is_evaluable(arg, ctx), SymbolicUtils.arguments(expr))
    end
end

"""
    $(TYPEDSIGNATURES)

Add to `values` the value of each array whose elements are all unknowns, so that the
unscalarized equations, which use the whole array, can be evaluated.
"""
function add_array_values!(values::AbstractDict, unks)
    seen = Set{SymbolicT}()
    for u in unks
        SU.iscall(u) && SU.operation(u) === getindex || continue
        arr = SU.arguments(u)[1]
        arr in seen && continue
        push!(seen, arr)
        sh = SU.shape(arr)
        sh isa SU.ShapeVecT && all(r -> first(r) == 1, sh) || continue
        idxs = CartesianIndices(Tuple(length.(sh)))
        all(I -> haskey(values, arr[Tuple(I)...]), idxs) || continue
        values[arr] = map(I -> SU.unwrap_const(values[arr[Tuple(I)...]]), idxs)
    end
    return values
end

# Printing it shows the error, as `showerror` does.
struct ShownError
    err::Any
end
Base.print(io::IO, e::ShownError) = showerror(io, e.err)

"""
    $(TYPEDSIGNATURES)

The value of `expr` at the failing state: a number, an array of numbers, or `nothing` if it
does not evaluate to either.

Evaluating it can throw, e.g. `sqrt` of a negative number, also in an observed equation that
the dynamics do not use. The analysis must not replace the solver's result with an exception,
so the error is reported as a finding in `eq` and `nothing` is returned. Only the evaluation
of the model's own expressions is guarded: interrupts and running out of memory or stack
are rethrown, and errors in the analysis itself are not caught.
"""
function value_at_state(expr, ctx::SingularityAnalysis, eq::Equation)
    is_evaluable(expr, ctx) || return nothing
    val = try
        Symbolics.value(ctx.subber(expr))
    catch err
        err isa Union{InterruptException, OutOfMemoryError, StackOverflowError} && rethrow()
        evaluation_failed!(eq, expr, err, ctx)
        return nothing
    end
    val isa Number && return val
    val isa AbstractArray && all(x -> x isa Number, val) && return val
    return nothing
end

# The value of `expr` at the failing state if it is a real number, otherwise `nothing`.
function real_value_at_state(expr, ctx::SingularityAnalysis, eq::Equation)
    val = value_at_state(expr, ctx, eq)
    return val isa Real ? val : nothing
end

function evaluation_failed!(eq::Equation, expr, err, ctx::SingularityAnalysis)
    # the subexpressions of a failing expression often fail too: report only the first
    haskey(ctx.failed, eq) && return nothing
    ctx.failed[eq] = nothing
    add_finding!(ctx) do
        error_text = replace(abbreviated(ShownError(err)), '\n' => ' ')
        "in equation $(abbreviated(eq)): evaluating $(abbreviated(expr)) at the failing state threw $error_text"
    end
    return nothing
end

# An `IO` that accepts at most `limit` bytes, so that printing stops early instead of
# expanding a large expression in full.
struct BoundedIO <: IO
    buf::IOBuffer
    limit::Int
end
struct OutputLimitReached <: Exception end
function Base.unsafe_write(io::BoundedIO, p::Ptr{UInt8}, n::UInt)
    room = io.limit - position(io.buf)
    unsafe_write(io.buf, p, min(n, UInt(max(room, 0))))
    n > room && throw(OutputLimitReached())
    return n
end
Base.write(io::BoundedIO, b::UInt8) = unsafe_write(io, Ref(b), UInt(1))

function abbreviated(x)
    io = BoundedIO(IOBuffer(), DIAGNOSIS_MAX_EXPRESSION_LENGTH)
    try
        print(io, x)
    catch err
        err isa OutputLimitReached || rethrow()
        text = String(take!(io.buf))
        # the limit can fall inside a multi-byte character: drop its incomplete start
        final = thisind(text, ncodeunits(text))
        if final > 0 && !isvalid(text[final])
            text = text[1:prevind(text, final)]
        end
        return text * " …"
    end
    return String(take!(io.buf))
end

# `message()` builds the finding's text, only if it is listed.
function add_finding!(message, ctx::SingularityAnalysis)
    if length(ctx.diagnosis) < DIAGNOSIS_MAX_FINDINGS
        push!(ctx.diagnosis, message())
    else
        ctx.unlisted[] += 1
    end
    return nothing
end

function find_singular_subterms(eq, expr, ctx::SingularityAnalysis)
    expr = unwrap(expr)
    !SymbolicUtils.iscall(expr) && return ctx.diagnosis
    haskey(ctx.visited, expr) && return ctx.diagnosis
    ctx.budget[] <= 0 && return ctx.diagnosis
    ctx.budget[] -= 1
    ctx.visited[expr] = nothing
    op = SymbolicUtils.operation(expr)
    args = SymbolicUtils.arguments(expr)
    diagnosis = ctx.diagnosis

    if op === (/) #division, singular if we divide by small thing
        d = real_value_at_state(args[2], ctx, eq)
        if d !== nothing && abs(d) < 1.0e-10
            add_finding!(ctx) do
                "in equation $(abbreviated(eq)): division by very small value $(abbreviated(args[2])) ≈ $(@sprintf("%.4g", d)) leads to singularity."
            end
        end
    elseif op === log #singular if we log small thing
        x = real_value_at_state(args[1], ctx, eq)
        if x !== nothing && x <= 1.0e-10
            add_finding!(ctx) do
                "in equation $(abbreviated(eq)): log of $(abbreviated(args[1])) = $(@sprintf("%.4g", x)) near/at singularity (derivative blows up)."
            end
        end
    elseif op === sqrt
        x = real_value_at_state(args[1], ctx, eq)
        if x !== nothing && x < 1.0e-10
            add_finding!(ctx) do
                "in equation $(abbreviated(eq)): sqrt of $(abbreviated(args[1])) = $(@sprintf("%.4g", x)) near/at singularity (derivative blows up)."
            end
        end
    elseif op === (^)
        e = real_value_at_state(args[2], ctx, eq)
        b = e === nothing ? nothing : real_value_at_state(args[1], ctx, eq)
        if e !== nothing && b !== nothing #two cases
            if e < 0 && abs(b) < 1.0e-10
                add_finding!(ctx) do
                    "in equation $(abbreviated(eq)): ($(abbreviated(args[1]))) raised to power $e with base ≈ $(@sprintf("%.4g", b)) going to 0; result diverges."
                end
            elseif e > 0 && abs(b) > 1
                add_finding!(ctx) do
                    "in equation $(abbreviated(eq)): ($(abbreviated(args[1])) ≈ $(@sprintf("%.4g", b))) raised to power $e - base magnitude is large and being amplified."
                end
            end
        end
    elseif op === broadcast && SU.isconst(args[1])
        # an unscalarized array equation, e.g. `D(x) ~ 1 ./ x`
        find_singular_elements(eq, SU.unwrap_const(args[1]), @view(args[2:end]), ctx)
    end

    for arg in args
        find_singular_subterms(eq, arg, ctx)
    end
    return diagnosis
end

"""
    $(TYPEDSIGNATURES)

Check the elements of the broadcast of `f` over `fargs`, with the same thresholds as the
scalar operations in `find_singular_subterms`. The term is one finding, which lists the first
`DIAGNOSIS_MAX_ELEMENTS` failing elements.
"""
function find_singular_elements(eq, f, fargs, ctx::SingularityAnalysis)
    vals = Any[]
    for arg in fargs
        val = value_at_state(arg, ctx, eq)
        val === nothing && return nothing
        push!(vals, val)
    end
    if f === (/) && length(vals) == 2
        elements = failing_elements((n, d) -> abs(d) < 1.0e-10, (n, d) -> format_value(d), vals, ctx)
        elements === nothing || add_finding!(ctx) do
            "in equation $(abbreviated(eq)): division by very small values of $(abbreviated(fargs[2])) ($elements) leads to singularity."
        end
    elseif (f === log || f === sqrt) && length(vals) == 1
        threshold_reached = f === log ? (x -> x <= 1.0e-10) : (x -> x < 1.0e-10)
        elements = failing_elements(threshold_reached, format_value, vals, ctx)
        elements === nothing || add_finding!(ctx) do
            "in equation $(abbreviated(eq)): $f of $(abbreviated(fargs[1])) ($elements) near/at singularity (derivative blows up)."
        end
    elseif f === (^) && length(vals) == 2
        show_element = (b, e) -> "base ≈ $(format_value(b)), power $e"
        elements = failing_elements((b, e) -> e < 0 && abs(b) < 1.0e-10, show_element, vals, ctx)
        elements === nothing || add_finding!(ctx) do
            "in equation $(abbreviated(eq)): ($(abbreviated(fargs[1]))) raised to negative powers with base going to 0 ($elements); result diverges."
        end
        elements = failing_elements((b, e) -> e > 0 && abs(b) > 1, show_element, vals, ctx)
        elements === nothing || add_finding!(ctx) do
            "in equation $(abbreviated(eq)): ($(abbreviated(fargs[1]))) raised to positive powers with large base magnitude ($elements) - being amplified."
        end
    end
    return nothing
end

format_value(x) = @sprintf("%.4g", x)

"""
    $(TYPEDSIGNATURES)

The elements of the broadcast over `vals` (arrays, or scalars, which broadcast) for which
`singular` holds, as text such as `"[1]: 0; [3]: 1e-12"` with each element shown by `describe`,
or `nothing` if there are none. Each element counts against the analysis budget.
"""
function failing_elements(singular, describe, vals, ctx::SingularityAnalysis)
    singular_real = (xs...) -> all(x -> x isa Real, xs) && singular(xs...)
    if !any(v -> v isa AbstractArray, vals)
        ctx.budget[] -= 1
        return singular_real(vals...) ? describe(vals...) : nothing
    end
    # the arguments' elements, broadcast to a common shape
    elements = try
        broadcast(tuple, vals...)
    catch err
        err isa DimensionMismatch || rethrow()
        return nothing
    end
    ctx.budget[] -= length(elements)
    failing = [I for I in CartesianIndices(elements) if singular_real(elements[I]...)]
    isempty(failing) && return nothing
    shown = first(failing, DIAGNOSIS_MAX_ELEMENTS)
    text = join(("[$(join(Tuple(I), ", "))]: $(describe(elements[I]...))" for I in shown), "; ")
    nmore = length(failing) - length(shown)
    return nmore > 0 ? "$text; and $nmore more" : text
end

function find_failing_subterms(cond, prev_map, curr_map, diagnosis)
    c = Symbolics.unwrap(cond)
    !SymbolicUtils.iscall(c) && return diagnosis
    op = SymbolicUtils.operation(c)
    args = SymbolicUtils.arguments(c)

    if (op === (<) || op === (>) || op === (<=) || op === (>=)) && length(args) == 2
        #compare using previous non-nan values to find violating subclauses, then output current values
        lhs = Symbolics.value(Symbolics.substitute(args[1], prev_map))
        rhs = Symbolics.value(Symbolics.substitute(args[2], prev_map))
        if lhs isa Number && rhs isa Number
            # small margin -> violated
            margin = (op === (<) || op === (<=)) ? rhs - lhs : lhs - rhs
            if margin <= 1.0e-6
                push!(diagnosis, "   subclause `$c` violated: $(clause_values(c, curr_map))")
            end
        end
    elseif op === (!=) && length(args) == 2
        lhs = Symbolics.value(Symbolics.substitute(args[1], prev_map))
        rhs = Symbolics.value(Symbolics.substitute(args[2], prev_map))
        if lhs isa Number && rhs isa Number && abs(lhs - rhs) <= 1.0e-6
            push!(diagnosis, "   subclause `$c` violated: $(clause_values(c, curr_map))")
        end
    elseif op === (==) && length(args) == 2
        lhs = Symbolics.value(Symbolics.substitute(args[1], prev_map))
        rhs = Symbolics.value(Symbolics.substitute(args[2], prev_map))
        if lhs isa Number && rhs isa Number && abs(lhs - rhs) > 1.0e-6
            push!(diagnosis, "   subclause `$c` violated: $(clause_values(c, curr_map))")
        end
    else #recurse
        for arg in args
            find_failing_subterms(arg, prev_map, curr_map, diagnosis)
        end
    end
    return diagnosis
end

function clause_values(c, curr_map)
    parts = String[]
    for v in Symbolics.get_variables(c)
        val = Symbolics.value(Symbolics.substitute(v, curr_map))
        push!(parts, val isa Number ? "$v = $(@sprintf("%.4g", val))" : "$v = $val")
    end
    return join(parts, ", ")
end
