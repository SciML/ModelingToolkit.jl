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
function _nan_condition(condition::Bool)
    return condition ? 0.0 : NaN
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
    singularities = String[]
    visited = IdDict{SymbolicT, Nothing}()
    subber = SymbolicUtils.IRSubstituter{true}(get_irstructure(sys), prev_substitution_map)
    for eq in full_equations(sys)
        find_singular_subterms(eq, eq.rhs, subber, singularities, visited)
    end
    if !isempty(singularities)
        push!(diagnosis, "\nSymbolic Analysis of MTK System:")
        append!(diagnosis, singularities)
    end

    return isempty(diagnosis) ? "" : join(diagnosis, "\n")
end

function find_singular_subterms(eq, expr, sub_map, diagnosis, visited)
    expr = unwrap(expr)
    !SymbolicUtils.iscall(expr) && return diagnosis
    op = SymbolicUtils.operation(expr)
    args = SymbolicUtils.arguments(expr)
    haskey(visited, expr) && return diagnosis  
    visited[expr] = nothing

    if op === (/) #division, singular if we divide by small thing
        d = Symbolics.value(sub_map(args[2]))
        if d isa Number && abs(d) < 1e-10
            push!(diagnosis, "in equation $eq: division by very small value $(args[2]) ≈ $(@sprintf("%.4g", d)) leads to singularity.")
        end
    elseif op === log #singular if we log small thing
        x = Symbolics.value(sub_map(args[1]))
        if x isa Number && x <= 1e-10
            push!(diagnosis, "in equation $eq: log of $(args[1]) = $(@sprintf("%.4g", x)) near/at singularity (derivative blows up).")
        end
    elseif op === sqrt 
        x = Symbolics.value(sub_map(args[1]))
        if x isa Number && x < 1e-10
            push!(diagnosis, "in equation $eq: sqrt of $(args[1]) = $(@sprintf("%.4g", x)) near/at singularity (derivative blows up).")
        end
    elseif op === (^)
        e = Symbolics.value(sub_map(args[2]))
        b = Symbolics.value(sub_map(args[1]))
        if e isa Number && b isa Number #two cases
            if e < 0 && abs(b) < 1e-10
                push!(diagnosis, "in equation $eq: ($(args[1])) raised to power $e with base ≈ $(@sprintf("%.4g", b)) going to 0; result diverges.")
            elseif e > 0 && abs(b) > 1
                push!(diagnosis, "in equation $eq: ($(args[1]) ≈ $(@sprintf("%.4g", b))) raised to power $e - base magnitude is large and being amplified.")
            end
        end
    end

    for arg in args
        find_singular_subterms(eq, arg, sub_map, diagnosis, visited)
    end
    return diagnosis
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
            if margin <= 1e-6
                push!(diagnosis, "   subclause `$c` violated: $(clause_values(c, curr_map))")
            end
        end
    elseif op === (!=) && length(args) == 2
        lhs = Symbolics.value(Symbolics.substitute(args[1], prev_map))
        rhs = Symbolics.value(Symbolics.substitute(args[2], prev_map))
        if lhs isa Number && rhs isa Number && abs(lhs - rhs) <= 1e-6
            push!(diagnosis, "   subclause `$c` violated: $(clause_values(c, curr_map))")
        end
    elseif op === (==) && length(args) == 2
        lhs = Symbolics.value(Symbolics.substitute(args[1], prev_map))
        rhs = Symbolics.value(Symbolics.substitute(args[2], prev_map))
        if lhs isa Number && rhs isa Number && abs(lhs - rhs) > 1e-6
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
