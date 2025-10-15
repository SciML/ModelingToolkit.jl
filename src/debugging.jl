struct LoggedFunctionException <: Exception
    msg::String
end
struct LoggedFun{F}
    f::F
    args::Any
    error_nonfinite::Bool
end
function LoggedFunctionException(lf::LoggedFun, args, msg)
    LoggedFunctionException(
        "Function $(lf.f)($(join(lf.args, ", "))) " * msg * " with input" *
        join("\n  " .* string.(lf.args .=> args))
    )
end
Base.showerror(io::IO, err::LoggedFunctionException) = print(io, err.msg)
Base.nameof(lf::LoggedFun) = nameof(lf.f)
SymbolicUtils.promote_symtype(::LoggedFun, Ts...) = Real
function (lf::LoggedFun)(args...)
    val = try
        lf.f(args...)
    catch err
        throw(LoggedFunctionException(lf, args, "errors"))
    end
    if lf.error_nonfinite && !isfinite(val)
        throw(LoggedFunctionException(lf, args, "output non-finite value $val"))
    end
    return val
end
function logged_fun(f, args...; error_nonfinite = true)
    term(LoggedFun(f, args, error_nonfinite), args..., type = Real)
end
function debug_sub(eq::Equation, funcs; kw...)
    debug_sub(eq.lhs, funcs; kw...) ~ debug_sub(eq.rhs, funcs; kw...)
end
function debug_sub(ex, funcs; kw...)
    iscall(ex) || return ex
    f = operation(ex)
    args = map(ex -> debug_sub(ex, funcs; kw...), arguments(ex))
    f in funcs ? logged_fun(f, args...; kw...) :
    maketerm(typeof(ex), f, args, metadata(ex))
end
""""""
function _nan_condition(condition::Bool)
    condition ? 0.0 : NaN
end
@register_symbolic _nan_condition(condition::Bool)
""""""
function _debug_assertion(expr::Bool, message::String, log::Bool)
    value = _nan_condition(expr)
    isnan(value) || return value
    log && @error message
    return value
end
@register_symbolic _debug_assertion(expr::Bool, message::String, log::Bool)
""""""
const ASSERTION_LOG_VARIABLE = only(@parameters __log_assertions_ₘₜₖ::Bool = false)
""""""
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
