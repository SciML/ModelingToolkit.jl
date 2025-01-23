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
        join("\n  " .* string.(lf.args .=> args)) # one line for each "var => val" for readability
    )
end
Base.showerror(io::IO, err::LoggedFunctionException) = print(io, err.msg)
Base.nameof(lf::LoggedFun) = nameof(lf.f)
SymbolicUtils.promote_symtype(::LoggedFun, Ts...) = Real
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
