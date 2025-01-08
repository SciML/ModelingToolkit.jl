const LOGGED_FUN = Set([log, sqrt, (^), /, inv])

struct LoggedFunctionException <: Exception
    msg::String
end
struct LoggedFun{F}
    f::F
    args::Any
end
Base.showerror(io::IO, err::LoggedFunctionException) = print(io, err.msg)
Base.nameof(lf::LoggedFun) = nameof(lf.f)
SymbolicUtils.promote_symtype(::LoggedFun, Ts...) = Real
function (lf::LoggedFun)(args...)
    try
        return lf.f(args...) # try to call with numerical input, as usual
    catch err
        throw(LoggedFunctionException(
            "Function $(lf.f)($(join(lf.args, ", "))) errors with input" *
            join("\n  " .* string.(lf.args .=> args)) # one line for each "var => val" for readability
        )) # Julia automatically attaches original error message
    end
end

function logged_fun(f, args...)
    # Currently we don't really support complex numbers
    term(LoggedFun(f, args), args..., type = Real)
end

debug_sub(eq::Equation) = debug_sub(eq.lhs) ~ debug_sub(eq.rhs)
function debug_sub(ex)
    iscall(ex) || return ex
    f = operation(ex)
    args = map(debug_sub, arguments(ex))
    f in LOGGED_FUN ? logged_fun(f, args...) :
    maketerm(typeof(ex), f, args, metadata(ex))
end
