const LOGGED_FUN = Set([log, sqrt, (^), /, inv])
is_legal(::typeof(/), a, b) = is_legal(inv, b)
is_legal(::typeof(inv), a) = !iszero(a)
is_legal(::Union{typeof(log), typeof(sqrt)}, a) = a isa Complex || a >= zero(a)
is_legal(::typeof(^), a, b) = a isa Complex || b isa Complex || isinteger(b) || a >= zero(a)

struct LoggedFun{F}
    f::F
    args::Any
end
Base.nameof(lf::LoggedFun) = nameof(lf.f)
SymbolicUtils.promote_symtype(::LoggedFun, Ts...) = Real
function (lf::LoggedFun)(args...)
    f = lf.f
    symbolic_args = lf.args
    if is_legal(f, args...)
        f(args...)
    else
        args_str = join(string.(symbolic_args .=> args), ", ", ", and ")
        throw(DomainError(args, "$(lf.f) errors with input(s): $args_str"))
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
