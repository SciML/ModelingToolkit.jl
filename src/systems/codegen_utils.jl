""""""
function isdelay(var, iv)
    iv === nothing && return false
    if iscall(var) && ModelingToolkit.isoperator(var, Differential)
        return isdelay(arguments(var)[1], iv)
    end
    isvariable(var) || return false
    isparameter(var) && return false
    if iscall(var) && !ModelingToolkit.isoperator(var, Symbolics.Operator)
        args = arguments(var)
        length(args) == 1 || return false
        arg = args[1]
        isequal(arg, iv) && return false
        iscall(arg) || return true
        issym(operation(arg)) && !iscalledparameter(arg) && return false
        return true
    end
    return false
end
""""""
const DDE_HISTORY_FUN = SSym(:___history___; type = SU.FnType{Tuple{Any, <:Real}, Vector{Real}}, shape = SU.Unknown(1))
const BVP_SOLUTION = SSym(:__sol__; type = Symbolics.FnType{Tuple{<:Real}, Vector{Real}}, shape = SU.Unknown(1))
""""""
function delay_to_function(
        sys::AbstractSystem, eqs = full_equations(sys); param_arg = MTKPARAMETERS_ARG, histfn = DDE_HISTORY_FUN)
    delay_to_function(eqs,
        get_iv(sys),
        Dict{Any, Int}(operation(s) => i for (i, s) in enumerate(unknowns(sys))),
        parameters(sys),
        histfn; param_arg)
end
function delay_to_function(eqs::Vector, iv, sts, ps, h; param_arg = MTKPARAMETERS_ARG)
    delay_to_function.(eqs, (iv,), (sts,), (ps,), (h,); param_arg)
end
function delay_to_function(eq::Equation, iv, sts, ps, h; param_arg = MTKPARAMETERS_ARG)
    delay_to_function(eq.lhs, iv, sts, ps, h; param_arg) ~ delay_to_function(
        eq.rhs, iv, sts, ps, h; param_arg)
end
function delay_to_function(expr, iv, sts, ps, h; param_arg = MTKPARAMETERS_ARG)
    if isdelay(expr, iv)
        v = operation(expr)
        time = arguments(expr)[1]
        idx = sts[v]
        return term(getindex, h(param_arg, time), idx, type = Real)
    elseif iscall(expr)
        return maketerm(typeof(expr),
            operation(expr),
            map(x -> delay_to_function(x, iv, sts, ps, h; param_arg), arguments(expr)),
            metadata(expr))
    else
        return expr
    end
end
struct GeneratedFunctionWrapper{P, O, I} <: Function
    f_oop::O
    f_iip::I
end
