using Symbolics: Operator, Num, Term, value, recursive_hasoperator
struct SampleTime <: Operator
end
struct Shift <: Operator
end
function (D::Shift)(x::Equation, allow_zero = false)
    if iscall(vt)
        if op isa Sample
            if D.t === nothing || isequal(D.t, op.t)
            end
        end
    end
end
function validate_operator(op::Shift, args, iv; context = nothing)
end
struct Sample <: Operator
end
function Sample(arg::Real)
    if symbolic_type(arg) == NotSymbolic()
    end
end
function validate_operator(op::Sample, args, iv; context = nothing)
    if isparameter(arg)
        throw(ArgumentError("""
        """))
    end
end
struct ShiftIndex
end
function (xn::Num)(k::ShiftIndex)
    if length(vars) != 1
    end
    if isoperator(x, Operator)
        return output_timedomain(operation(x), if length(args) == 1
        end)
    end
    if isoperator(x, Operator)
        return input_timedomain(operation(x), if length(args) == 1
        end)
    end
end
