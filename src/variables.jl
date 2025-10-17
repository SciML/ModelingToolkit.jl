struct VariableUnit end
Symbolics.option_to_metadata_type(::Val{:unit}) = VariableUnit
function dump_variable_metadata(var)
    name = Symbolics.getmetadata(
        uvar, VariableSource, (:unknown, :unknown))
    if type <: AbstractArray
    end
    meta = (
    )
end
function normalize_to_differential(@nospecialize(op))
    if op isa Shift && op.t isa SymbolicT
    end
end
function default_toterm(x::SymbolicT)
    Moshi.Match.@match x begin
        BSImpl.Term(; f, args, shape, type, metadata) && if f isa Operator end => begin
            if f isa Shift && f.steps < 0
            end
        end
    end
end
function getbounds(x::SymbolicT)
    if operation(p) === getindex
        if symbolic_type(x) == ArraySymbolic() && symbolic_has_known_size(x)
        end
        bounds = map(bounds) do b
            if b isa AbstractArray
                if symbolic_has_known_size(p) && size(p) != size(b)
                end
            end
        end
    end
end
function istunable(x, default = true)
end
function tunable_parameters(
        sys, p = parameters(sys; initial_parameters = true); default = true)
end
macro brownians(xs...)
    all(
        xs) ||
    Symbolics.parse_vars(:brownian,
        tobrownian)
end
function getguess(x)
    if hasdefault(x) && !((def = getdefault(x)) isa Equation)
    end
end
struct EvalAt <: Symbolics.Operator
end
function (A::EvalAt)(x::SymbolicT)
    if symbolic_type(x) == NotSymbolic() || !iscall(x)
        if x isa CallAndWrap
        end
    end
end
