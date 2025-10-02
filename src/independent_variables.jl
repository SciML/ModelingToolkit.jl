"""
    @independent_variables t₁ t₂ ...

Define one or more independent variables. For example:

    @independent_variables t
    @variables x(t)
"""
macro independent_variables(ts...)
    Symbolics.parse_vars(:independent_variables,
        Real,
        ts,
        toiv)
end

toiv(s::SymbolicT) = GlobalScope(setmetadata(s, MTKVariableTypeCtx, PARAMETER))
toiv(s::Symbolics.Arr) = wrap(toiv(value(s)))
toiv(s::Num) = Num(toiv(value(s)))
