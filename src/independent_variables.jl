""""""
macro independent_variables(ts...)
    Symbolics.parse_vars(:independent_variables,
        Real,
        ts,
        toiv)
end
toiv(s::SymbolicT) = GlobalScope(setmetadata(s, MTKVariableTypeCtx, PARAMETER))
toiv(s::Symbolics.Arr) = wrap(toiv(value(s)))
toiv(s::Num) = Num(toiv(value(s)))
