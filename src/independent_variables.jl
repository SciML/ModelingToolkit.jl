"""
    @independent_variables t₁ t₂ ...

Define one or more independent variables. For example:

    @independent_variables t
    @variables x(t)
"""
macro independent_variables(ts...)
    :(@parameters $(ts...)) |> esc # TODO: treat independent variables separately from variables and parameters
end

toiv(s::Symbolic) = setmetadata(s, MTKVariableTypeCtx, PARAMETER)
toiv(s::Num) = Num(toiv(value(s)))
