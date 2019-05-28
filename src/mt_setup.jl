macro mt_setup(expr)
    MacroTools.postwalk(expr) do x
        x == :ODEFunction ? :generate_function : x
    end
    expr
end
