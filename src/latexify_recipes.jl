@latexrecipe function f(eq::ModelingToolkit.Operation)
    # Set default option values.
    cdot --> false
    eq = postwalk(x -> x isa ModelingToolkit.Constant ? x.value : x, eq)
    eq = postwalk(x -> x isa Expr && length(x.args) == 1 ? x.args[1] : x, eq)
    eq = postwalk(x -> x isa Expr && x.args[1] == :derivative && length(x.args[2].args) == 2 ? :($(Symbol(:d, x.args[2]))/($(Symbol(:d, x.args[2].args[2])))) : x, eq)
    eq = postwalk(x -> x isa Expr && x.args[1] == :derivative ? "\\frac{d\\left($(Latexify.latexraw(x.args[2]))\\right)}{d$(Latexify.latexraw(x.args[3]))}" : x, eq)
    return lhs, rhs
end

@latexrecipe function f(eq::ModelingToolkit.Constant)
    return eq.value
end

@latexrecipe function f(eq::ModelingToolkit.Equation)
    return eq.lhs, eq.rhs
end

@latexrecipe function f(sys::ModelingToolkit.AbstractSystem)
    return latexify(equations(sys))
end
