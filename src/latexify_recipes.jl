@latexrecipe function f(eqs::Vector{ModelingToolkit.Equation}; iv=:t)
    # Set default option values.
    env --> :align

    # Convert both the left and right hand side to expressions of basic types
    # that latexify can deal with.

    rhs = getfield.(eqs, :rhs)
    rhs = convert.(Expr, rhs)
    rhs = [postwalk(x -> x isa ModelingToolkit.Constant ? x.value : x, eq) for eq in rhs]
    rhs = [postwalk(x -> x isa Expr && length(x.args) == 1 ? x.args[1] : x, eq) for eq in rhs]
    rhs = [postwalk(x -> x isa Expr && x.args[1] == :Differential && length(x.args[2].args) == 2 ? :($(Symbol(:d, x.args[2]))/($(Symbol(:d, x.args[2].args[2])))) : x, eq) for eq in rhs]
    rhs = [postwalk(x -> x isa Expr && x.args[1] == :Differential ? "\\frac{d\\left($(Latexify.latexraw(x.args[2]))\\right)}{d$iv}" : x, eq) for eq in rhs]

    lhs = getfield.(eqs, :lhs)
    lhs = convert.(Expr, lhs)
    lhs = [postwalk(x -> x isa ModelingToolkit.Constant ? x.value : x, eq) for eq in lhs]
    lhs = [postwalk(x -> x isa Expr && length(x.args) == 1 ? x.args[1] : x, eq) for eq in lhs]
    lhs = [postwalk(x -> x isa Expr && x.args[1] == :Differential && length(x.args[2].args) == 2 ? :($(Symbol(:d, x.args[2]))/($(Symbol(:d, x.args[2].args[2])))) : x, eq) for eq in lhs]
    lhs = [postwalk(x -> x isa Expr && x.args[1] == :Differential ? "\\frac{d\\left($(Latexify.latexraw(x.args[2]))\\right)}{d$iv}" : x, eq) for eq in lhs]

    return lhs, rhs
end

@latexrecipe function f(eqs::Vector{ModelingToolkit.DiffEq}; iv=:t)
    # Set default option values.
    env --> :align

    # Convert both the left and right hand side to expressions of basic types
    # that latexify can deal with.

    rhs = getfield.(eqs, :rhs)
    rhs = convert.(Expr, rhs)
    rhs = [postwalk(x -> x isa ModelingToolkit.Constant ? x.value : x, eq) for eq in rhs]
    rhs = [postwalk(x -> x isa Expr && length(x.args) == 1 ? x.args[1] : x, eq) for eq in rhs]
    rhs = [postwalk(x -> x isa Expr && x.args[1] == :Differential && length(x.args[2].args) == 2 ? :($(Symbol(:d, x.args[2]))/($(Symbol(:d, x.args[2].args[2])))) : x, eq) for eq in rhs]
    rhs = [postwalk(x -> x isa Expr && x.args[1] == :Differential ? "\\frac{d\\left($(Latexify.latexraw(x.args[2]))\\right)}{d$iv}" : x, eq) for eq in rhs]

    var = getfield.(getfield.(eqs, :x),:name)
    ns  = getfield.(eqs, :n)
    lhs = [ns[i] == 1 ? Latexify.LaTeXString("\\frac{d$(Latexify.latexraw(var[i]))}{d$iv}") : LaTeXString("\\frac{d^{$(ns[i])}$(Latexify.latexraw(var[i]))}{d$iv^{$(ns[i])}}") for i in 1:length(var)]
    return lhs, rhs
end

@latexrecipe function f(sys::ModelingToolkit.ODESystem; iv=:t)
    sys.eqs
end
