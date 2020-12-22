"""
    pretty_expr

Convert equations to prettified expressions for latexification.
"""
function pretty_expr end

pretty_expr(expr; kwargs...) = expr
pretty_expr(s::Sym; kwargs...) = nameof(s)
pretty_expr(eq::Equation; kwargs...) = Expr(:(=), pretty_expr(eq.lhs; kwargs...), pretty_expr(eq.rhs; kwargs...))
pretty_expr(eq::AbstractArray; kwargs...) = pretty_expr.(eq; kwargs...)
pretty_expr(x::Integer; kwargs...) = x
pretty_expr(x::AbstractFloat; kwargs...) = x
pretty_expr(x::Num; kwargs...) = pretty_expr(value(x); kwargs...)
pretty_expr(expr::Expr; kwargs...) = Expr(expr.head, pretty_expr.(expr.args; kwargs...)...)
pretty_expr(f::Function; kwargs...) = nameof(f)
pretty_expr(f::Term; kwargs...) = pretty_expr(f.f, f.args; kwargs...)

pretty_expr(f, args; kwargs...) = Expr(:call, pretty_expr(f; kwargs...), pretty_expr.(args; kwargs...)...)
pretty_expr(f::Function, args; kwargs...) = Expr(:call, nameof(f), pretty_expr.(args; kwargs...)...)

function pretty_expr(f::Sym, args; show_iv=false, var_syms=false, kwargs...)
    isempty(args) && return f.name
    if show_iv
        var_syms ? Symbol(pretty_expr(f), "(", join(pretty_expr.(args), ","), ")") : Expr(:call, pretty_expr(f), pretty_expr.(args)...)
    else
        return f.name
    end
end

function pretty_expr(f::Differential, args; kwargs...)
    ## higher-order derivative?
    if hasfield(typeof(args[1]), :f) && args[1].f isa Differential
        arg, diffs = incept_derivative(f, args, Symbol[])
        diffvars = unique(diffs)
        exponents = [count(==(x), diffs) for x in diffvars]
        return Expr(:call, 
            :/, 
            Expr(:latexifymerge, "d^{$(length(diffs))}", pretty_expr(arg; kwargs...)), 
            join(Symbol.(:d, string.(diffvars) .* "^{", string.(exponents), "}"))
            )
    end
    Expr(:call, :/, Expr(:latexifymerge, :d, pretty_expr(args[1]; kwargs...)), Symbol(:d, f.x))
end
incept_derivative(f::Differential, args, diffs) = incept_derivative(args[1], vcat(diffs, f.x))
incept_derivative(x::Term, diffs) = incept_derivative(x.f, x.args, diffs)
incept_derivative(f, args, diffs) = (Term(f, args), diffs)

@latexrecipe function f(eqs::Vector{ModelingToolkit.Equation}; show_iv=true, var_syms=true)
    # Set default option values.
    env --> :align
    cdot --> false

    # Convert both the left- and right-hand sides to expressions of basic types
    # that latexify can deal with

    rhs = getfield.(eqs, :rhs)
    rhs = pretty_expr.(rhs; show_iv=show_iv, var_syms=var_syms)

    lhs = getfield.(eqs, :lhs)
    lhs = pretty_expr.(lhs; show_iv=show_iv, var_syms=var_syms)
    return lhs, rhs
end

@latexrecipe function f(sys::ModelingToolkit.AbstractSystem)
    ((args,), kw) = Latexify.apply_recipe(equations(sys); kwargs...)
    kwargs = merge(kw, kwargs)
    return args
end 
