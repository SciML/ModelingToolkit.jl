"""
    pretty_expr

Convert equations to prettified expressions for latexification.
"""
function pretty_expr end

# single-argument
pretty_expr(x; kw...) = x
pretty_expr(s::Sym; kw...) = nameof(s)
pretty_expr(f::Function; kw...) = nameof(f)
pretty_expr(x::Num; kw...) = pretty_expr(value(x); kw...)
pretty_expr(ex::Expr; kw...) = Expr(ex.head, pretty_expr.(ex.args; kw...)...)
pretty_expr(eq::Equation; kw...) = Expr(:(=), pretty_expr(eq.lhs; kw...), pretty_expr(eq.rhs; kw...))
pretty_expr(f::Term; kw...) = pretty_expr(f.f, f.args; kw...)

# double-argument
pretty_expr(f, args; kw...) = Expr(:call, pretty_expr(f; kw...), pretty_expr.(args; kw...)...)
pretty_expr(f::Function, args; kw...) = Expr(:call, nameof(f), pretty_expr.(args; kw...)...)

function pretty_expr(f::Sym, args; show_iv=false, var_syms=false, kw...)
    isempty(args) && return f.name
    if show_iv # show independent variable?
        ## squish the variables and independent variable into a symbol to avoid being interpreted as a function call in latexify?
        var_syms ? Symbol(pretty_expr(f), "(", join(pretty_expr.(args), ","), ")") : Expr(:call, pretty_expr(f), pretty_expr.(args)...)
    else
        f.name
    end
end

function pretty_expr(f::Differential, args; kw...)
    ## higher-order derivative?
    if hasfield(typeof(args[1]), :f) && args[1].f isa Differential
        arg, diffs = incept_derivative(f, args, Symbol[])
        diffvars = unique(diffs)
        exponents = [count(==(x), diffs) for x in diffvars]
        return Expr(:call, 
            :/, 
            Expr(:latexifymerge, "d^{$(length(diffs))}", pretty_expr(arg; kw...)), 
            join(Symbol.(:d, string.(diffvars) .* "^{", string.(exponents), "}"))
            )
    end
    Expr(:call, :/, Expr(:latexifymerge, :d, pretty_expr(args[1]; kw...)), Symbol(:d, f.x))
end

"""
    incept_derivative

recursively figure out what's being derivated and with respect to what.
"""
function incept_derivative end
incept_derivative(f::Differential, args, diffs) = incept_derivative(args[1], vcat(diffs, f.x))
incept_derivative(x::Term, diffs) = incept_derivative(x.f, x.args, diffs)
incept_derivative(f, args, diffs) = (Term(f, args), diffs)

@latexrecipe function f(eqs::Vector{ModelingToolkit.Equation}; show_iv=true, var_syms=true)
    # Set default option values.
    env --> :align
    cdot --> false

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
