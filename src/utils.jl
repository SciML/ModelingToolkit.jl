export @dvars, @idvars, @params, @diffs, @consts

function expr_arr_to_block(exprs)
  block = :(begin end)
  foreach(expr -> push!(block.args, expr), exprs)
  block
end

# Build variables more easily
for funs in ((:dvars, :DependentVariable), (:idvars, :IndependentVariable),
             (:params, :Parameter))
    @eval begin
        macro ($(funs[1]))(x...)
            ex = Expr(:block)
            for var in x
                @assert var isa Symbol "@$($funs[1]) expects a tuple of symbols!\nE.g. `@$($funs[1]) x y z`"
                expr = :( $(esc(var)) = $($funs[2])( Symbol($(String(var))) ) )
                push!(ex.args, expr)
            end
            push!(ex.args, Expr(:tuple, esc.(x)...))
            ex
        end
    end
end

function _const_assign(x)
    ex = Expr(:block)
    lhss = Symbol[]
    for eq in x
        @assert eq isa Expr && eq.head == :(=) "@cons expects a tuple of assignments!\nE.g. `@cons D=t W=g`"
        lhs = eq.args[1]
        push!(lhss, lhs)
        rhs = eq.args[2]
        expr = :($lhs = Constant($rhs))
        push!(ex.args,  expr)
    end
    push!(ex.args, Expr(:tuple, lhss...))
    ex
end

macro consts(x...)
    esc(_const_assign(x))
end

function count_order(x)
    @assert !(x isa Symbol) "The variable $x must have a order of differentiation that is greater or equal to 1!"
    n = 1
    while !(x.args[1] isa Symbol)
        n = n+1
        x = x.args[1]
    end
    n, x.args[1]
end

function _differetial_macro(x)
    ex = Expr(:block)
    lhss = Symbol[]
    for di in x
        @assert di isa Expr && di.args[1] == :~ "@diffs expects a form that looks like `@diffs D''~t E'~t`"
        lhs = di.args[2]
        rhs = di.args[3]
        order, lhs = count_order(lhs)
        push!(lhss, lhs)
        expr = :($lhs = Differential($rhs, $order))
        push!(ex.args,  expr)
    end
    push!(ex.args, Expr(:tuple, lhss...))
    ex
end

macro diffs(x...)
    esc(_differetial_macro(x))
end
