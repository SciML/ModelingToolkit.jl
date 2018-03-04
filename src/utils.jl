export @dvars, @idvars, @paras, @diffs, @cons

function expr_arr_to_block(exprs)
  block = :(begin end)
  foreach(expr -> push!(block.args, expr), exprs)
  block
end

# Build variables more easily
for funs in ((:dvars, :DependentVariable), (:idvars, :IndependentVariable),
             (:paras, :Parameter))
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

function assign(macroname, x, typ)
    ex = Expr(:block)
    lhss = Symbol[]
    for eq in x
        @assert eq isa Expr && eq.head == :(=) "@$macroname expects a tuple of assignments!\nE.g. `@$macroname D=t W=g`"
        lhs = eq.args[1]
        push!(lhss, lhs)
        rhs = eq.args[2]
        expr = :($lhs = $typ($rhs))
        push!(ex.args,  expr)
    end
    push!(ex.args, Expr(:tuple, lhss...))
    ex
end

for funs in ((:cons, :Constant), (:diffs, :Differential))
    @eval begin
        macro ($(funs[1]))(x...)
            esc(assign(String($funs[1]), x, $funs[2]))
        end
    end
end
