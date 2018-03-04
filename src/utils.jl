export @dvars, @idvars, @paras

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
