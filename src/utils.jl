function expr_arr_to_block(exprs)
  block = :(begin end)
  foreach(expr -> push!(block.args, expr), exprs)
  block
end
