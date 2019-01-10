using MacroTools


function Base.convert(::Type{Expression}, ex::Expr)
    ex.head === :call || throw(ArgumentError("internal representation does not support non-call Expr"))
    f = ex.args[1]
    operands = ex.args[2:end]
    return convert(Expression, f, convert.(Expression, operands))
end
Base.convert(::Type{Expression}, sym::Symbol, args) = Operation(eval(sym), args)
Base.convert(::Type{Expression}, x::Expression) = x
Base.convert(::Type{Expression}, x::Number) = Constant(x)


function expr_arr_to_block(exprs)
  block = :(begin end)
  foreach(expr -> push!(block.args, expr), exprs)
  block
end

# used in parsing
isblock(x) = length(x) == 1 && x[1] isa Expr && x[1].head == :block
function flatten_expr!(x)
    isb = isblock(x)
    if isb
        x = MacroTools.striplines(x[1])
        filter!(z->z isa Symbol || z.head != :line, x.args)
        x = (x.args...,)
    end
    x
end

toexpr(ex) = MacroTools.postwalk(x -> isa(x, Expression) ? convert(Expr, x) : x, ex)

is_constant(x::Variable) = x.subtype === :Constant
is_constant(::Any) = false

is_operation(::Operation) = true
is_operation(::Any) = false

has_dependent(t::Variable) = Base.Fix2(has_dependent, t)
has_dependent(x::Variable, t::Variable) =
    t ∈ x.dependents || any(has_dependent(t), x.dependents)
