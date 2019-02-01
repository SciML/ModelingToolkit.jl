using MacroTools


function Base.convert(::Type{Expression}, ex::Expr)
    ex.head === :call || throw(ArgumentError("internal representation does not support non-call Expr"))

    op = eval(ex.args[1])  # HACK
    args = convert.(Expression, ex.args[2:end])

    return Operation(op, args)
end
Base.convert(::Type{Expression}, x::Expression) = x
Base.convert(::Type{Expression}, x::Number) = Constant(x)

function build_expr(head::Symbol, args)
    ex = Expr(head)
    append!(ex.args, args)
    ex
end
expr_arr_to_block(exprs) = build_expr(:block, exprs)

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

function partition(f, xs)
    idxs = map(f, xs)
    not_idxs = eachindex(xs) .∉ (idxs,)
    return (xs[idxs], xs[not_idxs])
end

is_constant(::Constant) = true
is_constant(::Any) = false

is_operation(::Operation) = true
is_operation(::Any) = false

has_dependent(t::Operation) = Base.Fix2(has_dependent, t)
has_dependent(x::Operation, t::Operation) = x == t || any(has_dependent(t), x.args)
has_dependent(::Expression, ::Expression) = false
