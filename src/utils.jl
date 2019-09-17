using MacroTools


function Base.convert(::Type{Expression}, ex::Expr)
    ex.head === :if && (ex = Expr(:call, ifelse, ex.args...))
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

function build_function(rhss, vs, ps = (), args = (), conv = simplified_expr, expression = Val{true}; constructor=nothing)
    _vs = map(x-> x isa Operation ? x.op : x, vs)
    _ps = map(x-> x isa Operation ? x.op : x, ps)
    var_pairs   = [(u.name, :(u[$i])) for (i, u) ∈ enumerate(_vs)]
    param_pairs = [(p.name, :(p[$i])) for (i, p) ∈ enumerate(_ps)]
    (ls, rs) = zip(var_pairs..., param_pairs...)

    var_eqs = Expr(:(=), build_expr(:tuple, ls), build_expr(:tuple, rs))

    fname = gensym(:ModelingToolkitFunction)

    X = gensym(:MTIIPVar)
    ip_sys_exprs = [:($X[$i] = $(conv(rhs))) for (i, rhs) ∈ enumerate(rhss)]
    ip_let_expr = Expr(:let, var_eqs, build_expr(:block, ip_sys_exprs))

    sys_expr = build_expr(:tuple, [conv(rhs) for rhs ∈ rhss])
    let_expr = Expr(:let, var_eqs, sys_expr)

    fargs = ps == () ? :(u,$(args...)) : :(u,p,$(args...))

    oop_ex = :(
        ($(fargs.args...),) -> begin
            X = $let_expr
            T = promote_type(map(typeof,X)...)
            convert.(T,X)
            construct = $(constructor === nothing ? :(u isa ModelingToolkit.StaticArrays.StaticArray ? ModelingToolkit.StaticArrays.similar_type(typeof(u), eltype(X)) : x->(du=similar(u, T, $(size(rhss)...)); vec(du) .= x; du)) : constructor)
            construct(X)
        end
    )

    iip_ex = :(
        ($X,$(fargs.args...)) -> begin
            $ip_let_expr
            nothing
        end
    )

    if expression == Val{true}
        return oop_ex, iip_ex
    else
        return GeneralizedGenerated.mk_function(@__MODULE__,oop_ex), GeneralizedGenerated.mk_function(@__MODULE__,iip_ex)
    end
end

is_constant(::Constant) = true
is_constant(::Any) = false

is_operation(::Operation) = true
is_operation(::Any) = false

is_derivative(O::Operation) = isa(O.op, Differential)
is_derivative(::Any) = false

Base.occursin(t::Expression) = Base.Fix1(occursin, t)
Base.occursin(t::Expression, x::Operation ) = isequal(x, t) || any(occursin(t), x.args)
Base.occursin(t::Expression, x::Expression) = isequal(x, t)

clean(x::Variable) = x
clean(O::Operation) = isa(O.op, Variable) ? O.op : throw(ArgumentError("invalid variable: $(O.op)"))


vars(exprs) = foldl(vars!, exprs; init = Set{Variable}())
function vars!(vars, O)
    isa(O, Operation) || return vars
    for arg ∈ O.args
        if isa(arg, Operation)
            isa(arg.op, Variable) && push!(vars, arg.op)
            vars!(vars, arg)
        end
    end

    return vars
end
