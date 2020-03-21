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
Base.convert(::Type{Expression}, x::Bool) = Constant(x)
Expression(x::Bool) = Constant(x)

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

function build_function(rhss, vs, ps = (), args = (), conv = simplified_expr, expression = Val{true};
                        checkbounds = false, constructor=nothing, linenumbers = true)
    _vs = map(x-> x isa Operation ? x.op : x, vs)
    _ps = map(x-> x isa Operation ? x.op : x, ps)
    var_pairs   = [(u.name, :(u[$i])) for (i, u) ∈ enumerate(_vs)]
    param_pairs = [(p.name, :(p[$i])) for (i, p) ∈ enumerate(_ps)]
    (ls, rs) = zip(var_pairs..., param_pairs...)

    var_eqs = Expr(:(=), build_expr(:tuple, ls), build_expr(:tuple, rs))

    fname = gensym(:ModelingToolkitFunction)

    X = gensym(:MTIIPVar)
    if rhss isa SparseMatrixCSC
        ip_sys_exprs = [:($X.nzval[$i] = $(conv(rhs))) for (i, rhs) ∈ enumerate(rhss.nzval)]
    else
        ip_sys_exprs = [:($X[$i] = $(conv(rhs))) for (i, rhs) ∈ enumerate(rhss)]
    end

    ip_let_expr = Expr(:let, var_eqs, build_expr(:block, ip_sys_exprs))

    tuple_sys_expr = build_expr(:tuple, [conv(rhs) for rhs ∈ rhss])

    if rhss isa Matrix
        arr_sys_expr = build_expr(:vcat, [build_expr(:row,[conv(rhs) for rhs ∈ rhss[i,:]]) for i in 1:size(rhss,1)])
        # : x because ??? what to do in the general case?
        _constructor = constructor === nothing ? :(u isa ModelingToolkit.StaticArrays.StaticArray ? ModelingToolkit.StaticArrays.SMatrix{$(size(rhss)...)} :  x->(out = similar(typeof(u),$(size(rhss)...)); out .= x)) : constructor
    elseif typeof(rhss) <: Array && !(typeof(rhss) <: Vector)
        vector_form = build_expr(:vect, [conv(rhs) for rhs ∈ rhss])
        arr_sys_expr = :(reshape($vector_form,$(size(rhss)...)))
        _constructor = constructor === nothing ? :(u isa ModelingToolkit.StaticArrays.StaticArray ? ModelingToolkit.StaticArrays.SArray{$(size(rhss)...)} :  x->(out = similar(typeof(u),$(size(rhss)...)); out .= x)) : constructor
    elseif rhss isa SparseMatrixCSC
        vector_form = build_expr(:vect, [conv(rhs) for rhs ∈ nonzeros(rhss)])
        arr_sys_expr = :(SparseMatrixCSC{eltype(u),Int}($(size(rhss)...), $(rhss.colptr), $(rhss.rowval), $vector_form))
        # Static and sparse? Probably not a combo that will actually be hit, but give a default anyways
        _constructor = constructor === nothing ? :(u isa ModelingToolkit.StaticArrays.StaticArray ? ModelingToolkit.StaticArrays.SMatrix{$(size(rhss)...)} : x->x) : constructor
    else # Vector
        arr_sys_expr = build_expr(:vect, [conv(rhs) for rhs ∈ rhss])
        # Handle vector constructor separately using `typeof(u)` to support things like LabelledArrays
        _constructor = constructor === nothing ? :(u isa ModelingToolkit.StaticArrays.StaticArray ? ModelingToolkit.StaticArrays.similar_type(typeof(u), eltype(X)) : x->convert(typeof(u),x)) : constructor
    end

    let_expr = Expr(:let, var_eqs, tuple_sys_expr)
    arr_let_expr = Expr(:let, var_eqs, arr_sys_expr)
    bounds_block = checkbounds ? let_expr : :(@inbounds begin $let_expr end)
    arr_bounds_block = checkbounds ? arr_let_expr : :(@inbounds begin $arr_let_expr end)
    ip_bounds_block = checkbounds ? ip_let_expr : :(@inbounds begin $ip_let_expr end)

    fargs = ps == () ? :(u,$(args...)) : :(u,p,$(args...))

    oop_ex = :(
        ($(fargs.args...),) -> begin
            # If u is a weird non-StaticArray type and we want a sparse matrix, just do the optimized sparse anyways
            if $(fargs.args[1]) isa Array || (!(typeof($(fargs.args[1])) <: StaticArray) && $(rhss isa SparseMatrixCSC))
                return $arr_bounds_block
            else
                X = $bounds_block
                construct = $_constructor
                return construct(X)
            end
        end
    )

    iip_ex = :(
        ($X,$(fargs.args...)) -> begin
            $ip_bounds_block
            nothing
        end
    )

    if !linenumbers
        oop_ex = striplines(oop_ex)
        iip_ex = striplines(iip_ex)
    end

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

# variable extraction
is_singleton(e) = false
is_singleton(e::Operation) = e.op isa Variable

get_variables(e::ModelingToolkit.Constant, vars) = nothing
function get_variables(e::Expression, vars = Operation[])
    if is_singleton(e)
        push!(vars, e)
    else
        foreach(x -> get_variables(x, vars), e.args)
    end
    return unique(vars)
end

# variable substitution
function substitute_expr!(expr::Expression, s::Pair{Operation, Operation})
    if is_singleton(expr) || expr isa ModelingToolkit.Constant
        # do nothing
    else
        expr.args .= replace(expr.args, s)
        [substitute_expr!(arg, s) for arg in expr.args if !is_singleton(arg)]
    end
    return nothing
end
