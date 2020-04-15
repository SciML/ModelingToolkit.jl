Base.convert(::Type{Expression}, ex::Expr) = Expression(ex)
function Expression(ex;mod=Main)
    ex.head === :if && (ex = Expr(:call, ifelse, ex.args...))
    ex.head === :call || throw(ArgumentError("internal representation does not support non-call Expr"))

    op = getproperty(mod,ex.args[1])
    args = convert.(Expression, ex.args[2:end])

    return Operation(op, args)
end
Base.convert(::Type{Expression}, x::Expression) = x
Base.convert(::Type{Expression}, x::Number) = Constant(x)
Base.convert(::Type{Expression}, x::Bool) = Constant(x)
Base.convert(::Type{Expression}, x::Variable) = convert(Operation,x)
Base.convert(::Type{Expression}, x::Operation) = x
Base.convert(::Type{Expression}, x::Symbol) = Operation(Variable(x),[])
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

function detime_dvs(op::Operation)
  if op.op isa Variable
    Operation(Variable{vartype(op.op)}(op.op.name),Expression[])
  else
    Operation(op.op,detime_dvs.(op.args))
  end
end
detime_dvs(op::Constant) = op

function retime_dvs(op::Operation,dvs,iv)
  if op.op isa Variable && op.op ∈ dvs
    Operation(Variable{vartype(op.op)}(op.op.name),Expression[iv])
  else
    Operation(op.op,retime_dvs.(op.args,(dvs,),iv))
  end
end
retime_dvs(op::Constant,dvs,iv) = op

is_constant(::Constant) = true
is_constant(::Any) = false

is_operation(::Operation) = true
is_operation(::Any) = false

is_derivative(O::Operation) = isa(O.op, Differential)
is_derivative(::Any) = false

Base.occursin(t::Expression) = Base.Fix1(occursin, t)
Base.occursin(t::Expression, x::Operation ) = isequal(x, t) || any(occursin(t), x.args)
Base.occursin(t::Expression, x::Expression) = isequal(x, t)

vars(exprs) = foldl(vars!, exprs; init = Set{Variable}())
function vars!(vars, O)
    isa(O, Operation) || return vars
    O.op isa Variable && push!(vars, O.op)
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

"""
get_variables(O::Operation)

Returns the variables in the Operation
"""
get_variables(e::Constant, vars = Operation[]) = vars
function get_variables(e::Operation, vars = Operation[])
    if is_singleton(e)
        push!(vars, e)
    else
        foreach(x -> get_variables(x, vars), e.args)
    end
    return unique(vars)
end

# variable substitution
"""
substitute_expr!(expr::Operation, s::Pair{Operation, Operation})

Performs the substitution `Operation => Operation` on the `expr` Operation.
"""
substitute_expr!(expr::Constant, s::Pair{Operation, Operation}) = nothing
function substitute_expr!(expr::Operation, s::Pair{Operation, Operation})
    if !is_singleton(expr)
        expr.args .= replace(expr.args, s)
        for arg in expr.args
            substitute_expr!(arg, s)
        end
    end
    return nothing
end
