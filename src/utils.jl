function make_operation(@nospecialize(op), args)
    if op === (*)
        args = filter(!_isone, args)
        if isempty(args)
            return 1
        end
    elseif op === (+)
        args = filter(!_iszero, args)
        if isempty(args)
            return 0
        end
    end
    return op(args...)
end

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

"""
get_variables(O)

Returns the variables in the expression
"""
get_variables(e::Num, varlist=nothing) = get_variables(value(e), varlist)
get_variables!(vars, e, varlist=nothing) = vars
get_variables!(vars, e::Sym, varlist=nothing) = push!(vars, e)

function get_variables!(vars, e::Term, varlist=nothing)
    foreach(x -> get_variables!(vars, x, varlist), e.args)
    return (vars isa AbstractVector) ? unique!(vars) : vars
end

function get_variables!(vars, e::Equation, varlist=nothing)
  get_variables!(vars, e.rhs, varlist)
end

modified_states!(mstates, e::Equation, statelist=nothing) = get_variables!(mstates, e.lhs, statelist)

get_variables(e, varlist=nothing) = get_variables!([], e, varlist)

# variable substitution
# Piracy but mild
"""
substitute(expr::Operation, s::Pair)
substitute(expr::Operation, s::Dict)
substitute(expr::Operation, s::Vector)

Performs the substitution `Operation => val` on the `expr` Operation.
"""
substitute(expr::Num, s::Union{Pair, Vector, Dict}; kw...) = Num(substitute(value(expr), s; kw...))
substitute(expr::Term, s::Pair; kw...) = substitute(expr, Dict(s[1] => s[2]); kw...)
substitute(expr::Term, s::Vector; kw...) = substitute(expr, Dict(s); kw...)

@deprecate substitute_expr!(expr,s) substitute(expr,s)
