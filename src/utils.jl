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
get_variables!(vars, e::Constant, varlist=nothing) = vars
get_variables(e::Constant, varlist=nothing) = get_variables!(Operation[], e, varlist)

function get_variables!(vars, e::Operation, varlist=nothing)
    if is_singleton(e)
      (isnothing(varlist) ? true : (e.op in varlist)) && push!(vars, e)
    else
        foreach(x -> get_variables!(vars, x, varlist), e.args)
    end
    return (vars isa AbstractVector) ? unique!(vars) : vars
end
get_variables(e::Operation, varlist=nothing) = get_variables!(Operation[], e, varlist)

function get_variables!(vars, e::Equation, varlist=nothing)
  get_variables!(vars, e.rhs, varlist)
end
get_variables(e::Equation, varlist=nothing) = get_variables!(Operation[],e,varlist)

modified_states!(mstates, e::Equation, statelist=nothing) = get_variables!(mstates, e.lhs, statelist)

# variable substitution
"""
substitute(expr::Operation, s::Pair)
substitute(expr::Operation, s::Dict)
substitute(expr::Operation, s::Vector)

Performs the substitution `Operation => val` on the `expr` Operation.
"""
substitute(expr::Constant, s) = expr
substitute(expr::Operation, s::Pair) = substituter([s[1] => s[2]])(expr)
substitute(expr::Operation, s::Union{Vector, Dict}) = substituter(s)(expr)

function substituter(pairs)
    dict = Dict(to_symbolic(k) => to_symbolic(v)  for (k, v) in pairs)
    expr -> to_mtk(SymbolicUtils.simplify(SymbolicUtils.substitute(expr, dict)))
end

@deprecate substitute_expr!(expr,s) substitute(expr,s)
