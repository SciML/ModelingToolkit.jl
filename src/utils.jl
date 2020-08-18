function make_operation(@nospecialize(op), args)
    if op === (*)
        args = filter(!_isone, args)
        if isempty(args)
            return Constant(1)
        end
    elseif op === (+)
        args = filter(!_iszero, args)
        if isempty(args)
            return Constant(0)
        end
    end
    if !isempty(args) && all(x-> x isa Constant || !(x isa Expression), args)
        x = op(map(x->x isa Constant ? x.value : x, args)...)
        return x isa Expression ? x : Constant(x)
    else
        return op(args...)
    end
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
    expr -> to_mtk(SymbolicUtils.substitute(expr, dict))
end

@deprecate substitute_expr!(expr,s) substitute(expr,s)

# Really bad solve for vars
function solve_for(eqs, vars)
    @assert length(eqs) >= length(vars)
    @assert all(iszero(eq.lhs) for eq in eqs)
    neweqs = []
    for (i, eq) in enumerate(eqs)
        rhs = eq.rhs
        if rhs.op == (-)
            if any(isequal(rhs.args[1]), vars) && any(isequal(rhs.args[2]), vars)
                push!(neweqs, rhs.args[1] ~ rhs.args[2]) # pick one?
                @warn("todo")
            elseif any(isequal(rhs.args[1]), vars)
                push!(neweqs, rhs.args[1] ~ rhs.args[2])
            elseif any(isequal(rhs.args[2]), vars)
                push!(neweqs, rhs.args[2] ~ rhs.args[1])
            else
                @warn("may require unimplemented solve")
                #error("todo 2")
                push!(neweqs, eq)
            end
        elseif rhs.op == (+)
            eqs[i] = 0 ~ rhs.args[1] - (-rhs.args[2])
        else
            error("todo")
        end
    end
    if length(neweqs) >= length(vars)
        return neweqs
    else
        # substitute
        eqs′ = Equation.(0, substitute.(rhss(eqs), (Dict(lhss(neweqs) .=> rhss(neweqs),))))
        solve_for(eqs′, vars)
    end
end
