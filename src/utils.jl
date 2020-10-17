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

function detime_dvs(op::Term)
  if op.op isa Sym
      Sym{Real}(nameof(op.op))
  else
    Term(op.op,detime_dvs.(op.args))
  end
end
detime_dvs(op) = op

function retime_dvs(op::Sym,dvs,iv)
    Sym{FnType{Tuple{symtype(iv)}, Real}}(nameof(op))(iv)
end

function retime_dvs(op::Term, dvs, iv)
    similarterm(op, op.op, retime_dvs.(op.args,(dvs,),(iv,)))
end
retime_dvs(op,dvs,iv) = op

is_constant(::Constant) = true
is_constant(::Any) = false

is_operation(::Operation) = true
is_operation(::Any) = false

is_derivative(O::Term) = isa(O.op, Differential)
is_derivative(::Any) = false

"""
get_variables(O)

Returns the variables in the expression
"""
get_variables(e::Num, varlist=nothing) = get_variables(value(e), varlist)
get_variables!(vars, e, varlist=nothing) = vars

is_singleton(e::Term) = e.op isa Sym
is_singleton(e::Sym) = true
is_singleton(e) = false

get_variables!(vars, e::Number, varlist=nothing) = vars

function get_variables!(vars, e::Symbolic, varlist=nothing)
    if is_singleton(e)
        if isnothing(varlist) || any(isequal(e), varlist)
            push!(vars, e)
        end
    else
        foreach(x -> get_variables!(vars, x, varlist), e.args)
    end
    return (vars isa AbstractVector) ? unique!(vars) : vars
end

function get_variables!(vars, e::Equation, varlist=nothing)
  get_variables!(vars, e.rhs, varlist)
end

get_variables(e, varlist=nothing) = get_variables!([], e, varlist)

modified_states!(mstates, e::Equation, statelist=nothing) = get_variables!(mstates, e.lhs, statelist)


# variable substitution
# Piracy but mild
"""
substitute(expr::Operation, s::Pair)
substitute(expr::Operation, s::Dict)
substitute(expr::Operation, s::Vector)

Performs the substitution `Operation => val` on the `expr` Operation.
"""
substitute(expr::Num, s::Union{Pair, Vector, Dict}; kw...) = Num(substituter(s)(value(expr); kw...))
# TODO: move this to SymbolicUtils
substitute(expr::Term, s::Pair; kw...) = substituter([s[1] => s[2]])(expr; kw...)
substitute(expr::Term, s::Vector; kw...) = substituter(s)(expr; kw...)

substituter(pair::Pair) = substituter((pair,))
function substituter(pairs)
    dict = Dict(to_symbolic(k) => to_symbolic(v)  for (k, v) in pairs)
    (expr; kw...) -> SymbolicUtils.substitute(expr, dict; kw...)
end

macro showarr(x)
    n = string(x)
    quote
        y = $(esc(x))
        println($n, " = ", summary(y))
        Base.print_array(stdout, y)
        println()
        y
    end
end

@deprecate substitute_expr!(expr,s) substitute(expr,s)

function states_to_sym(states)
    function _states_to_sym(O)
        if O isa Equation
            Expr(:(=), _states_to_sym(O.lhs), _states_to_sym(O.rhs))
        elseif O isa Term
            if isa(O.op, Sym)
                any(isequal(O), states) && return O.op.name  # dependent variables
                return build_expr(:call, Any[O.op.name; _states_to_sym.(O.args)])
            else
                return build_expr(:call, Any[O.op; _states_to_sym.(O.args)])
            end
        elseif O isa Num
            return _states_to_sym(value(O))
        else
            return toexpr(O)
        end
    end
end

"""
    toparam(s::Sym) -> Sym{<:Parameter}

Maps the variable to a paramter.
"""
toparam(s::Sym) = Sym{Parameter{symtype(s)}}(s.name)
toparam(s::Sym{<:Parameter}) = s

"""
    tovar(s::Sym) -> Sym{Real}

Maps the variable to a variable (state).
"""
tovar(s::Sym{<:Parameter}) = Sym{symtype(s)}(s.name)
tovar(s::Sym) = s
