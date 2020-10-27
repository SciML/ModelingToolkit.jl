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
substitute(expr, s::Pair)
substitute(expr, s::Dict)
substitute(expr, s::Vector)

Performs the substitution `Num => val` on the `expr` Num.
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

function states_to_sym(states::Set)
    function _states_to_sym(O)
        if O isa Equation
            Expr(:(=), _states_to_sym(O.lhs), _states_to_sym(O.rhs))
        elseif O isa Term
            if isa(O.op, Sym)
                O in states && return tosymbol(O)
                # dependent variables
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
states_to_sym(states) = states_to_sym(Set(states))

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

Base.Symbol(x::Union{Num,Symbolic}) = tosymbol(x)
tosymbol(x; kwargs...) = x
tosymbol(x::Sym; kwargs...) = nameof(x)
tosymbol(t::Num; kwargs...) = tosymbol(value(t); kwargs...)

"""
    tosymbol(x::Union{Num,Symbolic}; states=nothing, escape=true) -> Symbol

Convert `x` to a symbol. `states` are the states of a system, and `escape`
means if the target has escapes like `val"y⦗t⦘"`. If `escape` then it will only
output `y` instead of `y⦗t⦘`.
"""
function tosymbol(t::Term; states=nothing, escape=true)
    if t.op isa Sym
        if states !== nothing && !(any(isequal(t), states))
            return nameof(t.op)
        end
        op = nameof(t.op)
        args = t.args
    elseif t.op isa Differential
        term = diff2term(t)
        op = Symbol(operation(term))
        args = arguments(term)
    else
        @goto err
    end

    return escape ? Symbol(op, "⦗", join(args, ", "), "⦘") : op
    @label err
    error("Cannot convert $t to a symbol")
end

makesym(t::Symbolic; kwargs...) = Sym{symtype(t)}(tosymbol(t; kwargs...))
makesym(t::Num; kwargs...) = makesym(value(t); kwargs...)

function lower_varname(var::Term, idv, order)
    order == 0 && return var
    name = string(nameof(var.op))
    underscore = 'ˍ'
    idx = findlast(underscore, name)
    append = string(idv)^order
    if idx === nothing
        newname = Symbol(name, underscore, append)
    else
        nidx = nextind(name, idx)
        newname = Symbol(name[1:idx], name[nidx:end], append)
    end
    return Sym{symtype(var.op)}(newname)(var.args[1])
end

function lower_varname(t::Term, iv)
    var, order = var_from_nested_derivative(t)
    lower_varname(var, iv, order)
end
lower_varname(t::Sym, iv) = t

function flatten_differential(O::Term)
    @assert is_derivative(O) "invalid differential: $O"
    is_derivative(O.args[1]) || return (O.args[1], O.op.x, 1)
    (x, t, order) = flatten_differential(O.args[1])
    isequal(t, O.op.x) || throw(ArgumentError("non-matching differentials on lhs: $t, $(O.op.x)"))
    return (x, t, order + 1)
end

"""
    diff2term(x::Term) -> Term
    diff2term(x) -> x

diff2term(D(D(x(t)))) -> xˍtt(t)
"""
function diff2term(O)
    isa(O, Term) || return O
    if is_derivative(O)
        (x, t, order) = flatten_differential(O)
        return lower_varname(x, t, order)
    end
    return Term(O.op, diff2term.(O.args))
end
