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

function detime_dvs(op)
    if !istree(op)
        op
    elseif operation(op) isa Sym
        Sym{Real}(nameof(operation(op)))
    else
        similarterm(op, operation(op),detime_dvs.(arguments(op)))
    end
end

function retime_dvs(op::Sym,dvs,iv)
    Sym{FnType{Tuple{symtype(iv)}, Real}}(nameof(op))(iv)
end

function retime_dvs(op, dvs, iv)
    istree(op) ?
        similarterm(op, operation(op), retime_dvs.(arguments(op),(dvs,),(iv,))) :
        op
end

modified_states!(mstates, e::Equation, statelist=nothing) = get_variables!(mstates, e.lhs, statelist)

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
        elseif istree(O)
            op = operation(O)
            args = arguments(O)
            if op isa Sym
                O in states && return tosymbol(O)
                # dependent variables
                return build_expr(:call, Any[nameof(op); _states_to_sym.(args)])
            else
                canonical, O = canonicalexpr(O)
                return canonical ? O : build_expr(:call, Any[op; _states_to_sym.(args)])
            end
        elseif O isa Num
            return _states_to_sym(value(O))
        else
            return toexpr(O)
        end
    end
end
states_to_sym(states) = states_to_sym(Set(states))

function todict(d)
    eltype(d) <: Pair || throw(ArgumentError("The variable-value mapping must be a Dict."))
    d isa Dict ? d : Dict(d)
end

_merge(d1, d2) = merge(todict(d1), todict(d2))

function indepvar2depvar(s::Sym, args...)
    T = FnType{NTuple{length(args)}, symtype(s)}
    ns = Sym{T}(nameof(s))(args...)
    @set! ns.metadata = s.metadata
end

function _readable_code(ex)
    ex isa Expr || return ex
    if ex.head === :call
        f, args = ex.args[1], ex.args[2:end]
        if f isa Function && (nf = nameof(f); Base.isoperator(nf))
            expr = Expr(:call, nf)
            for a in args
                push!(expr.args, _readable_code(a))
            end
            return expr
        end
    end
    expr = Expr(ex.head)
    for a in ex.args
        push!(expr.args, _readable_code(a))
    end
    expr
end
readable_code(expr) = JuliaFormatter.format_text(string(Base.remove_linenums!(_readable_code(expr))))

function check_parameters(ps, iv)
    for p in ps
        isequal(iv, p) && throw(ArgumentError("Independent variable $iv not allowed in parameters."))
    end
end

function check_variables(dvs, iv)
    for dv in dvs
        isequal(iv, dv) && throw(ArgumentError("Independent variable $iv not allowed in dependent variables."))
        isequal(iv, iv_from_nested_derivative(dv)) || throw(ArgumentError("Variable $dv is not a function of independent variable $iv."))
    end
end

"Get all the independent variables with respect to which differentials are taken."
function collect_differentials(eqs)
    vars = Set()
    ivs = Set()
    for eq in eqs
        vars!(vars, eq)
        for v in vars
            isdifferential(v) || continue
            collect_ivs_from_nested_differential!(ivs, v)
        end
        empty!(vars)
    end
    return ivs
end

"Assert that equations are well-formed when building ODE."
function check_equations(eqs, iv)
    ivs = collect_differentials(eqs)
    display = collect(ivs)
    length(ivs) <= 1 || throw(ArgumentError("Differential w.r.t. multiple variables $display are not allowed."))
    if length(ivs) == 1
        single_iv = pop!(ivs)
        isequal(single_iv, iv) || throw(ArgumentError("Differential w.r.t. variable ($single_iv) other than the independent variable ($iv) are not allowed."))
    end
end
"Get all the independent variables with respect to which differentials are taken."
function collect_ivs_from_nested_differential!(ivs, x::Term)
    op = operation(x)
    if op isa Differential
        push!(ivs, op.x)
        collect_ivs_from_nested_differential!(ivs, arguments(x)[1])
    end
end

iv_from_nested_derivative(x::Term) = operation(x) isa Differential ? iv_from_nested_derivative(arguments(x)[1]) : arguments(x)[1]
iv_from_nested_derivative(x::Sym) = x
iv_from_nested_derivative(x) = missing

hasdefault(v) = hasmetadata(v, Symbolics.VariableDefaultValue)
getdefault(v) = value(getmetadata(v, Symbolics.VariableDefaultValue))
setdefault(v, val) = val === nothing ? v : setmetadata(v, Symbolics.VariableDefaultValue, value(val))

function process_variables!(array_vars, defs, vars)
    collect_defaults!(defs, vars)
    collect_array_vars!(array_vars, vars)
    return nothing
end

function collect_defaults!(defs, vars)
    for v in vars; (haskey(defs, v) || !hasdefault(v)) && continue
        defs[v] = getdefault(v)
    end
    return defs
end

function collect_array_vars!(array_vars, vars)
    for v in vars
        x = v = value(v)
        while istree(v)
            op = operation(v)
            op === getindex && (v = arguments(v)[1]; break)
            v = op
        end

        if symtype(v) <: AbstractArray
            if v isa Symbolics.ArrayOp
                n = nameof(arguments(operation(arguments(x)[2]))[1])
            else
                n = nameof(v)
            end
            if !hasmetadata(x, Symbolics.GetindexParent)
                array_vars[n] = arguments(x)[1]
            else
                array_vars[n] = getmetadata(x, Symbolics.GetindexParent)
            end
        end
    end
    return array_vars
end

"Throw error when difference/derivative operation occurs in the R.H.S."
@noinline function throw_invalid_operator(opvar, eq, op::Type)
    if op === Difference
        optext = "difference"
    elseif op === Differential
        optext="derivative"
    end
    msg = "The $optext variable must be isolated to the left-hand " *
    "side of the equation like `$opvar ~ ...`.\n Got $eq."
    throw(InvalidSystemException(msg))
end

"Check if difference/derivative operation occurs in the R.H.S. of an equation"
function check_operator_variables(eq, op::Type, expr=eq.rhs)
    istree(expr) || return nothing
    if operation(expr) isa op
        throw_invalid_operator(expr, eq, op)
    end
    foreach(expr -> check_operator_variables(eq, op, expr), arguments(expr))
end
