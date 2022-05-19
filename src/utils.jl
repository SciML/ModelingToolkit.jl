get_iv(D::Differential) = D.x
get_iv(D::Difference) = D.t

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
        similarterm(op, operation(op), detime_dvs.(arguments(op)))
    end
end

function retime_dvs(op::Sym, dvs, iv)
    Sym{FnType{Tuple{symtype(iv)},Real}}(nameof(op))(iv)
end

function retime_dvs(op, dvs, iv)
    istree(op) ?
    similarterm(op, operation(op), retime_dvs.(arguments(op), (dvs,), (iv,))) :
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

@deprecate substitute_expr!(expr, s) substitute(expr, s)

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

function is_delay_var(iv, var)
    args = nothing
    try
        args = arguments(var)
    catch
        return false
    end
    length(args) > 1 && return false
    isequal(first(args), iv) && return false
    delay = iv - first(args)
    delay isa Integer ||
        delay isa AbstractFloat ||
        (delay isa Num && isreal(value(delay)))
end

function check_variables(dvs, iv)
    for dv in dvs
        isequal(iv, dv) && throw(ArgumentError("Independent variable $iv not allowed in dependent variables."))
        (is_delay_var(iv, dv) || occursin(iv, iv_from_nested_derivative(dv))) || throw(ArgumentError("Variable $dv is not a function of independent variable $iv."))
    end
end

function check_lhs(eq::Equation, op, dvs::Set)
    v = unwrap(eq.lhs)
    _iszero(v) && return
    (operation(v) isa op && only(arguments(v)) in dvs) && return
    error("$v is not a valid LHS. Please run structural_simplify before simulation.")
end
check_lhs(eqs, op, dvs::Set) =
    for eq in eqs
        check_lhs(eq, op, dvs)
    end

"""
    collect_ivs(eqs, op = Differential)

Get all the independent variables with respect to which differentials (`op`) are taken.
"""
function collect_ivs(eqs, op=Differential)
    vars = Set()
    ivs = Set()
    for eq in eqs
        vars!(vars, eq; op=op)
        for v in vars
            if isoperator(v, op)
                collect_ivs_from_nested_operator!(ivs, v, op)
            end
        end
        empty!(vars)
    end
    return ivs
end

"""
    check_equations(eqs, iv)

Assert that equations are well-formed when building ODE, i.e., only containing a single independent variable.
"""
function check_equations(eqs, iv)
    ivs = collect_ivs(eqs)
    display = collect(ivs)
    length(ivs) <= 1 || throw(ArgumentError("Differential w.r.t. multiple variables $display are not allowed."))
    if length(ivs) == 1
        single_iv = pop!(ivs)
        isequal(single_iv, iv) || throw(ArgumentError("Differential w.r.t. variable ($single_iv) other than the independent variable ($iv) are not allowed."))
    end
end
"Get all the independent variables with respect to which differentials/differences are taken."
function collect_ivs_from_nested_operator!(ivs, x::Term, target_op)
    op = operation(x)
    if op isa target_op
        push!(ivs, get_iv(op))
        collect_ivs_from_nested_operator!(ivs, arguments(x)[1], target_op)
    end
end

iv_from_nested_derivative(x::Term, op=Differential) = operation(x) isa op ? iv_from_nested_derivative(arguments(x)[1], op) : arguments(x)[1]
iv_from_nested_derivative(x::Sym, op=Differential) = x
iv_from_nested_derivative(x, op=Differential) = nothing

hasdefault(v) = hasmetadata(v, Symbolics.VariableDefaultValue)
getdefault(v) = value(getmetadata(v, Symbolics.VariableDefaultValue))
setdefault(v, val) = val === nothing ? v : setmetadata(v, Symbolics.VariableDefaultValue, value(val))

function process_variables!(var_to_name, defs, vars)
    collect_defaults!(defs, vars)
    collect_var_to_name!(var_to_name, vars)
    return nothing
end

function collect_defaults!(defs, vars)
    for v in vars
        (haskey(defs, v) || !hasdefault(v)) && continue
        defs[v] = getdefault(v)
    end
    return defs
end

function collect_var_to_name!(vars, xs)
    for x in xs
        x = unwrap(x)
        if hasmetadata(x, Symbolics.GetindexParent)
            xarr = getmetadata(x, Symbolics.GetindexParent)
            vars[Symbolics.getname(xarr)] = xarr
        else
            if istree(x) && operation(x) === getindex
                x = arguments(x)[1]
            end
            vars[Symbolics.getname(unwrap(x))] = x
        end
    end

end

"Throw error when difference/derivative operation occurs in the R.H.S."
@noinline function throw_invalid_operator(opvar, eq, op::Type)
    if op === Difference
        optext = "difference"
    elseif op === Differential
        optext = "derivative"
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
    foreach(expr -> check_operator_variables(eq, op, expr), SymbolicUtils.unsorted_arguments(expr))
end

isoperator(expr, op) = istree(expr) && operation(expr) isa op
isoperator(op) = expr -> isoperator(expr, op)

isdifferential(expr) = isoperator(expr, Differential)
isdiffeq(eq) = isdifferential(eq.lhs)

isdifference(expr) = isoperator(expr, Difference)
isdifferenceeq(eq) = isdifference(eq.lhs)

iv_from_nested_difference(x::Term) = operation(x) isa Difference ? iv_from_nested_difference(arguments(x)[1]) : arguments(x)[1]
iv_from_nested_difference(x::Sym) = x
iv_from_nested_difference(x) = nothing

var_from_nested_difference(x, i=0) = (nothing, nothing)
var_from_nested_difference(x::Term, i=0) = operation(x) isa Difference ? var_from_nested_difference(arguments(x)[1], i + 1) : (x, i)
var_from_nested_difference(x::Sym, i=0) = (x, i)


isvariable(x::Num) = isvariable(value(x))
function isvariable(x)
    x isa Symbolic || return false
    p = getparent(x, nothing)
    p === nothing || (x = p)
    hasmetadata(x, VariableSource)
end

"""
    vars(x; op=Differential)
Return a `Set` containing all variables in `x` that appear in
- differential equations if `op = Differential`
- difference equations if `op = Differential`

Example:
```
@variables t u(t) y(t)
D  = Differential(t)
v  = ModelingToolkit.vars(D(y) ~ u)
v == Set([D(y), u])
```
"""
vars(x::Sym; op=Differential) = Set([x])
vars(exprs::Symbolic; op=Differential) = vars([exprs]; op=op)
vars(exprs; op=Differential) = foldl((x, y) -> vars!(x, y; op=op), exprs; init=Set())
vars(eq::Equation; op=Differential) = vars!(Set(), eq; op=op)
vars!(vars, eq::Equation; op=Differential) = (vars!(vars, eq.lhs; op=op); vars!(vars, eq.rhs; op=op); vars)
function vars!(vars, O; op=Differential)
    if isvariable(O)
        return push!(vars, O)
    end
    !istree(O) && return vars

    operation(O) isa op && return push!(vars, O)

    if operation(O) === (getindex) &&
       isvariable(first(arguments(O)))

        return push!(vars, O)
    end

    isvariable(operation(O)) && push!(vars, O)
    for arg in arguments(O)
        vars!(vars, arg; op=op)
    end

    return vars
end

difference_vars(x) = vars(x; op=Difference)
difference_vars!(vars, O) = vars!(vars, O; op=Difference)

collect_operator_variables(sys::AbstractSystem, args...) = collect_operator_variables(equations(sys), args...)
collect_operator_variables(eq::Equation, args...) = collect_operator_variables([eq], args...)

"""
    collect_operator_variables(eqs::AbstractVector{Equation}, op)

Return a `Set` containing all variables that have Operator `op` applied to them.
See also [`collect_differential_variables`](@ref), [`collect_difference_variables`](@ref).
"""
function collect_operator_variables(eqs::AbstractVector{Equation}, op)
    vars = Set()
    diffvars = Set()
    for eq in eqs
        vars!(vars, eq; op=op)
        for v in vars
            isoperator(v, op) || continue
            push!(diffvars, arguments(v)[1])
        end
        empty!(vars)
    end
    return diffvars
end
collect_differential_variables(sys) = collect_operator_variables(sys, Differential)
collect_difference_variables(sys) = collect_operator_variables(sys, Difference)

"""
    collect_applied_operators(x, op)

Return  a `Set` with all applied operators in `x`, example:
```
@variables t u(t) y(t)
D = Differential(t)
eq = D(y) ~ u
ModelingToolkit.collect_applied_operators(eq, Differential) == Set([D(y)])
```
The difference compared to `collect_operator_variables` is that `collect_operator_variables` returns the variable without the operator applied.
"""
function collect_applied_operators(x, op)
    v = vars(x, op=op)
    filter(v) do x
        x isa Sym && return false
        istree(x) && return operation(x) isa op
        false
    end
end

find_derivatives!(vars, expr::Equation, f=identity) = (find_derivatives!(vars, expr.lhs, f); find_derivatives!(vars, expr.rhs, f); vars)
function find_derivatives!(vars, expr, f)
    !istree(O) && return vars
    operation(O) isa Differential && push!(vars, f(O))
    for arg in arguments(O)
        vars!(vars, arg)
    end
    return vars
end

function collect_vars!(states, parameters, expr, iv)
    if expr isa Sym
        collect_var!(states, parameters, expr, iv)
    else
        for var in vars(expr)
            if istree(var) && operation(var) isa Differential
                var, _ = var_from_nested_derivative(var)
            end
            collect_var!(states, parameters, var, iv)
        end
    end
    return nothing
end

function collect_vars_difference!(states, parameters, expr, iv)
    if expr isa Sym
        collect_var!(states, parameters, expr, iv)
    else
        for var in vars(expr)
            if istree(var) && operation(var) isa Difference
                var, _ = var_from_nested_difference(var)
            end
            collect_var!(states, parameters, var, iv)
        end
    end
    return nothing
end

function collect_var!(states, parameters, var, iv)
    isequal(var, iv) && return nothing
    if isparameter(var) || (istree(var) && isparameter(operation(var)))
        push!(parameters, var)
    else
        push!(states, var)
    end
    return nothing
end


function get_postprocess_fbody(sys)
    if has_preface(sys) && (pre = preface(sys); pre !== nothing)
        pre_ = let pre = pre
            ex -> Let(pre, ex, false)
        end
    else
        pre_ = ex -> ex
    end
    return pre_
end

"""
$(SIGNATURES)

find duplicates in an iterable object.
"""
function find_duplicates(xs, ::Val{Ret}=Val(false)) where {Ret}
    appeared = Set()
    duplicates = Set()
    for x in xs
        if x in appeared
            push!(duplicates, x)
        else
            push!(appeared, x)
        end
    end
    return Ret ? (appeared, duplicates) : duplicates
end

isarray(x) = x isa AbstractArray || x isa Symbolics.Arr

function empty_substitutions(sys)
    has_substitutions(sys) || return true
    subs = get_substitutions(sys)
    isnothing(subs) || isempty(subs.deps)
end

function get_substitutions_and_solved_states(sys; no_postprocess=false)
    if empty_substitutions(sys)
        sol_states = Code.LazyState()
        pre = no_postprocess ? (ex -> ex) : get_postprocess_fbody(sys)
    else
        @unpack subs = get_substitutions(sys)
        sol_states = Code.NameState(Dict(eq.lhs => Symbol(eq.lhs) for eq in subs))
        if no_postprocess
            pre = ex -> Let(Assignment[Assignment(eq.lhs, eq.rhs) for eq in subs], ex, false)
        else
            process = get_postprocess_fbody(sys)
            pre = ex -> Let(Assignment[Assignment(eq.lhs, eq.rhs) for eq in subs], process(ex), false)
        end
    end
    return pre, sol_states
end

function mergedefaults(defaults, varmap, vars)
    defs = if varmap isa Dict
        merge(defaults, varmap)
    elseif eltype(varmap) <: Pair
        merge(defaults, Dict(varmap))
    elseif eltype(varmap) <: Number
        merge(defaults, Dict(zip(vars, varmap)))
    else
        defaults
    end
end

function promote_to_concrete(vs; tofloat=true, use_union=false)
    if isempty(vs)
        return vs
    end
    T = eltype(vs)
    if Base.isconcretetype(T) && (!tofloat || T === float(T)) # nothing to do
        vs
    else
        C = typeof(first(vs))
        I = Int8
        has_int = false
        has_array = false
        array_T = nothing
        for v in vs
            if v isa AbstractArray
                has_array = true
                array_T = typeof(v)
            end
            E = eltype(v)
            C = promote_type(C, E)
            if E <: Integer
                has_int = true
                I = promote_type(I, E)
            end
        end
        if tofloat && !has_array
            C = float(C)
        elseif has_array || (use_union && has_int && C !== I)
            if has_array
                C = Union{C,array_T}
            end
            if has_int
                C = Union{C,I}
            end
            return copyto!(similar(vs, C), vs)
        end
        convert.(C, vs)
    end
end

"""
Counts the Differential operators for given variable x. This is used to determine
the order of a PDE.
"""
function count_differentials(term, x::Symbolics.Symbolic)
    S = Symbolics
    SU = SymbolicUtils
    if !S.istree(term)
        return 0
    else
        op = SU.operation(term)
        count_children = sum(map(arg -> count_differentials(arg, x), SU.arguments(term)))
        if op isa Differential && isequal(op.x, x)
            return 1 + count_children
        end
        return count_children
    end
end

"""
return list of differential orders in the equation
"""
function differential_order(eq, x::Symbolics.Symbolic)
    S = Symbolics
    SU = SymbolicUtils
    orders = Set{Int}()
    if S.istree(eq)
        op = SU.operation(eq)
        if op isa Differential
            push!(orders, count_differentials(eq, x))
        else
            for o in map(ch -> differential_order(ch, x), SU.arguments(eq))
                union!(orders, o)
            end
        end
    end
    return filter(!iszero, orders)
end

"""
find all the dependent variables given by depvar_ops in an expression
"""
function get_depvars(eq, depvar_ops)
    S = Symbolics
    SU = SymbolicUtils
    depvars = Set()
    if eq isa Num
        eq = eq.val
    end
    if S.istree(eq)
        if eq isa Term && any(u -> isequal(operation(eq), u), depvar_ops)
            push!(depvars, eq)
        else
            for o in map(x -> get_depvars(x, depvar_ops), SU.arguments(eq))
                union!(depvars, o)
            end
        end
    end
    return depvars
end

@inline function get_all_depvars(pdesys, depvar_ops)
    pdeeqs = pdesys.eqs # Vector
    return collect(mapreduce(x -> get_depvars(x.lhs, depvar_ops), union, pdeeqs) ∪ mapreduce(x -> get_depvars(x.rhs, depvar_ops), union, pdeeqs))
end
