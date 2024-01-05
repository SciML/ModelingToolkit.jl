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
    elseif issym(operation(op))
        Sym{Real}(nameof(operation(op)))
    else
        similarterm(op, operation(op), detime_dvs.(arguments(op)))
    end
end

function retime_dvs(op, dvs, iv)
    issym(op) && return Sym{FnType{Tuple{symtype(iv)}, Real}}(nameof(op))(iv)
    istree(op) ?
    similarterm(op, operation(op), retime_dvs.(arguments(op), (dvs,), (iv,))) :
    op
end

function modified_states!(mstates, e::Equation, statelist = nothing)
    get_variables!(mstates, e.lhs, statelist)
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

@deprecate substitute_expr!(expr, s) substitute(expr, s)

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

function rec_remove_macro_linenums!(expr)
    if expr isa Expr
        if expr.head === :macrocall
            expr.args[2] = nothing
            rec_remove_macro_linenums!(expr.args[3])
        else
            for ex in expr.args
                rec_remove_macro_linenums!(ex)
            end
        end
    end
    expr
end
function readable_code(expr)
    expr = Base.remove_linenums!(_readable_code(expr))
    rec_remove_macro_linenums!(expr)
    JuliaFormatter.format_text(string(expr), JuliaFormatter.SciMLStyle())
end

# System validation enums
const CheckNone = 0
const CheckAll = 1 << 0
const CheckComponents = 1 << 1
const CheckUnits = 1 << 2

function check_parameters(ps, iv)
    for p in ps
        isequal(iv, p) &&
            throw(ArgumentError("Independent variable $iv not allowed in parameters."))
        isparameter(p) ||
            throw(ArgumentError("$p is not a parameter."))
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
        isequal(iv, dv) &&
            throw(ArgumentError("Independent variable $iv not allowed in dependent variables."))
        (is_delay_var(iv, dv) || occursin(iv, dv)) ||
            throw(ArgumentError("Variable $dv is not a function of independent variable $iv."))
        isparameter(dv) &&
            throw(ArgumentError("$dv is not a state. It is a parameter."))
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
function collect_ivs(eqs, op = Differential)
    vars = Set()
    ivs = Set()
    for eq in eqs
        vars!(vars, eq; op = op)
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
    length(ivs) <= 1 ||
        throw(ArgumentError("Differential w.r.t. multiple variables $display are not allowed."))
    if length(ivs) == 1
        single_iv = pop!(ivs)
        isequal(single_iv, iv) ||
            throw(ArgumentError("Differential w.r.t. variable ($single_iv) other than the independent variable ($iv) are not allowed."))
    end
end
"""
Get all the independent variables with respect to which differentials/differences are taken.
"""
function collect_ivs_from_nested_operator!(ivs, x, target_op)
    if !istree(x)
        return
    end
    op = operation(unwrap(x))
    if op isa target_op
        push!(ivs, get_iv(op))
        x = if target_op <: Differential
            op.x
        elseif target_op <: Difference
            op.t
        else
            error("Unknown target op type in collect_ivs $target_op. Pass Difference or Differential")
        end
        collect_ivs_from_nested_operator!(ivs, x, target_op)
    end
end

function iv_from_nested_derivative(x, op = Differential)
    if istree(x) && operation(x) == getindex
        iv_from_nested_derivative(arguments(x)[1], op)
    elseif istree(x)
        operation(x) isa op ? iv_from_nested_derivative(arguments(x)[1], op) :
        arguments(x)[1]
    elseif issym(x)
        x
    else
        nothing
    end
end

hasdefault(v) = hasmetadata(v, Symbolics.VariableDefaultValue)
getdefault(v) = value(getmetadata(v, Symbolics.VariableDefaultValue))
function getdefaulttype(v)
    def = value(getmetadata(unwrap(v), Symbolics.VariableDefaultValue, nothing))
    def === nothing ? Float64 : typeof(def)
end
function setdefault(v, val)
    val === nothing ? v : wrap(setdefaultval(unwrap(v), value(val)))
end

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

"""
Throw error when difference/derivative operation occurs in the R.H.S.
"""
@noinline function throw_invalid_operator(opvar, eq, op::Type)
    if op === Difference
        optext = "difference"
    elseif op === Differential
        optext = "derivative"
    end
    msg = "The $optext variable must be isolated to the left-hand " *
          "side of the equation like `$opvar ~ ...`. You may want to use `structural_simplify` or the DAE form.\nGot $eq."
    throw(InvalidSystemException(msg))
end

"""
Check if difference/derivative operation occurs in the R.H.S. of an equation
"""
function _check_operator_variables(eq, op::T, expr = eq.rhs) where {T}
    istree(expr) || return nothing
    if operation(expr) isa op
        throw_invalid_operator(expr, eq, op)
    end
    foreach(expr -> _check_operator_variables(eq, op, expr),
        SymbolicUtils.unsorted_arguments(expr))
end
"""
Check if all the LHS are unique
"""
function check_operator_variables(eqs, op::T) where {T}
    ops = Set()
    tmp = Set()
    for eq in eqs
        _check_operator_variables(eq, op)
        vars!(tmp, eq.lhs)
        if length(tmp) == 1
            x = only(tmp)
            if op === Differential
                # Having a difference is fine for ODEs
                is_tmp_fine = isdifferential(x) || isdifference(x)
            else
                is_tmp_fine = istree(x) && !(operation(x) isa op)
            end
        else
            nd = count(x -> istree(x) && !(operation(x) isa op), tmp)
            is_tmp_fine = iszero(nd)
        end
        is_tmp_fine ||
            error("The LHS cannot contain nondifferentiated variables. Please run `structural_simplify` or use the DAE form.\nGot $eq")
        for v in tmp
            v in ops &&
                error("The LHS operator must be unique. Please run `structural_simplify` or use the DAE form. $v appears in LHS more than once.")
            push!(ops, v)
        end
        empty!(tmp)
    end
end

isoperator(expr, op) = istree(expr) && operation(expr) isa op
isoperator(op) = expr -> isoperator(expr, op)

isdifferential(expr) = isoperator(expr, Differential)
isdiffeq(eq) = isdifferential(eq.lhs)

isdifference(expr) = isoperator(expr, Difference)
isdifferenceeq(eq) = isdifference(eq.lhs)

function iv_from_nested_difference(x::Symbolic)
    istree(x) || return x
    operation(x) isa Difference ? iv_from_nested_difference(arguments(x)[1]) :
    arguments(x)[1]
end
iv_from_nested_difference(x) = nothing

var_from_nested_difference(x, i = 0) = (nothing, nothing)
function var_from_nested_difference(x::Symbolic, i = 0)
    istree(x) && operation(x) isa Difference ?
    var_from_nested_difference(arguments(x)[1], i + 1) :
    (x, i)
end

isvariable(x::Num)::Bool = isvariable(value(x))
function isvariable(x)::Bool
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
function vars(exprs::Symbolic; op = Differential)
    istree(exprs) ? vars([exprs]; op = op) : Set([exprs])
end
vars(exprs; op = Differential) = foldl((x, y) -> vars!(x, y; op = op), exprs; init = Set())
vars(eq::Equation; op = Differential) = vars!(Set(), eq; op = op)
function vars!(vars, eq::Equation; op = Differential)
    (vars!(vars, eq.lhs; op = op); vars!(vars, eq.rhs; op = op); vars)
end
function vars!(vars, O; op = Differential)
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
        vars!(vars, arg; op = op)
    end

    return vars
end

difference_vars(x) = vars(x; op = Difference)
difference_vars!(vars, O) = vars!(vars, O; op = Difference)

function collect_operator_variables(sys::AbstractSystem, args...)
    collect_operator_variables(equations(sys), args...)
end
function collect_operator_variables(eq::Equation, args...)
    collect_operator_variables([eq], args...)
end

"""
    collect_operator_variables(eqs::AbstractVector{Equation}, op)

Return a `Set` containing all variables that have Operator `op` applied to them.
See also [`collect_differential_variables`](@ref), [`collect_difference_variables`](@ref).
"""
function collect_operator_variables(eqs::AbstractVector{Equation}, op)
    vars = Set()
    diffvars = Set()
    for eq in eqs
        vars!(vars, eq; op = op)
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
    v = vars(x, op = op)
    filter(v) do x
        issym(x) && return false
        istree(x) && return operation(x) isa op
        false
    end
end

function find_derivatives!(vars, expr::Equation, f = identity)
    (find_derivatives!(vars, expr.lhs, f); find_derivatives!(vars, expr.rhs, f); vars)
end
function find_derivatives!(vars, expr, f)
    !istree(O) && return vars
    operation(O) isa Differential && push!(vars, f(O))
    for arg in arguments(O)
        vars!(vars, arg)
    end
    return vars
end

function collect_vars!(states, parameters, expr, iv)
    if issym(expr)
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
    if issym(expr)
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
    elseif !isconstant(var)
        push!(states, var)
    end
    # Add also any parameters that appear only as defaults in the var
    if hasdefault(var)
        collect_vars!(states, parameters, getdefault(var), iv)
    end
    return nothing
end

"""
Find all the symbolic constants of some equations or terms and return them as a vector.
"""
function collect_constants(x)
    constants = BasicSymbolic[]
    collect_constants!(constants, x)
    return constants
end

function collect_constants!(constants, arr::AbstractArray)
    for el in arr
        collect_constants!(constants, el)
    end
end

function collect_constants!(constants, eq::Equation)
    collect_constants!(constants, eq.lhs)
    collect_constants!(constants, eq.rhs)
end

function collect_constants!(constants, eq::Inequality)
    collect_constants!(constants, eq.lhs)
    collect_constants!(constants, eq.rhs)
end

collect_constants!(constants, x::Num) = collect_constants!(constants, unwrap(x))
collect_constants!(constants, x::Real) = nothing
collect_constants(n::Nothing) = BasicSymbolic[]

function collect_constants!(constants, expr::Symbolic)
    if issym(expr) && isconstant(expr)
        push!(constants, expr)
    else
        evars = vars(expr)
        if length(evars) == 1 && isequal(only(evars), expr)
            return nothing #avoid infinite recursion for vars(x(t)) == [x(t)]
        else
            for var in evars
                collect_constants!(constants, var)
            end
        end
    end
end

"""
Replace symbolic constants with their literal values
"""
function eliminate_constants(eqs, cs)
    cmap = Dict(x => getdefault(x) for x in cs)
    return substitute(eqs, cmap)
end

"""
Create a function preface containing assignments of default values to constants.
"""
function get_preprocess_constants(eqs)
    cs = collect_constants(eqs)
    pre = ex -> Let(Assignment[Assignment(x, getdefault(x)) for x in cs],
        ex, false)
    return pre
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
function find_duplicates(xs, ::Val{Ret} = Val(false)) where {Ret}
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

function get_cmap(sys)
    #Inject substitutions for constants => values
    cs = collect_constants([get_eqs(sys); get_observed(sys)]) #ctrls? what else?
    if !empty_substitutions(sys)
        cs = [cs; collect_constants(get_substitutions(sys).subs)]
    end
    # Swap constants for their values
    cmap = map(x -> x ~ getdefault(x), cs)
    return cmap, cs
end

function get_substitutions_and_solved_states(sys; no_postprocess = false)
    cmap, cs = get_cmap(sys)
    if empty_substitutions(sys) && isempty(cs)
        sol_states = Code.LazyState()
        pre = no_postprocess ? (ex -> ex) : get_postprocess_fbody(sys)
    else # Have to do some work
        if !empty_substitutions(sys)
            @unpack subs = get_substitutions(sys)
        else
            subs = []
        end
        subs = [cmap; subs] # The constants need to go first
        sol_states = Code.NameState(Dict(eq.lhs => Symbol(eq.lhs) for eq in subs))
        if no_postprocess
            pre = ex -> Let(Assignment[Assignment(eq.lhs, eq.rhs) for eq in subs], ex,
                false)
        else
            process = get_postprocess_fbody(sys)
            pre = ex -> Let(Assignment[Assignment(eq.lhs, eq.rhs) for eq in subs],
                process(ex), false)
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

@noinline function throw_missingvars_in_sys(vars)
    throw(ArgumentError("$vars are either missing from the variable map or missing from the system's states/parameters list."))
end

function promote_to_concrete(vs; tofloat = true, use_union = true)
    if isempty(vs)
        return vs
    end
    if vs isa Tuple #special rule, if vs is a Tuple, preserve types, container converted to Array
        tofloat = false
        use_union = true
        vs = Any[vs...]
    end
    T = eltype(vs)
    if Base.isconcretetype(T) && (!tofloat || T === float(T)) # nothing to do
        return vs
    else
        sym_vs = filter(x -> SymbolicUtils.issym(x) || SymbolicUtils.istree(x), vs)
        isempty(sym_vs) || throw_missingvars_in_sys(sym_vs)

        C = nothing
        for v in vs
            E = typeof(v)
            if E <: Number
                if tofloat
                    E = float(E)
                end
            end
            if C === nothing
                C = E
            end
            if use_union
                C = Union{C, E}
            else
                @assert C==E "`promote_to_concrete` can't make type $E uniform with $C"
                C = E
            end
        end

        y = similar(vs, C)
        for i in eachindex(vs)
            if (vs[i] isa Number) & tofloat
                y[i] = float(vs[i]) #needed because copyto! can't convert Int to Float automatically
            else
                y[i] = vs[i]
            end
        end

        return y
    end
end

struct BitDict <: AbstractDict{Int, Int}
    keys::Vector{Int}
    values::Vector{Union{Nothing, Int}}
end
BitDict(n::Integer) = BitDict(Int[], Union{Nothing, Int}[nothing for _ in 1:n])
struct BitDictKeySet <: AbstractSet{Int}
    d::BitDict
end

Base.keys(d::BitDict) = BitDictKeySet(d)
Base.in(v::Integer, s::BitDictKeySet) = s.d.values[v] !== nothing
Base.iterate(s::BitDictKeySet, state...) = iterate(s.d.keys, state...)
function Base.setindex!(d::BitDict, val::Integer, ind::Integer)
    if 1 <= ind <= length(d.values) && d.values[ind] === nothing
        push!(d.keys, ind)
    end
    d.values[ind] = val
end
function Base.getindex(d::BitDict, ind::Integer)
    if 1 <= ind <= length(d.values) && d.values[ind] === nothing
        return d.values[ind]
    else
        throw(KeyError(ind))
    end
end
function Base.iterate(d::BitDict, state...)
    r = Base.iterate(d.keys, state...)
    r === nothing && return nothing
    k, state = r
    (k => d.values[k]), state
end
function Base.empty!(d::BitDict)
    for v in d.keys
        d.values[v] = nothing
    end
    empty!(d.keys)
    d
end

abstract type AbstractSimpleTreeIter{T} end
Base.IteratorSize(::Type{<:AbstractSimpleTreeIter}) = Base.SizeUnknown()
Base.eltype(::Type{<:AbstractSimpleTreeIter{T}}) where {T} = childtype(T)
has_fast_reverse(::Type{<:AbstractSimpleTreeIter}) = true
has_fast_reverse(::T) where {T <: AbstractSimpleTreeIter} = has_fast_reverse(T)
reverse_buffer(it::AbstractSimpleTreeIter) = has_fast_reverse(it) ? nothing : eltype(it)[]
reverse_children!(::Nothing, cs) = Iterators.reverse(cs)
function reverse_children!(rev_buff, cs)
    Iterators.reverse(cs)
    empty!(rev_buff)
    for c in cs
        push!(rev_buff, c)
    end
    Iterators.reverse(rev_buff)
end

struct StatefulPreOrderDFS{T} <: AbstractSimpleTreeIter{T}
    t::T
end
function Base.iterate(it::StatefulPreOrderDFS,
        state = (eltype(it)[it.t], reverse_buffer(it)))
    stack, rev_buff = state
    isempty(stack) && return nothing
    t = pop!(stack)
    for c in reverse_children!(rev_buff, children(t))
        push!(stack, c)
    end
    return t, state
end
struct StatefulPostOrderDFS{T} <: AbstractSimpleTreeIter{T}
    t::T
end
function Base.iterate(it::StatefulPostOrderDFS,
        state = (eltype(it)[it.t], falses(1), reverse_buffer(it)))
    isempty(state[2]) && return nothing
    vstack, sstack, rev_buff = state
    while true
        t = pop!(vstack)
        isresume = pop!(sstack)
        isresume && return t, state
        push!(vstack, t)
        push!(sstack, true)
        for c in reverse_children!(rev_buff, children(t))
            push!(vstack, c)
            push!(sstack, false)
        end
    end
end

# Note that StatefulBFS also returns the depth.
struct StatefulBFS{T} <: AbstractSimpleTreeIter{T}
    t::T
end
Base.eltype(::Type{<:StatefulBFS{T}}) where {T} = Tuple{Int, childtype(T)}
function Base.iterate(it::StatefulBFS, queue = (eltype(it)[(0, it.t)]))
    isempty(queue) && return nothing
    lv, t = popfirst!(queue)
    nextlv = lv + 1
    for c in children(t)
        push!(queue, (nextlv, c))
    end
    return (lv, t), queue
end

function jacobian_wrt_vars(pf::F, p, input_idxs, chunk::C) where {F, C}
    E = eltype(p)
    tag = ForwardDiff.Tag(pf, E)
    T = typeof(tag)
    dualtype = ForwardDiff.Dual{T, E, ForwardDiff.chunksize(chunk)}
    p_big = similar(p, dualtype)
    copyto!(p_big, p)
    p_closure = let pf = pf,
        input_idxs = input_idxs,
        p_big = p_big

        function (p_small_inner)
            p_big[input_idxs] .= p_small_inner
            pf(p_big)
        end
    end
    p_small = p[input_idxs]
    cfg = ForwardDiff.JacobianConfig(p_closure, p_small, chunk, tag)
    ForwardDiff.jacobian(p_closure, p_small, cfg, Val(false))
end

function fold_constants(ex)
    if istree(ex)
        similarterm(ex, operation(ex), map(fold_constants, arguments(ex)),
            symtype(ex); metadata = metadata(ex))
    elseif issym(ex) && isconstant(ex)
        getdefault(ex)
    else
        ex
    end
end

# Symbolics needs to call unwrap on the substitution rules, but most of the time
# we don't want to do that in MTK.
const Eq = Union{Equation, Inequality}
function fast_substitute(eq::Eq, subs)
    if eq isa Inequality
        Inequality(fast_substitute(eq.lhs, subs), fast_substitute(eq.rhs, subs),
            eq.relational_op)
    else
        Equation(fast_substitute(eq.lhs, subs), fast_substitute(eq.rhs, subs))
    end
end
function fast_substitute(eq::T, subs::Pair) where {T <: Eq}
    T(fast_substitute(eq.lhs, subs), fast_substitute(eq.rhs, subs))
end
fast_substitute(eqs::AbstractArray{<:Eq}, subs) = fast_substitute.(eqs, (subs,))
fast_substitute(a, b) = substitute(a, b)
function fast_substitute(expr, pair::Pair)
    a, b = pair
    isequal(expr, a) && return b

    istree(expr) || return expr
    op = fast_substitute(operation(expr), pair)
    canfold = Ref(!(op isa Symbolic))
    args = let canfold = canfold
        map(SymbolicUtils.unsorted_arguments(expr)) do x
            x′ = fast_substitute(x, pair)
            canfold[] = canfold[] && !(x′ isa Symbolic)
            x′
        end
    end
    canfold[] && return op(args...)

    similarterm(expr,
        op,
        args,
        symtype(expr);
        metadata = metadata(expr))
end

normalize_to_differential(s) = s
