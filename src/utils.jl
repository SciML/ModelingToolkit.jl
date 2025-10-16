get_iv(D::Differential) = D.x
""""""
const CheckNone = 0
""""""
const CheckAll = 1 << 0
""""""
const CheckComponents = 1 << 1
""""""
const CheckUnits = 1 << 2
function check_independent_variables(ivs)
    for iv in ivs
        isparameter(iv) || @invokelatest warn_indepvar(iv)
    end
end
@noinline function warn_indepvar(iv::SymbolicT)
    @warn "Independent variable $iv should be defined with @independent_variables $iv."
end
function check_parameters(ps, iv)
    for p in ps
        isequal(iv, p) &&
            throw(ArgumentError("Independent variable $iv not allowed in parameters."))
    end
end
function is_delay_var(iv::SymbolicT, var::SymbolicT)
    Moshi.Match.@match var begin
        BSImpl.Term(; f, args) => begin
            length(args) > 1 && return false
            arg = args[1]
            isequal(arg, iv) && return false
            return symtype(arg) <: Real
        end
        _ => false
    end
end
function check_variables(dvs, iv)
    for dv in dvs
        isequal(iv, dv) &&
            throw(ArgumentError("Independent variable $iv not allowed in dependent variables."))
        (is_delay_var(iv, dv) || SU.query(isequal(iv), dv)) ||
            throw(ArgumentError("Variable $dv is not a function of independent variable $iv."))
    end
end
function check_lhs(eq::Equation, op, dvs::Set)
    v = unwrap(eq.lhs)
    _iszero(v) && return
    (operation(v) isa op && only(arguments(v)) in dvs) && return
    error("$v is not a valid LHS. Please run mtkcompile before simulation.")
end
check_lhs(eqs, op, dvs::Set) =
    for eq in eqs
        check_lhs(eq, op, dvs)
    end
""""""
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
struct IndepvarCheckPredicate
    iv::SymbolicT
end
function (icp::IndepvarCheckPredicate)(ex::SymbolicT)
    Moshi.Match.@match ex begin
        BSImpl.Term(; f) && if f isa Differential end => begin
            f = f::Differential
            isequal(f.x, icp.iv) || throw_multiple_iv(icp.iv, f.x)
            return false
        end
        _ => false
    end
end
@noinline function throw_multiple_iv(iv, newiv)
    throw(ArgumentError("Differential w.r.t. variable ($newiv) other than the independent variable ($iv) are not allowed."))
end
""""""
function check_equations(eqs::Vector{Equation}, iv::SymbolicT)
    icp = IndepvarCheckPredicate(iv)
    for eq in eqs
        SU.query(icp, eq.lhs)
        SU.query(icp, eq.rhs)
    end
end
""""""
function check_subsystems(systems)
    idxs = findall(!does_namespacing, systems)
    isempty(idxs) || throw_bad_namespacing(systems, idxs)
end
@noinline function throw_bad_namespacing(systems, idxs)
    names = join("  " .* string.(nameof.(systems[idxs])), "\n")
    throw(ArgumentError("All subsystems must have namespacing enabled. The following subsystems do not perform namespacing:\n$(names)"))
end
""""""
function collect_ivs_from_nested_operator!(ivs, x, target_op)
    if !iscall(x)
        return
    end
    op = operation(unwrap(x))
    if op isa target_op
        push!(ivs, get_iv(op))
        x = if target_op <: Differential
            op.x
        else
            error("Unknown target op type in collect_ivs $target_op. Pass Differential")
        end
        collect_ivs_from_nested_operator!(ivs, x, target_op)
    end
end
function iv_from_nested_derivative(x, op = Differential)
    if iscall(x) &&
       (operation(x) == getindex || operation(x) == real || operation(x) == imag)
        iv_from_nested_derivative(arguments(x)[1], op)
    elseif iscall(x)
        operation(x) isa op ? iv_from_nested_derivative(arguments(x)[1], op) :
        arguments(x)[1]
    elseif issym(x)
        x
    else
        nothing
    end
end
""""""
hasdefault(v) = hasmetadata(v, Symbolics.VariableDefaultValue)
""""""
getdefault(v) = value(Symbolics.getdefaultval(v))
""""""
function setdefault(v, val)
    val === nothing ? v : wrap(setdefaultval(unwrap(v), value(val)))
end
function process_variables!(var_to_name::Dict{Symbol, SymbolicT}, defs::SymmapT, guesses::SymmapT, vars::Vector{SymbolicT})
    collect_defaults!(defs, vars)
    collect_guesses!(guesses, vars)
    collect_var_to_name!(var_to_name, vars)
    return nothing
end
function process_variables!(var_to_name::Dict{Symbol, SymbolicT}, defs::SymmapT, vars::Vector{SymbolicT})
    collect_defaults!(defs, vars)
    collect_var_to_name!(var_to_name, vars)
    return nothing
end
function collect_defaults!(defs::SymmapT, vars::Vector{SymbolicT})
    for v in vars
        isconst(v) && continue
        haskey(defs, v) && continue
        def = Symbolics.getdefaultval(v, nothing)
        if def !== nothing
            defs[v] = SU.Const{VartypeT}(def)
            continue
        end
        Moshi.Match.@match v begin
            BSImpl.Term(; f, args) && if f === getindex end => begin
                haskey(defs, args[1]) && continue
                def = Symbolics.getdefaultval(args[1], nothing)
                def === nothing && continue
                defs[args[1]] = def
            end
            _ => nothing
        end
    end
    return defs
end
function collect_guesses!(guesses::SymmapT, vars::Vector{SymbolicT})
    for v in vars
        isconst(v) && continue
        symbolic_type(v) == NotSymbolic() && continue
        haskey(guesses, v) && continue
        def = getguess(v)
        if def !== nothing
            guesses[v] = SU.Const{VartypeT}(def)
            continue
        end
        Moshi.Match.@match v begin
            BSImpl.Term(; f, args) && if f === getindex end => begin
                haskey(guesses, args[1]) && continue
                def = Symbolics.getdefaultval(args[1], nothing)
                def === nothing && continue
                guesses[args[1]] = def
            end
            _ => nothing
        end
    end
    return guesses
end
function collect_var_to_name!(vars::Dict{Symbol, SymbolicT}, xs::Vector{SymbolicT})
    for x in xs
        x = Moshi.Match.@match x begin
            BSImpl.Const(;) => continue
            BSImpl.Term(; f, args) && if f === getindex end => args[1]
            _ => x
        end
        hasname(x) || continue
        vars[getname(x)] = x
    end
end
""""""
@noinline function throw_invalid_operator(opvar, eq, op::Type)
    if op === Differential
        optext = "derivative"
    end
    msg = "The $optext variable must be isolated to the left-hand " *
          "side of the equation like `$opvar ~ ...`. You may want to use `mtkcompile` or the DAE form.\nGot $eq."
    throw(InvalidSystemException(msg))
end
""""""
function _check_operator_variables(eq, op::T, expr = eq.rhs) where {T}
    iscall(expr) || return nothing
    if operation(expr) isa op
        throw_invalid_operator(expr, eq, op)
    end
    foreach(expr -> _check_operator_variables(eq, op, expr),
        SymbolicUtils.arguments(expr))
end
""""""
function check_operator_variables(eqs, op::T) where {T}
    ops = Set{SymbolicT}()
    tmp = Set{SymbolicT}()
    for eq in eqs
        _check_operator_variables(eq, op)
        SU.search_variables!(tmp, eq.lhs; is_atomic = OperatorIsAtomic{Differential}())
        if length(tmp) == 1
            x = only(tmp)
            if op === Differential
                is_tmp_fine = isdifferential(x)
            else
                is_tmp_fine = iscall(x) && !(operation(x) isa op)
            end
        else
            nd = count(x -> iscall(x) && !(operation(x) isa op), tmp)
            is_tmp_fine = iszero(nd)
        end
        is_tmp_fine ||
            error("The LHS cannot contain nondifferentiated variables. Please run `mtkcompile` or use the DAE form.\nGot $eq")
        for v in tmp
            v in ops &&
                error("The LHS operator must be unique. Please run `mtkcompile` or use the DAE form. $v appears in LHS more than once.")
            push!(ops, v)
        end
        empty!(tmp)
    end
end
function isoperator(expr::SymbolicT, ::Type{op}) where {op <: SU.Operator}
    Moshi.Match.@match expr begin
        BSImpl.Term(; f) => f isa op
        _ => false
    end
end
isoperator(::Type{op}) where {op <: SU.Operator} = Base.Fix2(isoperator, op)
isdifferential(expr) = isoperator(expr, Differential)
isdiffeq(eq) = isdifferential(eq.lhs) || isoperator(eq.lhs, Shift)
isvariable(x::Num)::Bool = isvariable(value(x))
function isvariable(x)
    x isa SymbolicT || return false
    hasmetadata(x, VariableSource) || iscall(x) && operation(x) === getindex && isvariable(arguments(x)[1])::Bool
end
function collect_operator_variables(sys::AbstractSystem, args...)
    collect_operator_variables(equations(sys), args...)
end
function collect_operator_variables(eq::Equation, args...)
    collect_operator_variables([eq], args...)
end
""""""
function collect_operator_variables(eqs::Vector{Equation}, ::Type{op}) where {op}
    vars = Set{SymbolicT}()
    diffvars = Set{SymbolicT}()
    for eq in eqs
        SU.search_variables!(vars, eq; is_atomic = OperatorIsAtomic{op}())
        for v in vars
            isoperator(v, op) || continue
            push!(diffvars, arguments(v)[1])
        end
        empty!(vars)
    end
    return diffvars
end
collect_differential_variables(sys) = collect_operator_variables(sys, Differential)
""""""
function collect_applied_operators(x, ::Type{op}) where {op}
    v = Set{SymbolicT}()
    SU.search_variables!(v, x; is_atomic = OnlyOperatorIsAtomic{op}())
    return v
end
""""""
function collect_scoped_vars!(unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, sys::AbstractSystem, iv::Union{SymbolicT, Nothing}; depth = 1, op = Differential)
    if has_eqs(sys)
        for eq in equations(sys)
            eqtype_supports_collect_vars(eq) || continue
            if eq isa Equation
                symtype(eq.lhs) <: Number || continue
            end
            collect_vars!(unknowns, parameters, eq, iv; depth, op)
        end
    end
    if has_jumps(sys)
        for eq in jumps(sys)
            eqtype_supports_collect_vars(eq) || continue
            collect_vars!(unknowns, parameters, eq, iv; depth, op)
        end
    end
    if has_constraints(sys)
        for eq in constraints(sys)
            eqtype_supports_collect_vars(eq) || continue
            collect_vars!(unknowns, parameters, eq, iv; depth, op)
        end
    end
end
""""""
function validate_operator(op, args, iv; context = nothing)
    error("`$validate_operator` is not implemented for operator `$op` in $context.")
end
function validate_operator(op::Differential, args, iv; context = nothing)
    isequal(op.x, iv) || throw(OperatorIndepvarMismatchError(op, iv, context))
    arg = unwrap(only(args))
    if !is_variable_floatingpoint(arg)
        throw(ContinuousOperatorDiscreteArgumentError(op, arg, context))
    end
end
struct ContinuousOperatorDiscreteArgumentError <: Exception
    op::Any
    arg::Any
    context::Any
end
function Base.showerror(io::IO, err::ContinuousOperatorDiscreteArgumentError)
end
struct OperatorIndepvarMismatchError <: Exception
    op::Any
    iv::Any
    context::Any
end
function Base.showerror(io::IO, err::OperatorIndepvarMismatchError)
end
struct OnlyOperatorIsAtomic{O} end
function (::OnlyOperatorIsAtomic{O})(ex::SymbolicT) where {O}
    Moshi.Match.@match ex begin
        BSImpl.Term(; f) && if f isa O end => true
        _ => false
    end
end
struct OperatorIsAtomic{O} end
function (::OperatorIsAtomic{O})(ex::SymbolicT) where {O}
    SU.default_is_atomic(ex) && Moshi.Match.@match ex begin
        BSImpl.Term(; f) && if f isa Operator end => f isa O
        _ => true
    end
end
""""""
function collect_vars!(unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, expr::SymbolicT, iv::Union{SymbolicT, Nothing}; depth = 0, op = Symbolics.Operator)
    Moshi.Match.@match expr begin
        BSImpl.Const(;) => return
        BSImpl.Sym(;) => return collect_var!(unknowns, parameters, expr, iv; depth)
        _ => nothing
    end
    vars = Set{SymbolicT}()
    SU.search_variables!(vars, expr; is_atomic = OperatorIsAtomic{op}())
    for var in vars
        while iscall(var) && operation(var) isa op
            validate_operator(operation(var), arguments(var), iv; context = expr)
            var = arguments(var)[1]
        end
        collect_var!(unknowns, parameters, var, iv; depth)
    end
    return nothing
end
function collect_vars!(unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, expr::AbstractArray, iv::Union{SymbolicT, Nothing}; depth = 0, op = Symbolics.Operator)
    for var in expr
        collect_vars!(unknowns, parameters, var, iv; depth, op)
    end
    return nothing
end
""""""
eqtype_supports_collect_vars(eq) = false
eqtype_supports_collect_vars(eq::Equation) = true
eqtype_supports_collect_vars(eq::Inequality) = true
eqtype_supports_collect_vars(eq::Pair) = true
function collect_vars!(unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, eq::Union{Equation, Inequality}, iv::Union{SymbolicT, Nothing};
        depth = 0, op = Symbolics.Operator)
    collect_vars!(unknowns, parameters, eq.lhs, iv; depth, op)
    collect_vars!(unknowns, parameters, eq.rhs, iv; depth, op)
    return nothing
end
function collect_vars!(
        unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, p::Pair, iv::Union{SymbolicT, Nothing}; depth = 0, op = Symbolics.Operator)
    collect_vars!(unknowns, parameters, p[1], iv; depth, op)
    collect_vars!(unknowns, parameters, p[2], iv; depth, op)
    return nothing
end
function collect_vars!(
        unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, expr, iv::Union{SymbolicT, Nothing}; depth = 0, op = Symbolics.Operator)
    return nothing
end
""""""
function collect_var!(unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, var::SymbolicT, iv::Union{SymbolicT, Nothing}; depth = 0)
    isequal(var, iv) && return nothing
    if Symbolics.iswrapped(var)
        error("""
        Internal Error. Please open an issue with an MWE.
        Encountered a wrapped value in `collect_var!`. This function should only ever \
        receive unwrapped symbolic variables. This is likely a bug in the code generating \
        an expression passed to `collect_vars!` or `collect_scoped_vars!`. A common cause \
        is using `substitute` with rules where the values are \
        wrapped symbolic variables.
        """)
    end
    check_scope_depth(getmetadata(var, SymScope, LocalScope())::AllScopes, depth) || return nothing
    var = setmetadata(var, SymScope, LocalScope())
    if iscalledparameter(var)
        callable = getcalledparameter(var)
        push!(parameters, callable)
        collect_vars!(unknowns, parameters, arguments(var), iv)
    elseif isparameter(var) || (iscall(var) && isparameter(operation(var)))
        push!(parameters, var)
    else
        push!(unknowns, var)
    end
    if hasdefault(var) && (def = getdefault(var)) !== missing
        collect_vars!(unknowns, parameters, def, iv)
    end
    return nothing
end
""""""
function check_scope_depth(scope, depth)
    if scope isa LocalScope
        return depth == 0
    elseif scope isa ParentScope
        return depth > 0 && check_scope_depth(scope.parent, depth - 1)::Bool
    elseif scope isa GlobalScope
        return depth == -1
    end
end
isarray(x) = x isa AbstractArray || x isa Symbolics.Arr
""""""
function empty_substitutions(sys)
    isempty(observed(sys))
end
""""""
function get_substitutions(sys)
    obs, _ = unhack_observed(observed(sys), equations(sys))
    Dict([eq.lhs => eq.rhs for eq in obs])
end
@noinline function throw_missingvars_in_sys(vars)
    throw(ArgumentError("$vars are either missing from the variable map or missing from the system's unknowns/parameters list."))
end
function _with_unit(f, x, t, args...)
    x = f(x, args...)
    if hasmetadata(x, VariableUnit) && (t isa SymbolicT && hasmetadata(t, VariableUnit))
        xu = getmetadata(x, VariableUnit)
        tu = getmetadata(t, VariableUnit)
        x = setmetadata(x, VariableUnit, xu / tu)
    end
    return x
end
diff2term_with_unit(x, t) = _with_unit(diff2term, x, t)
lower_varname_with_unit(var, iv, order) = _with_unit(lower_varname, var, iv, iv, order)
shift2term_with_unit(x, t) = _with_unit(shift2term, x, t)
lower_shift_varname_with_unit(var, iv) = _with_unit(lower_shift_varname, var, iv, iv)
""""""
function is_variable_floatingpoint(sym)
    sym = unwrap(sym)
    T = symtype(sym)
    is_floatingpoint_symtype(T)
end
""""""
function is_floatingpoint_symtype(T)
    return T === Real || T === Number || T === Complex || T <: AbstractFloat ||
           T <: AbstractArray && is_floatingpoint_symtype(eltype(T))
end
abstract type ObservedGraphCacheKey end
struct ObservedGraphCache
    graph::DiCMOBiGraph{false, Int, BipartiteGraph{Int, Nothing},
        Matching{Unassigned, Vector{Union{Unassigned, Int}}}}
    obsvar_to_idx::Dict{Any, Int}
end
function similar_variable(var::BasicSymbolic, name = :anon; use_gensym = true)
    if use_gensym
        name = gensym(name)
    end
    stype = symtype(var)
    sym = Symbolics.variable(name; T = stype)
    if size(var) !== ()
        sym = setmetadata(sym, Symbolics.ArrayShapeCtx, map(Base.OneTo, size(var)))
    end
    return sym
end
""""""
function _eq_unordered(a::AbstractArray, b::AbstractArray)
    a = vec(a)
    b = vec(b)
    length(a) === length(b) || return false
    n = length(a)
    idxs = Set(1:n)
    for x in a
        idx = findfirst(isequal(x), b)
        while idx !== nothing && !(idx in idxs)
            idx = findnext(isequal(x), b, idx + 1)
        end
        idx === nothing && return false
        delete!(idxs, idx)
    end
    return true
end
_eq_unordered(a, b) = isequal(a, b)
""""""
function flatten_equations(eqs::Vector{Equation})
    _eqs = Equation[]
    for eq in eqs
        if !SU.is_array_shape(SU.shape(eq.lhs))
            push!(_eqs, eq)
            continue
        end
        lhs = vec(collect(eq.lhs)::Array{SymbolicT})::Vector{SymbolicT}
        rhs = vec(collect(eq.rhs)::Array{SymbolicT})::Vector{SymbolicT}
        for (l, r) in zip(lhs, rhs)
            push!(_eqs, l ~ r)
        end
    end
    return _eqs
end
const JumpType = Union{VariableRateJump, ConstantRateJump, MassActionJump}
struct NotPossibleError <: Exception end
function Base.showerror(io::IO, ::NotPossibleError)
    print(io, """
    This should not be possible. Please open an issue in ModelingToolkit.jl with an MWE.
    """)
end
