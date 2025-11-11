"""
    union_nothing(x::Union{T1, Nothing}, y::Union{T2, Nothing}) where {T1, T2}

Unite x and y gracefully when they could be nothing. If neither is nothing, x and y are united normally. If one is nothing, the other is returned unmodified. If both are nothing, nothing is returned.
"""
function union_nothing(x::Union{T1, Nothing}, y::Union{T2, Nothing}) where {T1, T2}
    isnothing(x) && return y # y can be nothing or something
    isnothing(y) && return x # x can be nothing or something
    return union(x, y) # both x and y are something and can be united normally
end

get_iv(D::Differential) = D.x

"""
    $(TYPEDSIGNATURES)

Turn `x(t)` into `x`
"""
function detime_dvs(op)
    if !iscall(op)
        op
    elseif issym(operation(op))
        SSym(nameof(operation(op)); type = Real, shape = SU.shape(op))
    else
        maketerm(typeof(op), operation(op), detime_dvs.(arguments(op)),
            metadata(op))
    end
end

"""
    $(TYPEDSIGNATURES)

Reverse `detime_dvs` for the given `dvs` using independent variable `iv`.
"""
function retime_dvs(op, dvs, iv)
    issym(op) && return SSym(nameof(op); type = FnType{Tuple{symtype(iv)}, Real}, shape = SU.ShapeVecT())(iv)
    iscall(op) ?
    maketerm(typeof(op), operation(op), retime_dvs.(arguments(op), (dvs,), (iv,)),
        metadata(op)) :
    op
end

function modified_unknowns!(munknowns, e::Equation, unknownlist = nothing)
    get_variables!(munknowns, e.lhs, unknownlist)
end

function todict(d)
    eltype(d) <: Pair || throw(ArgumentError("The variable-value mapping must be a Dict."))
    d isa Dict ? d : Dict(d)
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
    return string(expr)
end

# System validation enums
"""
    const CheckNone

Value that can be provided to the `check` keyword of `System` to disable checking of input.
"""
const CheckNone = 0
"""
    const CheckAll

Value that can be provided to the `check` keyword of `System` to enable all input
validation.
"""
const CheckAll = 1 << 0
"""
    const CheckComponents

Value that can be provided to the `check` keyword of `System` to only enable checking of
basic components of the system, such as equations, variables, etc.
"""
const CheckComponents = 1 << 1
"""
    const CheckUnits

Value that can be provided to the `check` keyword of `System` to enable checking of units.
"""
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

function check_lhs(eq::Equation, ::Type{Differential}, dvs::Set)
    v = unwrap(eq.lhs)
    _iszero(v) && return
    op = operation(v)
    op isa Differential && isone(op.order) && only(arguments(v)) in dvs && return
    error("$v is not a valid LHS. Please run mtkcompile before simulation.")
end
function check_lhs(eqs::Vector{Equation}, ::Type{Differential}, dvs::Set)
    for eq in eqs
        check_lhs(eq, Differential, dvs)
    end
end

"""
    collect_ivs(eqs, op = Differential)

Get all the independent variables with respect to which differentials (`op`) are taken.
"""
function collect_ivs(eqs, ::Type{op} = Differential) where {op}
    vars = Set{SymbolicT}()
    ivs = Set{SymbolicT}()
    for eq in eqs
        SU.search_variables!(vars, eq; is_atomic = OperatorIsAtomic{op}())
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

"""
    check_equations(eqs, iv)

Assert that equations are well-formed when building ODE, i.e., only containing a single independent variable.
"""
function check_equations(eqs::Vector{Equation}, iv::SymbolicT)
    icp = IndepvarCheckPredicate(iv)
    for eq in eqs
        SU.query(icp, eq.lhs)
        SU.query(icp, eq.rhs)
    end
end

"""
    $(TYPEDSIGNATURES)

Assert that the subsystems have the appropriate namespacing behavior.
"""
function check_subsystems(systems)
    idxs = findall(!does_namespacing, systems)
    isempty(idxs) || throw_bad_namespacing(systems, idxs)
end

@noinline function throw_bad_namespacing(systems, idxs)
    names = join("  " .* string.(nameof.(systems[idxs])), "\n")
    throw(ArgumentError("All subsystems must have namespacing enabled. The following subsystems do not perform namespacing:\n$(names)"))
end

"""
Get all the independent variables with respect to which differentials are taken.
"""
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

"""
    $TYPEDSIGNATURES

Check the validity of `bindings` given the list of parameters `ps`. This method assumes
that there are no discrete values in `ps`.
"""
function check_bindings(ps::Vector{SymbolicT}, bindings::SymmapT)
    atomic_ps = AtomicArraySet()
    for p in ps
        push_as_atomic_array!(atomic_ps, p)
    end
    check_bindings(atomic_ps, bindings)
end

function check_bindings_is_atomic(x::SymbolicT)
    SU.default_is_atomic(x) && Moshi.Match.@match x begin
        BSImpl.Term(; f) && if f isa Operator end => f isa Initial
        BSImpl.Term(; f) && if f === getindex end => false
        _ => true
    end
end

"""
    $TYPEDSIGNATURES

Check if `bindings` are valid, given a list of parameters `atomic_ps`. Assumes no values in
`atomic_ps` are discretes.
"""
function check_bindings(atomic_ps::AtomicArraySet{Dict{SymbolicT, Nothing}}, bindings::SymmapT)
    varsbuf = Set{SymbolicT}()
    for p in atomic_ps
        val = get(bindings, p, COMMON_NOTHING)
        val === COMMON_NOTHING && continue
        if val === COMMON_MISSING
            if !is_variable_floatingpoint(p)
                throw(ArgumentError("""
                `missing` bindings are only valid for solvable parameters! Non-floating \
                point parameters cannot be solved for, and thus do not accept a binding \
                of `missing`. Found invalid parameter $p of symtype $(symtype(p)).
                """))
            end
        end
        empty!(varsbuf)
        SU.search_variables!(varsbuf, val; is_atomic = check_bindings_is_atomic)
        setdiff!(varsbuf, atomic_ps)
        filter!(x -> getmetadata(x, SymScope, LocalScope()) isa LocalScope, varsbuf)
        filter!(!isinitial, varsbuf)
        isempty(varsbuf) && continue

        throw(ArgumentError("""
        Bindings for parameters can only be functions of other parameters. For parameter \
        $p, encountered binding $val which contains non-parameter symbolics $varsbuf. If \
        you intended $p to be a discrete variable, pass it as an unknown of the system.
        """))
    end
end

function check_no_parameter_equations_recurse(ex::SymbolicT)
    iscall(ex) && !check_bindings_is_atomic(ex)
end

"""
    $(TYPEDSIGNATURES)

Validate that all equations of the system involve the unknowns/observables.
"""
function check_no_parameter_equations(sys::AbstractSystem)
    if !isempty(get_systems(sys))
        throw(ArgumentError("Expected flattened system"))
    end
    varsbuf = Set{SymbolicT}()
    pareqs = Equation[]
    allowed_vars = as_atomic_array_set(unknowns(sys))
    foreach(Base.Fix1(push_as_atomic_array!, allowed_vars), observables(sys))
    foreach(Base.Fix1(push_as_atomic_array!, allowed_vars), get_all_discretes_fast(sys))
    for eq in equations(sys)
        empty!(varsbuf)
        Symbolics.get_variables!(varsbuf, eq, allowed_vars; is_atomic = check_bindings_is_atomic, recurse = check_no_parameter_equations_recurse)
        isempty(varsbuf) && push!(pareqs, eq)
    end

    if !isempty(pareqs)
        error("""
        The equations of a system must involve the unknowns/observables. The following \
        equations were found to have no unknowns/observables:
        $(join(string.(pareqs), "\n"))
        """)
    end
end

"""
    $TYPEDSIGNATURES

Verify that bound parameters have not been provided initial conditions. Requires the \
existence of an up-to-data `parameter_bindings_graph`.
"""
function check_no_bound_initial_conditions(sys::AbstractSystem)
    bound_ps = (get_parameter_bindings_graph(sys)::ParameterBindingsGraph).bound_ps
    ics = initial_conditions(sys)
    bound_ics = intersect(bound_ps, keys(ics))
    isempty(bound_ics) || throw(BoundInitialConditionsError(collect(bound_ics)))
end

struct BoundInitialConditionsError <: Exception
    bound_pars::Vector{SymbolicT}
end

function Base.showerror(io::IO, err::BoundInitialConditionsError)
    print(io, """
    Bound parameters cannot have initial conditions. The following bound parameters \
    were found to have initial conditions:
    $(join(string.(collect(err.bound_pars)), "\n"))
    """)
end

"""
    $(TYPEDSIGNATURES)

Check if the symbolic variable `v` has a default value.
"""
hasdefault(v) = hasmetadata(v, Symbolics.VariableDefaultValue)
"""
    $(TYPEDSIGNATURES)

Return the default value of symbolic variable `v`.
"""
getdefault(v) = value(Symbolics.getdefaultval(v))
"""
    $(TYPEDSIGNATURES)

Set the default value of symbolic variable `v` to `val`.
"""
function setdefault(v, val)
    val === nothing ? v : wrap(setdefaultval(unwrap(v), value(val)))
end

function process_variables!(var_to_name::Dict{Symbol, SymbolicT}, initial_conditions::SymmapT, bindings::SymmapT, guesses::SymmapT, vars::Vector{SymbolicT})
    collect_defaults!(initial_conditions, bindings, vars)
    collect_guesses!(guesses, vars)
    collect_var_to_name!(var_to_name, vars)
    return nothing
end

function process_variables!(var_to_name::Dict{Symbol, SymbolicT}, initial_conditions::SymmapT, bindings::SymmapT, vars::Vector{SymbolicT})
    collect_defaults!(initial_conditions, bindings, vars)
    collect_var_to_name!(var_to_name, vars)
    return nothing
end

function collect_defaults!(initial_conditions::SymmapT, bindings::SymmapT, v::SymbolicT)
    if hasname(v) && occursin(NAMESPACE_SEPARATOR, string(getname(v)))
        return
    end
    Moshi.Match.@match v begin
        BSImpl.Const(;) => return
        BSImpl.Term(; f, args) && if f === getindex end => begin
            collect_defaults!(initial_conditions, bindings, args[1])
        end
        _ => begin
            def = Symbolics.getdefaultval(v, nothing)
            def === nothing && return
            def = BSImpl.Const{VartypeT}(def)
            Moshi.Match.@match def begin
                # `get!` here is just shorthand for "if the key doesn't exist, add this
                # value".
                BSImpl.Const(;) => if def === COMMON_MISSING
                    get!(bindings, v, def)
                else
                    get!(initial_conditions, v, def)
                end
                _ => get!(bindings, v, def)
            end
        end
    end
end

function collect_defaults!(initial_conditions::SymmapT, bindings::SymmapT, vars::Vector{SymbolicT})
    for v in vars
        collect_defaults!(initial_conditions, bindings, v)
    end
end

function collect_guesses!(guesses::SymmapT, v::SymbolicT)
    Moshi.Match.@match v begin
        BSImpl.Const(;) => return
        BSImpl.Term(; f, args) && if f === getindex end => begin
            collect_guesses!(guesses, args[1])
        end
        _ => begin
            def = getguess(v)
            def === nothing && return
            get!(guesses, v, BSImpl.Const{VartypeT}(def))
        end
    end
end
function collect_guesses!(guesses::SymmapT, vars::Vector{SymbolicT})
    for v in vars
        collect_guesses!(guesses, v)
    end
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

"""
Throw error when difference/derivative operation occurs in the R.H.S.
"""
@noinline function throw_invalid_operator(opvar, eq, op::Type)
    if op === Differential
        optext = "derivative"
    end
    msg = "The $optext variable must be isolated to the left-hand " *
          "side of the equation like `$opvar ~ ...`. You may want to use `mtkcompile` or the DAE form.\nGot $eq."
    throw(InvalidSystemException(msg))
end

"""
Check if difference/derivative operation occurs in the R.H.S. of an equation
"""
function _check_operator_variables(eq, op::T, expr = eq.rhs) where {T}
    iscall(expr) || return nothing
    if operation(expr) isa op
        throw_invalid_operator(expr, eq, op)
    end
    foreach(expr -> _check_operator_variables(eq, op, expr),
        SymbolicUtils.arguments(expr))
end
"""
Check if all the LHS are unique
"""
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

isoperator(::Any, ::Type{T}) where {T} = false
isoperator(ex::Union{Num, Arr, CallAndWrap}, ::Type{op}) where {op} = isoperator(unwrap(ex), op)
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

"""
    collect_operator_variables(eqs::Vector{Equation}, ::Type{op}) where {op}

Return a `Set` containing all variables that have Operator `op` applied to them.
See also [`collect_differential_variables`](@ref).
"""
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

"""
    collect_applied_operators(x, op)

Return  a `Set` with all applied operators in `x`, example:

```
@independent_variables t
@variables u(t) y(t)
D = Differential(t)
eq = D(y) ~ u
ModelingToolkit.collect_applied_operators(eq, Differential) == Set([D(y)])
```

The difference compared to `collect_operator_variables` is that `collect_operator_variables` returns the variable without the operator applied.
"""
function collect_applied_operators(x, ::Type{op}) where {op}
    v = Set{SymbolicT}()
    SU.search_variables!(v, x; is_atomic = OnlyOperatorIsAtomic{op}())
    return v
end

"""
    $(TYPEDSIGNATURES)

Search through equations and parameter dependencies of `sys`, where sys is at a depth of
`depth` from the root system, looking for variables scoped to the root system. Also
recursively searches through all subsystems of `sys`, increasing the depth if it is not
`-1`. A depth of `-1` indicates searching for variables with `GlobalScope`.
"""
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

"""
    $(TYPEDSIGNATURES)

Check whether the usage of operator `op` is valid in a system with independent variable
`iv`. If the system is time-independent, `iv` should be `nothing`. Throw an appropriate
error if `op` is invalid. `args` are the arguments to `op`.

# Keyword arguments

- `context`: The place where the operator occurs in the system/expression, or any other
  relevant information. Useful for providing extra information in the error message.
"""
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
    print(io, """
    Operator $(err.op) expects continuous arguments, with a `symtype` such as `Number`,
    `Real`, `Complex` or a subtype of `AbstractFloat`. Found $(err.arg) with a symtype of
    $(symtype(err.arg))$(err.context === nothing ? "." : "in $(err.context).")
    """)
end

struct OperatorIndepvarMismatchError <: Exception
    op::Any
    iv::Any
    context::Any
end

function Base.showerror(io::IO, err::OperatorIndepvarMismatchError)
    print(io, """
    Encountered operator `$(err.op)` which has different independent variable than the \
    one used in the system `$(err.iv)`.
    """)
    if err.context !== nothing
        println(io)
        print(io, "Context:\n$(err.context)")
    end
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

"""
    $(TYPEDSIGNATURES)

Search through `expr` for all symbolic variables present in it. Populate `dvs` with
unknowns and `ps` with parameters present. `iv` should be the independent variable of the
system or `nothing` for time-independent systems. Expressions where the operator `isa op`
go through `validate_operator`.

`depth` is a keyword argument which indicates how many levels down `expr` is from the root
of the system hierarchy. This is used to resolve scoping operators. The scope of a variable
can be checked using `check_scope_depth`.

This function should return `nothing`.
"""
function collect_vars!(unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, expr::SymbolicT, iv::Union{SymbolicT, Nothing}; depth = 0, op = Symbolics.Operator)
    Moshi.Match.@match expr begin
        BSImpl.Const(;) => return
        BSImpl.Sym(;) => return collect_var!(unknowns, parameters, expr, iv; depth)
        _ => nothing
    end
    vars = OrderedSet{SymbolicT}()
    SU.search_variables!(vars, expr; is_atomic = OperatorIsAtomic{op}())
    for var in vars
        Moshi.Match.@match var begin
            BSImpl.Term(; f, args) && if f isa op end => begin
                validate_operator(f, args, iv; context = expr)
                isempty(args) && continue
                push!(vars, args[1])
            end
            _ => collect_var!(unknowns, parameters, var, iv; depth)
        end
    end
    return nothing
end

function collect_vars!(unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, expr::AbstractArray, iv::Union{SymbolicT, Nothing}; depth = 0, op = Symbolics.Operator)
    for var in expr
        collect_vars!(unknowns, parameters, var, iv; depth, op)
    end
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Indicate whether the given equation type (Equation, Pair, etc) supports `collect_vars!`.
Can be dispatched by higher-level libraries to indicate support.
"""
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
        unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, ex::Union{Num, Arr, CallAndWrap}, iv::Union{SymbolicT, Nothing}; depth = 0, op = Symbolics.Operator)
    collect_vars!(unknowns, parameters, unwrap(ex), iv; depth, op)
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

"""
    $(TYPEDSIGNATURES)

Identify whether `var` belongs to the current system using `depth` and scoping information.
Add `var` to `unknowns` or `parameters` appropriately, and search through any expressions
in known metadata of `var` using `collect_vars!`.
"""
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
    # Add also any parameters that appear only as defaults in the var
    if hasdefault(var) && (def = getdefault(var)) !== missing
        collect_vars!(unknowns, parameters, def, iv)
    end
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Check if the given `scope` is at a depth of `depth` from the root system. Only
returns `true` for `scope::GlobalScope` if `depth == -1`.
"""
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

"""
    $(TYPEDSIGNATURES)

Check if any variables were eliminated from the system as part of `mtkcompile`.
"""
function empty_substitutions(sys)
    isempty(observed(sys))
end

"""
    $(TYPEDSIGNATURES)

Get a dictionary mapping variables eliminated from the system during `mtkcompile` to the
expressions used to calculate them.
"""
function get_substitutions(sys)
    obs, _ = unhack_observed(observed(sys), equations(sys))
    Dict([eq.lhs => eq.rhs for eq in obs])
end

@noinline function throw_missingvars_in_sys(vars)
    throw(ArgumentError("$vars are either missing from the variable map or missing from the system's unknowns/parameters list."))
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

    # return early if there is nothing to do
    #Base.isconcretetype(T) && (!tofloat || T === float(T)) && return vs # TODO: disabled float(T) to restore missing errors in https://github.com/SciML/ModelingToolkit.jl/issues/2873
    Base.isconcretetype(T) && !tofloat && return vs

    sym_vs = filter(x -> SymbolicUtils.issym(x) || SymbolicUtils.iscall(x), vs)
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
            C2 = promote_type(C, E)
            @assert C2 == E||C2 == C "`promote_to_concrete` can't make type $E uniform with $C"
            C = C2
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

"""
    $(TYPEDSIGNATURES)

Check if `sym` represents a symbolic floating point number or array of such numbers.
"""
function is_variable_floatingpoint(sym)
    sym = unwrap(sym)
    T = symtype(sym)
    is_floatingpoint_symtype(T)
end

"""
    $(TYPEDSIGNATURES)

Check if `T` is an appropriate symtype for a symbolic variable representing a floating
point number or array of such numbers.
"""
function is_floatingpoint_symtype(T)
    return T === Real || T === Number || T === Complex || T <: AbstractFloat ||
           T <: AbstractArray && is_floatingpoint_symtype(eltype(T))
end

"""
    $(TYPEDSIGNATURES)

Check if `sym` represents a symbolic number or array of numbers.
"""
function is_variable_numeric(sym)
    sym = unwrap(sym)
    T = symtype(sym)
    is_numeric_symtype(T)
end

"""
    $(TYPEDSIGNATURES)

Check if `T` is an appropriate symtype for a symbolic variable representing a number or
array of numbers.
"""
function is_numeric_symtype(T::Type)
    return T <: Number || T <: AbstractArray && is_numeric_symtype(eltype(T))
end

"""
    $(TYPEDSIGNATURES)

Return the `DiCMOBiGraph` denoting the dependencies between observed equations `eqs`.
"""
function observed_dependency_graph(eqs::Vector{Equation})
    for eq in eqs
        if symbolic_type(eq.lhs) == NotSymbolic()
            error("All equations must be observed equations of the form `var ~ expr`. Got $eq")
        end
    end
    graph, assigns = observed2graph(eqs, getproperty.(eqs, (:lhs,)))
    matching = complete(Matching(Vector{Union{Unassigned, Int}}(assigns)))
    return DiCMOBiGraph{false}(graph, matching)
end

abstract type ObservedGraphCacheKey end

struct ObservedGraphCache
    graph::DiCMOBiGraph{false, Int, BipartiteGraph{Int, Nothing},
        Matching{Unassigned, Vector{Union{Unassigned, Int}}}}
    obsvar_to_idx::Dict{Any, Int}
end

"""
    $(TYPEDSIGNATURES)

Return the indexes of observed equations of `sys` used by expression `exprs`.

Keyword arguments:
- `involved_vars`: A collection of the variables involved in `exprs`. This is the set of
  variables which will be explored to find dependencies on observed equations. Typically,
  providing this keyword is not necessary and is only useful to avoid repeatedly calling
  `vars(exprs)`
- `obs`: the list of observed equations.
- `available_vars`: If `exprs` involves a variable `x[1]`, this function will look for
  observed equations whose LHS is `x[1]` OR `x`. Sometimes, the latter is not required
  since `x[1]` might already be present elsewhere in the generated code (e.g. an argument
  to the function) but other elements of `x` are part of the observed equations, thus
  requiring them to be obtained from the equation for `x`. Any variable present in
  `available_vars` will not be searched for in the observed equations.
"""
function observed_equations_used_by(sys::AbstractSystem, exprs;
        involved_vars = nothing, obs = observed(sys), available_vars = Set{SymbolicT}())
    if involved_vars === nothing
        involved_vars = Set{SymbolicT}()
        SU.search_variables!(involved_vars, exprs; is_atomic = OperatorIsAtomic{Union{Shift, Differential, Initial}}())
    elseif !(involved_vars isa Set{SymbolicT})
        involved_vars = Set{SymbolicT}(involved_vars)
    end
    if !(available_vars isa Set)
        available_vars = Set(available_vars)
    end
    if iscomplete(sys) && obs == observed(sys)
        cache = getmetadata(sys, MutableCacheKey, nothing)
        obs_graph_cache = get!(cache, ObservedGraphCacheKey) do
            obsvar_to_idx = Dict{Any, Int}([eq.lhs => i for (i, eq) in enumerate(obs)])
            graph = observed_dependency_graph(obs)
            return ObservedGraphCache(graph, obsvar_to_idx)
        end
        @unpack obsvar_to_idx, graph = obs_graph_cache
    else
        obsvar_to_idx = Dict([eq.lhs => i for (i, eq) in enumerate(obs)])
        graph = observed_dependency_graph(obs)
    end

    obsidxs = BitSet()
    for sym in involved_vars
        sym in available_vars && continue
        arrsym = iscall(sym) && operation(sym) === getindex ? arguments(sym)[1] : nothing
        idx = @something(get(obsvar_to_idx, sym, nothing),
            get(obsvar_to_idx, arrsym, nothing),
            Some(nothing))
        idx === nothing && continue
        idx in obsidxs && continue
        parents = dfs_parents(graph, idx)
        for i in eachindex(parents)
            parents[i] == 0 && continue
            push!(obsidxs, i)
        end
    end

    obsidxs = collect(obsidxs)
    sort!(obsidxs)
    return obsidxs
end

"""
    $(TYPEDSIGNATURES)

Given an expression `expr`, return a dictionary mapping subexpressions of `expr` that do
not involve variables in `vars` to anonymous symbolic variables. Also return the modified
`expr` with the substitutions indicated by the dictionary. If `expr` is a function
of only `vars`, then all of the returned subexpressions can be precomputed.

Note that this will only process subexpressions floating point value. Additionally,
array variables must be passed in both scalarized and non-scalarized forms in `vars`.
"""
function subexpressions_not_involving_vars(expr, vars)
    expr = unwrap(expr)
    vars = map(unwrap, vars)
    state = Dict()
    newexpr = subexpressions_not_involving_vars!(expr, vars, state)
    return state, newexpr
end

"""
    $(TYPEDSIGNATURES)

Mutating version of `subexpressions_not_involving_vars` which writes to `state`. Only
returns the modified `expr`.
"""
function subexpressions_not_involving_vars!(expr, vars, state::Dict{Any, Any})
    expr = unwrap(expr)
    if symbolic_type(expr) == NotSymbolic()
        if is_array_of_symbolics(expr)
            return map(expr) do el
                subexpressions_not_involving_vars!(el, vars, state)
            end
        end
        return expr
    end
    any(isequal(expr), vars) && return expr
    iscall(expr) || return expr
    symbolic_has_known_size(expr) || return expr
    haskey(state, expr) && return state[expr]
    op = operation(expr)
    args = arguments(expr)
    # if this is a `getindex` and the getindex-ed value is a `Sym`
    # or it is not a called parameter
    # OR
    # none of `vars` are involved in `expr`
    if op === getindex && (issym(args[1]) || !iscalledparameter(args[1])) ||
       (vs = SU.search_variables(expr); intersect!(vs, vars); isempty(vs))
        sym = gensym(:subexpr)
        var = similar_variable(expr, sym)
        state[expr] = var
        return var
    end

    if (op == (+) || op == (*)) && symbolic_type(expr) !== ArraySymbolic()
        indep_args = SymbolicT[]
        dep_args = SymbolicT[]
        for arg in args
            _vs = SU.search_variables(arg)
            intersect!(_vs, vars)
            if !isempty(_vs)
                push!(dep_args, subexpressions_not_involving_vars!(arg, vars, state))
            else
                push!(indep_args, arg)
            end
        end
        indep_term = reduce(op, indep_args; init = Int(op == (*)))
        indep_term = subexpressions_not_involving_vars!(indep_term, vars, state)
        dep_term = reduce(op, dep_args; init = Int(op == (*)))
        return op(indep_term, dep_term)
    end
    newargs = map(args) do arg
        subexpressions_not_involving_vars!(arg, vars, state)
    end
    return maketerm(typeof(expr), op, newargs, metadata(expr))
end

"""
    $(TYPEDSIGNATURES)

Create an anonymous symbolic variable of the same shape, size and symtype as `var`, with
name `gensym(name)`. Does not support unsized array symbolics.

If `use_gensym == false`, will not `gensym` the name.
"""
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

"""
    $(TYPEDSIGNATURES)

If `sym isa Symbol`, try and convert it to a symbolic by matching against symbolic
variables in `allsyms`. If `sym` is not a `Symbol` or no match was found, return
`sym` as-is.
"""
function symbol_to_symbolic(sys::AbstractSystem, sym; allsyms = all_symbols(sys))
    sym isa Symbol || return sym
    idx = findfirst(x -> (hasname(x) ? getname(x) : Symbol(x)) == sym, allsyms)
    idx === nothing && return sym
    sym = allsyms[idx]
    if iscall(sym) && operation(sym) == getindex
        sym = arguments(sym)[1]
    end
    return sym
end

"""
    $(TYPEDSIGNATURES)

Check if `var` is present in `varlist`. `iv` is the independent variable of the system,
and should be `nothing` if not applicable.
"""
function var_in_varlist(var, varlist::AbstractSet, iv)
    var = unwrap(var)
    # simple case
    return var in varlist ||
           # indexed array symbolic, unscalarized array present
           (iscall(var) && operation(var) === getindex && arguments(var)[1] in varlist) ||
           # unscalarized sized array symbolic, all scalarized elements present
           (symbolic_type(var) == ArraySymbolic() && symbolic_has_known_size(var) &&
            all(x -> x in varlist, collect(var))) ||
           # delayed variables
           (isdelay(var, iv) && var_in_varlist(operation(var)(iv), varlist, iv))
end

"""
    $(TYPEDSIGNATURES)

Check if `a` and `b` contain identical elements, regardless of order. This is not
equivalent to `issetequal` because the latter does not account for identical elements that
have different multiplicities in `a` and `b`.
"""
function _eq_unordered(a::AbstractArray, b::AbstractArray)
    # a and b may be multidimensional
    # e.g. comparing noiseeqs of SDEs
    a = vec(a)
    b = vec(b)
    length(a) === length(b) || return false
    n = length(a)
    idxs = Set(1:n)
    for x in a
        idx = findfirst(isequal(x), b)
        # loop since there might be multiple identical entries in a/b
        # and while we might have already matched the first there could
        # be a second that is equal to x
        while idx !== nothing && !(idx in idxs)
            idx = findnext(isequal(x), b, idx + 1)
        end
        idx === nothing && return false
        delete!(idxs, idx)
    end
    return true
end

_eq_unordered(a, b) = isequal(a, b)

"""
    $(TYPEDSIGNATURES)

Given a list of equations where some may be array equations, flatten the array equations
without scalarizing occurrences of array variables and return the new list of equations.
"""
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

"""
    $(TYPEDSIGNATURES)

Given a vector of variables in the system, return the corresponding `Differential` form of variable if possible.
Else returns the variable as-is.
"""
function underscore_to_D(v::AbstractVector, sys)
    maps = isscheduled(sys) ? get_schedule(sys).dummy_sub : Dict()
    inv_maps = Dict{valtype(maps), Vector{Base.keytype(maps)}}()

    for (k, v) in maps
        push!(get!(() -> valtype(inv_maps)[], inv_maps, v), k)
    end
    iv = get_iv(sys)
    map(x -> underscore_to_D(x, iv, inv_maps), v)
end

function underscore_to_D(v, iv, inv_map)
    if haskey(inv_map, v)
        only(get(inv_map, v, [v]))
    else
        v = ModelingToolkit.detime_dvs(v)
        s = split(string(getname(v)), 'ˍ')
        if length(s) > 1
            n, suffix = s
        else
            n, suffix = first(s), ""
        end
        repeats = length(suffix) ÷ length(string(iv))
        D = Differential(iv)
        v = SSym(Symbol(n); type = FnType{Tuple, Real, Nothing}, shape = SymbolicUtils.ShapeVecT())(iv)
        wrap_with_D(v, D, repeats)
    end
end

function wrap_with_D(n, D, repeats)
    if repeats <= 0
        return n
    else
        wrap_with_D(D(n), D, repeats - 1)
    end
end

const DEFAULT_STABLE_INDEX = SU.StableIndex(Int[])

"""
    $TYPEDSIGNATURES

Given a symbolic variable `x`, check whether it is an indexed array symbolic. If it is,
return the array and `true`. Otherwise, return `x, false`.
"""
function split_indexed_var(x::SymbolicT)
    Moshi.Match.@match x begin
        BSImpl.Term(; f, args) && if f === getindex end => (args[1], true)
        BSImpl.Term(; f, args) && if f isa Operator && length(args) == 1 end => begin
            arr, isarr = split_indexed_var(args[1])
            isarr || return x, false
            return f(arr), isarr
        end
        _ => return x, false
    end
end

"""
    $TYPEDSIGNATURES

Given a symbolic variable `x`, assume `split_indexed_var(x)[2]` is `true`. Return the
corresponding `SymbolicUtils.StableIndex`.
"""
function get_stable_index(x::SymbolicT)
    Moshi.Match.@match x begin
        BSImpl.Term(; f, args) && if f === getindex end => return SU.StableIndex{Int}(x)
        BSImpl.Term(; f, args) && if f isa Operator end => return get_stable_index(args[1])
        _ => throw(ArgumentError("Invalid variable $x for `get_stable_index`."))
    end
end

"""
    $TYPEDSIGNATURES

Merge `b` into `a`, but error if `a` already contains that key. Return the modified `a`.
"""
function no_override_merge!(a::AbstractDict, b::AbstractDict)
    for (k, v) in b
        if haskey(a, k)
            throw(ArgumentError("Cannot merge without overriding: common key $k."))
        end
        a[k] = v
    end
    return a
end

"""
    $TYPEDSIGNATURES

Merge `b` into `a`, modifying `a`. For all keys common to `a` and `b`,
prefer the value in `a`.
"""
function left_merge!(a::AbstractDict, b::AbstractDict)
    mergewith!(first ∘ tuple, a, b)
end
