"""
    struct CondRewriter

Callable struct used to transform symbolic conditions into conditions involving discrete
variables.
"""
struct CondRewriter
    """
    The independent variable which the discrete variables depend on.
    """
    iv::BasicSymbolic
    """
    A mapping from a discrete variables to a `NamedTuple` containing the condition
    determining whether the discrete variable needs to be evaluated and the symbolic
    expression the discrete variable represents. The expression is used as a rootfinding
    function, and zero-crossings trigger re-evaluation of the condition (if `dependency`
    is `true`). `expression < 0` is evaluated on an up-crossing and `expression <= 0` is
    evaluated on a down-crossing to get the updated value of the condition variable.
    """
    conditions::Dict{Any, @NamedTuple{dependency, expression}}
end

function CondRewriter(iv)
    return CondRewriter(iv, Dict())
end

"""
    $(TYPEDSIGNATURES)

Given a symbolic condition `expr` and the condition `dep` it depends on, update the
mapping in `cw` and generate a new discrete variable if necessary.
"""
function new_cond_sym(cw::CondRewriter, expr, dep)
    if !iscall(expr) || operation(expr) != Base.:(<) || !iszero(arguments(expr)[2])
        throw(ArgumentError("`expr` passed to `new_cond_sym` must be of the form `f(args...) < 0`. Got $expr."))
    end
    # check if the same expression exists in the mapping
    existing_var = findfirst(p -> isequal(p.expression, expr), cw.conditions)
    if existing_var !== nothing
        # cache hit
        (existing_dep, _) = cw.conditions[existing_var]
        # update the dependency condition
        cw.conditions[existing_var] = (dependency = (dep | existing_dep), expression = expr)
        return existing_var
    end
    # generate a new condition variable
    cvar = gensym("cond")
    st = symtype(expr)
    iv = cw.iv
    cv = unwrap(first(@parameters $(cvar)(iv)::st = true)) # TODO: real init 
    cw.conditions[cv] = (dependency = dep, expression = expr)
    return cv
end

"""
Utility function for boolean implication.
"""
implies(a, b) = !a & b

"""
    $(TYPEDSIGNATURES)

Recursively rewrite conditions into discrete variables. `expr` is the condition to rewrite,
`dep` is a boolean expression/value which determines when the `expr` is to be evaluated. For
example, if `expr = expr1 | expr2` and `dep = dep1`, then `expr` should only be evaluated if
`dep1` evaluates to `true`. Recursively, `expr1` should only be evaluated if `dep1` is `true`,
and `expr2` should only be evaluated if `dep & !expr1`.

Returns a 3-tuple of the substituted expression, a condition describing when `expr` evaluates
to `true`, and a condition describing when `expr` evaluates to `false`.

This expects that all expressions with discontinuities or with discontinuous derivatives have
been rewritten into the form of `ifelse(rootfunc(args...) < 0, left(args...), right(args...))`.
The transformation is performed via `discontinuities_to_ifelse` using `Symbolics.rootfunction`
and family.
"""
function (cw::CondRewriter)(expr, dep)
    # single variable, trivial case
    if issym(expr) || iscall(expr) && issym(operation(expr))
        return (expr, expr, !expr)
        # literal boolean or integer
    elseif expr isa Bool
        return (expr, expr, !expr)
    elseif expr isa Int
        return (expr, true, true)
        # other singleton symbolic variables
    elseif !iscall(expr)
        @warn "Automatic conversion of if statments to events requires use of a limited conditional grammar; see the documentation. Skipping due to $expr"
        return (expr, true, true) # error case => conservative assumption is that both true and false have to be evaluated
    elseif operation(expr) == Base.:(|) # OR of two conditions
        a, b = arguments(expr)
        (rw_conda, truea, falsea) = cw(a, dep)
        # only evaluate second if first is false
        (rw_condb, trueb, falseb) = cw(b, dep & falsea)
        return (rw_conda | rw_condb, truea | trueb, falsea & falseb)

    elseif operation(expr) == Base.:(&) # AND of two conditions
        a, b = arguments(expr)
        (rw_conda, truea, falsea) = cw(a, dep)
        # only evaluate second if first is true
        (rw_condb, trueb, falseb) = cw(b, dep & truea)
        return (rw_conda & rw_condb, truea & trueb, falsea | falseb)
    elseif operation(expr) == ifelse
        c, a, b = arguments(expr)
        (rw_cond, ctrue, cfalse) = cw(c, dep)
        # only evaluate if condition is true
        (rw_conda, truea, falsea) = cw(a, dep & ctrue)
        # only evaluate if condition is false
        (rw_condb, trueb, falseb) = cw(b, dep & cfalse)
        # expression is true if condition is true and THEN branch is true, or condition is false
        # and ELSE branch is true
        # similarly for expression being false
        return (ifelse(rw_cond, rw_conda, rw_condb),
            implies(ctrue, truea) | implies(cfalse, trueb),
            implies(ctrue, falsea) | implies(cfalse, falseb))
    elseif operation(expr) == Base.:(!) # NOT of expression
        (a,) = arguments(expr)
        (rw, ctrue, cfalse) = cw(a, dep)
        return (!rw, cfalse, ctrue)
    elseif operation(expr) == Base.:(<)
        if !isequal(arguments(expr)[2], 0)
            throw(ArgumentError("Expected comparison to be written as `f(args...) < 0`. Found $expr."))
        end

        # if the comparison does not include time-dependent variables,
        # don't create a callback for it

        # Calling `expression_is_time_dependent` is `O(d)` where `d` is the depth of the
        # expression tree. We only call this in this here to avoid turning this into
        # an `O(d^2)` time complexity recursion, which would happen if it were called
        # at the beginning of the function. Now, it only happens near the leaves of
        # the recursive tree.
        if !expression_is_time_dependent(expr, cw.iv)
            return (expr, expr, !expr)
        end
        cv = new_cond_sym(cw, expr, dep)
        return (cv, cv, !cv)
    elseif operation(expr) == (==)
        # we don't touch equality since it's a point discontinuity. It's basically always
        # false for continuous variables. In case it's an equality between discrete
        # quantities, we don't need to transform it.
        return (expr, expr, !expr)
    elseif !expression_is_time_dependent(expr, cw.iv)
        return (expr, expr, !expr)
    end
    error("""
    Unsupported expression form in decision variable computation $expr. If the expression
    involves a registered function, declare the discontinuity using
    `Symbolics.@register_discontinuity`. If this is not meant to be transformed via
    `IfLifting`, wrap the parent expression in `ModelingToolkit.no_if_lift`.
    """)
end

"""
    $(TYPEDSIGNATURES)

Acts as the identity function, and prevents transformation of conditional expressions inside it. Useful
if specific `ifelse` or other functions with discontinuous derivatives shouldn't be transformed into
callbacks.
"""
no_if_lift(s) = s
@register_symbolic no_if_lift(s)

"""
    $(TYPEDEF)

A utility struct to search through an expression specifically for `ifelse` terms, and find
all variables used in the condition of such terms. The variables are stored in a field of
the struct.
"""
struct VarsUsedInCondition
    """
    Stores variables used in conditions of `ifelse` statements in the expression.
    """
    vars::Set{Any}
end

VarsUsedInCondition() = VarsUsedInCondition(Set())

function (v::VarsUsedInCondition)(expr)
    expr = Symbolics.unwrap(expr)
    if symbolic_type(expr) == NotSymbolic()
        is_array_of_symbolics(expr) || return
        foreach(v, expr)
        return
    end
    iscall(expr) || return
    op = operation(expr)

    # do not search inside no_if_lift to avoid discovering
    # redundant variables
    op == no_if_lift && return

    args = arguments(expr)
    if op == ifelse
        cond, branch_a, branch_b = arguments(expr)
        vars!(v.vars, cond)
        v(branch_a)
        v(branch_b)
    end
    foreach(v, args)
    return
end

"""
    $(TYPEDSIGNATURES)

Check if `expr` depends on the independent variable `iv`. Return `true` if `iv` is present
in the expression, `Differential(iv)` is in the expression, or a dependent variable such
as `@variables x(iv)` is in the expression.
"""
function expression_is_time_dependent(expr, iv)
    any(vars(expr)) do sym
        sym = unwrap(sym)
        isequal(sym, iv) && return true
        iscall(sym) || return false
        op = operation(sym)
        args = arguments(sym)
        op isa Differential && op == Differential(iv) ||
            issym(op) && length(args) == 1 && expression_is_time_dependent(args[1], iv)
    end
end

"""
    $(TYPEDSIGNATURES)

Given an expression `expr` which is to be evaluated if `dep` evaluates to `true`, transform
the conditions of all all `ifelse` statements in `expr` into functions of new discrete
variables. `cw` is used to store the information relevant to these newly introduced variables.
"""
function rewrite_ifs(cw::CondRewriter, expr, dep)
    expr = unwrap(expr)
    if symbolic_type(expr) == NotSymbolic()
        # non-symbolic expression might still be an array of symbolic expressions
        is_array_of_symbolics(expr) || return expr
        return map(ex -> rewrite_ifs(cw, ex, dep), expr)
    end

    iscall(expr) || return expr
    op = operation(expr)
    args = arguments(expr)
    # do not search into `no_if_lift`
    op == no_if_lift && return expr

    # transform `ifelse`
    if op == ifelse
        cond, iftrue, iffalse = args

        (rw_cond, deptrue, depfalse) = cw(cond, dep)
        rw_iftrue = rewrite_ifs(cw, iftrue, deptrue)
        rw_iffalse = rewrite_ifs(cw, iffalse, depfalse)
        return maketerm(
            typeof(expr), ifelse, [unwrap(rw_cond), rw_iftrue, rw_iffalse], metadata(expr))
    end

    # recurse into the rest of the cases
    args = map(ex -> rewrite_ifs(cw, ex, dep), args)
    return maketerm(typeof(expr), op, args, metadata(expr))
end

"""
    $(TYPEDSIGNATURES)

Return a modified `expr` where functions with known discontinuities or discontinuous
derivatives are transformed into `ifelse` statements. Utilizes the discontinuity API
in Symbolics. See [`Symbolics.rootfunction`](@ref),
[`Symbolics.left_continuous_function`](@ref), [`Symbolics.right_continuous_function`](@ref).

`iv` is the independent variable of the system. Only subexpressions of `expr` which
depend on `iv` are transformed.
"""
function discontinuities_to_ifelse(expr, iv)
    expr = unwrap(expr)
    if symbolic_type(expr) == NotSymbolic()
        # non-symbolic expression might still be an array of symbolic expressions
        is_array_of_symbolics(expr) || return expr
        return map(ex -> discontinuities_to_ifelse(ex, iv), expr)
    end

    iscall(expr) || return expr
    op = operation(expr)
    args = arguments(expr)
    # do not search into `no_if_lift`
    op == no_if_lift && return expr

    # Case I: the operation is symbolic.
    # We don't actually care if this is a callable parameter or not.
    # If it is, we want to search inside and perform if-lifting there.
    # If it isn't, either it's `x(t)` in which case this recursion is
    # effectively a no-op OR it's `x(f(t))` for DDEs and we want to
    # perform if-lifting inside.
    #
    # Case II: the operation is not symbolic.
    # We anyway want to recursively apply the transformation.
    #
    # Thus, we can do this here regardless of the subsequent checks
    args = map(ex -> discontinuities_to_ifelse(ex, iv), args)

    # if the operation is a known discontinuity
    if hasmethod(Symbolics.rootfunction, Tuple{typeof(op)})
        rootfn = Symbolics.rootfunction(op)
        leftfn = Symbolics.left_continuous_function(op)
        rightfn = Symbolics.right_continuous_function(op)
        rootexpr = rootfn(args...) < 0
        leftexpr = leftfn(args...)
        rightexpr = rightfn(args...)
        return maketerm(
            typeof(expr), ifelse, [rootexpr, leftexpr, rightexpr], metadata(expr))
    end

    return maketerm(typeof(expr), op, args, metadata(expr))
end

"""
    $(TYPEDSIGNATURES)

Generate the symbolic condition for discrete variable `sym`, which represents the condition
of an `ifelse` statement created through [`IfLifting`](@ref). This condition is used to
trigger a callback which updates the value of the condition appropriately.
"""
function generate_condition(cw::CondRewriter, sym)
    (dep, expr) = cw.conditions[sym]

    # expr is `f(args...) < 0`, `f(args...)` is the zero-crossing expression
    zero_crossing = arguments(expr)[1]

    # if we're meant to evaluate the condition, evaluate it. Otherwise, return `NaN`.
    # the solvers don't treat the transition from a number to NaN or back as a zero-crossing,
    # so it can be used to effectively disable the affect when the condition is not meant to
    # be evaluated.
    return ifelse(dep, zero_crossing, NaN) ~ 0
end

"""
    $(TYPEDSIGNATURES)

Generate the upcrossing and downcrossing affect functions for discrete variable `sym` involved
in `ifelse` statements that are lifted to callbacks using [`IfLifting`](@ref). `syms` is a
condition variable introduced by `cw`, and is thus a key in `cw.conditions`. `new_cond_vars`
is the list of all such new condition variables, corresponding to the order of vertices in
`new_cond_vars_graph`. `new_cond_vars_graph` is a directed graph where edges denote the
condition variables involved in the dependency expression of the source vertex.
"""
function generate_affects(cw::CondRewriter, sym, new_cond_vars, new_cond_vars_graph)
    sym_idx = findfirst(isequal(sym), new_cond_vars)
    if sym_idx === nothing
        throw(ArgumentError("Expected variable $sym to be a condition variable in $new_cond_vars."))
    end
    # use reverse direction of edges because instead of finding the variables it depends
    # on, we want the variables that depend on it
    parents = bfs_parents(new_cond_vars_graph, sym_idx; dir = :in)
    cond_vars_to_update = [new_cond_vars[i]
                           for i in eachindex(parents) if !iszero(parents[i])]
    update_syms = Symbol.(cond_vars_to_update)
    modified = NamedTuple{(update_syms...,)}(cond_vars_to_update)

    upcrossing_update_exprs = [arguments(last(cw.conditions[sym]))[1] < 0
                               for sym in cond_vars_to_update]
    upcrossing = ImperativeAffect(
        modified, observed = NamedTuple{(update_syms...,)}(upcrossing_update_exprs),
        skip_checks = true) do x, o, c, i
        return o
    end
    downcrossing_update_exprs = [arguments(last(cw.conditions[sym]))[1] <= 0
                                 for sym in cond_vars_to_update]
    downcrossing = ImperativeAffect(
        modified, observed = NamedTuple{(update_syms...,)}(downcrossing_update_exprs),
        skip_checks = true) do x, o, c, i
        return o
    end

    return upcrossing, downcrossing
end

const CONDITION_SIMPLIFIER = Rewriters.Fixpoint(Rewriters.Postwalk(Rewriters.Chain([
                                                                                    # simple boolean laws
                                                                                    (@rule (!!(~x)) => (~x))
                                                                                    (@rule ((~x) & true) => (~x))
                                                                                    (@rule ((~x) & false) => false)
                                                                                    (@rule ((~x) | true) => true)
                                                                                    (@rule ((~x) | false) => (~x))
                                                                                    (@rule ((~x) & !(~x)) => false)
                                                                                    (@rule ((~x) | !(~x)) => true)
                                                                                    # reversed order of the above, because it matters and `@acrule` refuses to do its job
                                                                                    (@rule (true & (~x)) => (~x))
                                                                                    (@rule (false & (~x)) => false)
                                                                                    (@rule (true | (~x)) => true)
                                                                                    (@rule (false | (~x)) => (~x))
                                                                                    (@rule (!(~x) & (~x)) => false)
                                                                                    (@rule (!(~x) | (~x)) => true)
                                                                                    # idempotent
                                                                                    (@rule ((~x) & (~x)) => (~x))
                                                                                    (@rule ((~x) | (~x)) => (~x))
                                                                                    # ifelse with determined branches
                                                                                    (@rule ifelse((~x), true, false) => (~x))
                                                                                    (@rule ifelse((~x), false, true) => !(~x))
                                                                                    # ifelse with identical branches
                                                                                    (@rule ifelse((~x), (~y), (~y)) => (~y))
                                                                                    (@rule ifelse((~x), (~y), !(~y)) => ((~x) &
                                                                                                                         (~y)))
                                                                                    (@rule ifelse((~x), !(~y), (~y)) => ((~x) &
                                                                                                                         !(~y)))
                                                                                    # ifelse with determined condition
                                                                                    (@rule ifelse(true, (~x), (~y)) => (~x))
                                                                                    (@rule ifelse(false, (~x), (~y)) => (~y))])))

"""
If lifting converts (nested) if statements into a series of continous events + a logically equivalent if statement + parameters.

Lifting proceeds through the following process:
* rewrite comparisons to be of the form eqn [op] 0; subtract the RHS from the LHS 
* replace comparisons with generated parameters; for each comparison eqn [op] 0, generate an event (dependent on op) that sets the parameter
"""
function IfLifting(sys::ODESystem)
    cw = CondRewriter(get_iv(sys))

    eqs = copy(equations(sys))
    obs = copy(observed(sys))

    # get variables used by `eqs`
    syms = vars(eqs)
    # get observed equations used by `eqs`
    obs_idxs = observed_equations_used_by(sys, eqs; involved_vars = syms)
    # and the variables used in those equations
    for i in obs_idxs
        vars!(syms, obs[i])
    end

    # get all integral variables used in conditions
    # this is used when performing the transformation on observed equations
    # since they are transformed differently depending on whether they are
    # discrete variables involved in a condition or not
    condition_vars = Set()
    # searcher struct
    # we can use the same one since it avoids iterating over duplicates
    vars_in_condition! = VarsUsedInCondition()
    for i in eachindex(eqs)
        eq = eqs[i]
        vars_in_condition!(eq.rhs)
        # also transform the equation
        eqs[i] = eq.lhs ~ rewrite_ifs(cw, discontinuities_to_ifelse(eq.rhs, cw.iv), true)
    end
    # also search through relevant observed equations
    for i in obs_idxs
        vars_in_condition!(obs[i].rhs)
    end
    # add to `condition_vars` after filtering out differential, parameter, independent and
    # non-integral variables
    for v in vars_in_condition!.vars
        v = unwrap(v)
        stype = symtype(v)
        if isdifferential(v) || isparameter(v) || isequal(v, get_iv(sys))
            continue
        end
        stype <: Union{Integer, AbstractArray{Integer}} && push!(condition_vars, v)
    end
    # transform observed equations
    for i in obs_idxs
        obs[i] = if obs[i].lhs in condition_vars
            obs[i].lhs ~ first(cw(discontinuities_to_ifelse(obs[i].rhs, cw.iv), true))
        else
            obs[i].lhs ~ rewrite_ifs(cw, discontinuities_to_ifelse(obs[i].rhs, cw.iv), true)
        end
    end

    # `rewrite_ifs` and calling `cw` generate a lot of redundant code, simplify it
    eqs = map(eqs) do eq
        eq.lhs ~ CONDITION_SIMPLIFIER(eq.rhs)
    end
    obs = map(obs) do eq
        eq.lhs ~ CONDITION_SIMPLIFIER(eq.rhs)
    end
    # also simplify dependencies
    for (k, v) in cw.conditions
        cw.conditions[k] = map(CONDITION_SIMPLIFIER ∘ unwrap, v)
    end

    # get directed graph where nodes are the new condition variables and edges from each
    # node denote the condition variables used in it's dependency expression

    # so we have an ordering for the vertices
    new_cond_vars = collect(keys(cw.conditions))
    # "observed" equations
    new_cond_dep_eqs = [v ~ cw.conditions[v] for v in new_cond_vars]
    # construct the graph as a `DiCMOBiGraph`
    new_cond_vars_graph = observed_dependency_graph(new_cond_dep_eqs)

    new_callbacks = continuous_events(sys)
    new_defaults = defaults(sys)
    new_ps = Vector{SymbolicParam}(parameters(sys))

    for var in new_cond_vars
        condition = generate_condition(cw, var)
        up_affect, down_affect = generate_affects(
            cw, var, new_cond_vars, new_cond_vars_graph)
        cb = SymbolicContinuousCallback([condition], up_affect; affect_neg = down_affect,
            initialize = up_affect, rootfind = SciMLBase.RightRootFind)

        push!(new_callbacks, cb)
        new_defaults[var] = getdefault(var)
        push!(new_ps, var)
    end

    @set! sys.defaults = new_defaults
    @set! sys.eqs = eqs
    # do not need to topsort because we didn't modify the order
    @set! sys.observed = obs
    @set! sys.continuous_events = new_callbacks
    @set! sys.ps = new_ps
    return sys
end
