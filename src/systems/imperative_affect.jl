
"""
    ImperativeAffect(f::Function; modified::NamedTuple, observed::NamedTuple, ctx)

`ImperativeAffect` is a helper for writing affect functions that will compute observed values and
ensure that modified values are correctly written back into the system. The affect function `f` needs to have
the signature 

```
    f(modified::NamedTuple, observed::NamedTuple, ctx, integrator)::NamedTuple
```

The function `f` will be called with `observed` and `modified` `NamedTuple`s that are derived from their respective `NamedTuple` definitions.
Each  declaration`NamedTuple` should map an expression to a symbol; for example if we pass `observed=(; x = a + b)` this will alias the result of executing `a+b` in the system as `x`
so the value of `a + b` will be accessible as `observed.x` in `f`. `modified` currently restricts symbolic expressions to only bare variables, so only tuples of the form
`(; x = y)` or `(; x)` (which aliases `x` as itself) are allowed.

The argument NamedTuples (for instance `(;x=y)`) will be populated with the declared values on function entry; if we require `(;x=y)` in `observed` and `y=2`, for example,
then the NamedTuple `(;x=2)` will be passed as `observed` to the affect function `f`. 

The NamedTuple returned from `f` includes the values to be written back to the system after `f` returns. For example, if we want to update the value of `x` to be the result of `x + y` we could write

    ImperativeAffect(observed=(; x_plus_y = x + y), modified=(; x)) do m, o
        @set! m.x = o.x_plus_y
    end

Where we use Setfield to copy the tuple `m` with a new value for `x`, then return the modified value of `m`. All values updated by the tuple must have names originally declared in
`modified`; a runtime error will be produced if a value is written that does not appear in `modified`. The user can dynamically decide not to write a value back by not including it
in the returned tuple, in which case the associated field will not be updated.
"""
struct ImperativeAffect
    f::Any
    observed::NamedTuple
    modified::NamedTuple
    ctx::Any
    skip_checks::Bool
end

function ImperativeAffect(f;
        observed::NamedTuple = NamedTuple{()}(()),
        modified::NamedTuple = NamedTuple{()}(()),
        ctx = nothing,
        skip_checks = false)
    ImperativeAffect(f, observed, modified, ctx, skip_checks)
end
function ImperativeAffect(f, modified::NamedTuple;
        observed::NamedTuple = NamedTuple{()}(()), ctx = nothing, skip_checks = false)
    ImperativeAffect(
        f, observed = observed, modified = modified, ctx = ctx, skip_checks = skip_checks)
end
function ImperativeAffect(
        f, modified::NamedTuple, observed::NamedTuple; ctx = nothing, skip_checks = false)
    ImperativeAffect(
        f, observed = observed, modified = modified, ctx = ctx, skip_checks = skip_checks)
end
function ImperativeAffect(
        f, modified::NamedTuple, observed::NamedTuple, ctx; skip_checks = false)
    ImperativeAffect(
        f, observed = observed, modified = modified, ctx = ctx, skip_checks = skip_checks)
end
function ImperativeAffect(; f, kwargs...)
    ImperativeAffect(f; kwargs...)
end

function Base.show(io::IO, mfa::ImperativeAffect)
    obs = mfa.observed
    mod = mfa.modified
    affect = mfa.f
    print(io,
        "ImperativeAffect(observed: [$(obs)], modified: [$(mod)], affect:$affect)")
end
func(f::ImperativeAffect) = f.f
context(a::ImperativeAffect) = a.ctx
function discretes(a::ImperativeAffect)
    Iterators.filter(ModelingToolkit.isparameter,
        Iterators.flatten(Iterators.map(
            x -> symbolic_type(x) == NotSymbolic() && x isa AbstractArray ? x : [x],
            a.modified)))
end

function Base.:(==)(a1::ImperativeAffect, a2::ImperativeAffect)
    isequal(a1.f, a2.f) && isequal(a1.observed, a2.observed) &&
        isequal(a1.modified, a2.modified) &&
        isequal(a1.ctx, a2.ctx)
end

function Base.hash(a::ImperativeAffect, s::UInt)
    s = hash(a.f, s)
    s = hash(a.observed, s)
    s = hash(a.modified, s)
    hash(a.ctx, s)
end

namespace_affects(af::ImperativeAffect, s) = namespace_affect(af, s)

function _namespace_nt(nt::NamedTuple, s::AbstractSystem)
    return NamedTuple{keys(nt)}(_namespace_nt(values(nt), s))
end

function _namespace_nt(nt::Union{AbstractArray, Tuple}, s::AbstractSystem)
    return map(nt) do v
        if symbolic_type(v) == NotSymbolic()
            _namespace_nt(v, s)
        else
            renamespace(s, v)
        end
    end
end

function namespace_affect(affect::ImperativeAffect, s)
    obs = _namespace_nt(affect.observed, s)
    mod = _namespace_nt(affect.modified, s)
    ImperativeAffect(affect.f, obs, mod, affect.ctx, affect.skip_checks)
end

function invalid_variables(sys, expr)
    filter(x -> !any(isequal(x), all_symbols(sys)), reduce(vcat, vars(expr); init = []))
end

function unassignable_variables(sys, expr)
    assignable_syms = reduce(
        vcat, Symbolics.scalarize.(vcat(
            unknowns(sys), parameters(sys; initial_parameters = true)));
        init = [])
    written = reduce(vcat, Symbolics.scalarize.(vars(expr)); init = [])
    return filter(
        x -> !any(isequal(x), assignable_syms), written)
end

function check_assignable(sys, sym)
    if symbolic_type(sym) == ScalarSymbolic()
        is_variable(sys, sym) || is_parameter(sys, sym)
    elseif symbolic_type(sym) == ArraySymbolic()
        is_variable(sys, sym) || is_parameter(sys, sym) ||
            all(x -> check_assignable(sys, x), collect(sym))
    elseif sym isa Union{AbstractArray, Tuple}
        all(x -> check_assignable(sys, x), sym)
    else
        false
    end
end

function _nt_check_valid(nt::NamedTuple, s::AbstractSystem, isobserved::Bool)
    _nt_check_valid(values(nt), s, isobserved)
end

function _nt_check_valid(
        nt::Union{Tuple, AbstractArray}, s::AbstractSystem, isobserved::Bool)
    for v in nt
        if symbolic_type(v) == NotSymbolic()
            _nt_check_valid(v, s, isobserved)
            continue
        end
        if !isobserved && !check_assignable(s, v)
            error("""
            Expression $v cannot be assigned to; currently only unknowns and parameters may \
            be updated by an affect.
            """)
        end
        invalid = invalid_variables(s, v)
        isempty(invalid) && continue
        name = isobserved ? "Observed" : "Modified"
        error("""
        $name expression $(v) in affect refers to missing variable(s) $(invalid); \
        the variables may not have been added (e.g. if a component is missing).
        """)
    end
end

function _nt_check_overlap(nta::NamedTuple, ntb::NamedTuple)
    common = intersect(keys(nta), keys(ntb))
    isempty(common) && return
    @warn """
    The symbols $common are declared as both observed and modified; this is a code smell \
    because it becomes easy to confuse them and assign/not assign a value.
    """
end

function compile_functional_affect(
        affect::ImperativeAffect, sys; reset_jumps = false, kwargs...)
    #=
    Implementation sketch:
        generate observed function (oop), should save to a component array under obs_syms
        do the same stuff as the normal FA for pars_syms
        call the affect method
        unpack and apply the resulting values
    =#

    if !affect.skip_checks
        _nt_check_valid(affect.observed, sys, true)
        _nt_check_valid(affect.modified, sys, false)
        _nt_check_overlap(affect.observed, affect.modified)
    end

    # sanity checks done! now build the data and update function for observed values
    let user_affect = func(affect), ctx = context(affect),
        obs_getter = isempty(affect.observed) ? Returns((;)) : getsym(sys, affect.observed),
        mod_getter = isempty(affect.modified) ? Returns((;)) : getsym(sys, affect.modified),
        mod_setter = isempty(affect.modified) ? Returns((;)) : setsym(sys, affect.modified),
        reset_jumps = reset_jumps

        @inline function (integ)
            mod = mod_getter(integ)
            obs = obs_getter(integ)

            # let the user do their thing
            upd_vals = user_affect(mod, obs, ctx, integ)
            mod_setter(integ, upd_vals)

            reset_jumps && reset_aggregated_jumps!(integ)
        end
    end
end

scalarize_affects(affects::ImperativeAffect) = affects

function _vars_nt!(vars, nt::NamedTuple, op)
    _vars_nt!(vars, values(nt), op)
end

function _vars_nt!(vars, nt::Union{AbstractArray, Tuple}, op)
    for v in nt
        if symbolic_type(v) == NotSymbolic()
            _vars_nt!(vars, v, op)
            continue
        end
        vars!(vars, v; op)
    end
end

function vars!(vars, aff::ImperativeAffect; op = Differential)
    _vars_nt!(vars, aff.observed, op)
    _vars_nt!(vars, aff.modified, op)
    return vars
end
