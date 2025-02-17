
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
@kwdef struct ImperativeAffect
    f::Any
    obs::Vector
    obs_syms::Vector{Symbol}
    modified::Vector
    mod_syms::Vector{Symbol}
    ctx::Any
    skip_checks::Bool
end

function ImperativeAffect(f;
        observed::NamedTuple = NamedTuple{()}(()),
        modified::NamedTuple = NamedTuple{()}(()),
        ctx = nothing,
        skip_checks = false)
    ImperativeAffect(f,
        collect(values(observed)), collect(keys(observed)),
        collect(values(modified)), collect(keys(modified)),
        ctx, skip_checks)
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

function Base.show(io::IO, mfa::ImperativeAffect)
    obs_vals = join(map((ob, nm) -> "$ob => $nm", mfa.obs, mfa.obs_syms), ", ")
    mod_vals = join(map((md, nm) -> "$md => $nm", mfa.modified, mfa.mod_syms), ", ")
    affect = mfa.f
    print(io,
        "ImperativeAffect(observed: [$obs_vals], modified: [$mod_vals], affect:$affect)")
end
func(f::ImperativeAffect) = f.f
context(a::ImperativeAffect) = a.ctx
observed(a::ImperativeAffect) = a.obs
observed_syms(a::ImperativeAffect) = a.obs_syms
function discretes(a::ImperativeAffect)
    Iterators.filter(ModelingToolkit.isparameter,
        Iterators.flatten(Iterators.map(
            x -> symbolic_type(x) == NotSymbolic() && x isa AbstractArray ? x : [x],
            a.modified)))
end
modified(a::ImperativeAffect) = a.modified
modified_syms(a::ImperativeAffect) = a.mod_syms

function Base.:(==)(a1::ImperativeAffect, a2::ImperativeAffect)
    isequal(a1.f, a2.f) && isequal(a1.obs, a2.obs) && isequal(a1.modified, a2.modified) &&
        isequal(a1.obs_syms, a2.obs_syms) && isequal(a1.mod_syms, a2.mod_syms) &&
        isequal(a1.ctx, a2.ctx)
end

function Base.hash(a::ImperativeAffect, s::UInt)
    s = hash(a.f, s)
    s = hash(a.obs, s)
    s = hash(a.obs_syms, s)
    s = hash(a.modified, s)
    s = hash(a.mod_syms, s)
    hash(a.ctx, s)
end

namespace_affects(af::ImperativeAffect, s) = namespace_affect(af, s)
function namespace_affect(affect::ImperativeAffect, s)
    ImperativeAffect(func(affect),
        namespace_expr.(observed(affect), (s,)),
        observed_syms(affect),
        renamespace.((s,), modified(affect)),
        modified_syms(affect),
        context(affect),
        affect.skip_checks)
end

function compile_affect(affect::ImperativeAffect, cb, sys, dvs, ps; kwargs...)
    compile_user_affect(affect, cb, sys, dvs, ps; kwargs...)
end

function compile_user_affect(affect::ImperativeAffect, cb, sys, dvs, ps; kwargs...)
    #=
    Implementation sketch:
        generate observed function (oop), should save to a component array under obs_syms
        do the same stuff as the normal FA for pars_syms
        call the affect method
        unpack and apply the resulting values
    =#
    function check_dups(syms, exprs) # = (syms_dedup, exprs_dedup)
        seen = Set{Symbol}()
        syms_dedup = []
        exprs_dedup = []
        for (sym, exp) in Iterators.zip(syms, exprs)
            if !in(sym, seen)
                push!(syms_dedup, sym)
                push!(exprs_dedup, exp)
                push!(seen, sym)
            elseif !affect.skip_checks
                @warn "Expression $(expr) is aliased as $sym, which has already been used. The first definition will be used."
            end
        end
        return (syms_dedup, exprs_dedup)
    end

    obs_exprs = observed(affect)
    if !affect.skip_checks
        for oexpr in obs_exprs
            invalid_vars = invalid_variables(sys, oexpr)
            if length(invalid_vars) > 0
                error("Observed equation $(oexpr) in affect refers to missing variable(s) $(invalid_vars); the variables may not have been added (e.g. if a component is missing).")
            end
        end
    end
    obs_syms = observed_syms(affect)
    obs_syms, obs_exprs = check_dups(obs_syms, obs_exprs)

    mod_exprs = modified(affect)
    if !affect.skip_checks
        for mexpr in mod_exprs
            if !check_assignable(sys, mexpr)
                @warn ("Expression $mexpr cannot be assigned to; currently only unknowns and parameters may be updated by an affect.")
            end
            invalid_vars = unassignable_variables(sys, mexpr)
            if length(invalid_vars) > 0
                error("Modified equation $(mexpr) in affect refers to missing variable(s) $(invalid_vars); the variables may not have been added (e.g. if a component is missing) or they may have been reduced away.")
            end
        end
    end
    mod_syms = modified_syms(affect)
    mod_syms, mod_exprs = check_dups(mod_syms, mod_exprs)

    overlapping_syms = intersect(mod_syms, obs_syms)
    if length(overlapping_syms) > 0 && !affect.skip_checks
        @warn "The symbols $overlapping_syms are declared as both observed and modified; this is a code smell because it becomes easy to confuse them and assign/not assign a value."
    end

    # sanity checks done! now build the data and update function for observed values
    mkzero(sz) =
        if sz === ()
            0.0
        else
            zeros(sz)
        end
    obs_fun = build_explicit_observed_function(
        sys, Symbolics.scalarize.(obs_exprs);
        mkarray = (es, _) -> MakeTuple(es))
    obs_sym_tuple = (obs_syms...,)

    # okay so now to generate the stuff to assign it back into the system
    mod_pairs = mod_exprs .=> mod_syms
    mod_names = (mod_syms...,)
    mod_og_val_fun = build_explicit_observed_function(
        sys, Symbolics.scalarize.(first.(mod_pairs));
        mkarray = (es, _) -> MakeTuple(es))

    upd_funs = NamedTuple{mod_names}((setu.((sys,), first.(mod_pairs))...,))

    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        save_idxs = get(ic.callback_to_clocks, cb, Int[])
    else
        save_idxs = Int[]
    end

    let user_affect = func(affect), ctx = context(affect)
        function (integ)
            # update the to-be-mutated values; this ensures that if you do a no-op then nothing happens
            modvals = mod_og_val_fun(integ.u, integ.p, integ.t)
            upd_component_array = NamedTuple{mod_names}(modvals)

            # update the observed values
            obs_component_array = NamedTuple{obs_sym_tuple}(obs_fun(
                integ.u, integ.p, integ.t))

            # let the user do their thing
            upd_vals = user_affect(upd_component_array, obs_component_array, ctx, integ)

            # write the new values back to the integrator
            _generated_writeback(integ, upd_funs, upd_vals)

            for idx in save_idxs
                SciMLBase.save_discretes!(integ, idx)
            end
        end
    end
end

scalarize_affects(affects::ImperativeAffect) = affects

function vars!(vars, aff::ImperativeAffect; op = Differential)
    for var in Iterators.flatten((observed(aff), modified(aff)))
        if symbolic_type(var) == NotSymbolic()
            if var isa AbstractArray
                for v in var
                    v = unwrap(v)
                    vars!(vars, v)
                end
            end
        else
            var = unwrap(var)
            vars!(vars, var)
        end
    end
    return vars
end
