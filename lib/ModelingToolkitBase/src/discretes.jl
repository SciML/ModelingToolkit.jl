function todiscrete_validate(s::SymbolicT)
    if !iscall(s)
        error(
            """
            `@discretes` cannot create time-independent variables. Encountered $s. Use \
            `@parameters` for this purpose.
            """
        )
    end
    return toparam(s)
end
function todiscrete_validate(s::Union{Num, Symbolics.Arr, Symbolics.CallAndWrap})
    return typeof(s)(todiscrete_validate(unwrap(s)))
end

"""
$(SIGNATURES)

Define one or more discrete variables, for use in events of continuous systems. Every
symbolic declared with this macro must be a dependent variable; declaring a
time-independent one is an error, since `@parameters` already covers that case.

Discrete variables are parameters of the system that events are allowed to write to, so
they are passed in the parameter list of [`System`](@ref).

# Examples

```julia
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@variables x(t)
@discretes c(t)
event = [x ~ 1.0] => [c ~ c + 1]
@named sys = System([D(x) ~ -c * x], t, [x], [c]; continuous_events = [event])
```

See also [`@independent_variables`](@ref),
[`@variables`](https://docs.sciml.ai/Symbolics/stable/manual/variables/#Symbolics.@variables)
and [`@constants`](@ref).
"""
macro discretes(xs...)
    return Symbolics.parse_vars(
        :discretes,
        Real,
        xs,
        todiscrete_validate
    )
end
