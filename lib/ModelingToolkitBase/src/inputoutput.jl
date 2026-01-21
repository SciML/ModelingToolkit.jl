using Symbolics: get_variables
"""
    inputs(sys)

Return all variables that mare marked as inputs. See also [`unbound_inputs`](@ref)
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref)
"""
inputs(sys) = collect(get_inputs(sys))

"""
    outputs(sys)

Return all variables that mare marked as outputs. See also [`unbound_outputs`](@ref)
See also [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
function outputs(sys)
    return collect(get_outputs(sys))
end

"""
    bound_inputs(sys)

Return inputs that are bound within the system, i.e., internal inputs
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
bound_inputs(sys) = filter(x -> is_bound(sys, x), inputs(sys))

"""
    unbound_inputs(sys)

Return inputs that are not bound within the system, i.e., external inputs
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
unbound_inputs(sys) = filter(x -> !is_bound(sys, x), inputs(sys))

"""
    bound_outputs(sys)

Return outputs that are bound within the system, i.e., internal outputs
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
bound_outputs(sys) = filter(x -> is_bound(sys, x), outputs(sys))

"""
    unbound_outputs(sys)

Return outputs that are not bound within the system, i.e., external outputs
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
unbound_outputs(sys) = filter(x -> !is_bound(sys, x), outputs(sys))

function _is_atomic_inside_operator(ex::SymbolicT)
    return SU.default_is_atomic(ex) && Moshi.Match.@match ex begin
        BSImpl.Term(; f) && if f isa Operator end => false
        _ => true
    end
end

struct IsBoundValidator
    eqs_vars::Vector{Set{SymbolicT}}
    obs_vars::Vector{Set{SymbolicT}}
    stack::OrderedSet{SymbolicT}
end

function IsBoundValidator(sys::System)
    eqs_vars = Set{SymbolicT}[]
    for eq in equations(sys)
        vars = Set{SymbolicT}()
        SU.search_variables!(vars, eq.rhs; is_atomic = _is_atomic_inside_operator)
        SU.search_variables!(vars, eq.lhs; is_atomic = _is_atomic_inside_operator)
        push!(eqs_vars, vars)
    end
    obs_vars = Set{SymbolicT}[]
    for eq in observed(sys)
        vars = Set{SymbolicT}()
        SU.search_variables!(vars, eq.rhs; is_atomic = _is_atomic_inside_operator)
        SU.search_variables!(vars, eq.lhs; is_atomic = _is_atomic_inside_operator)
        push!(obs_vars, vars)
    end
    return IsBoundValidator(eqs_vars, obs_vars, OrderedSet{SymbolicT}())
end

function (ibv::IsBoundValidator)(u::SymbolicT)
    #=
    For observed quantities, we check if a variable is connected to something that is bound to something further out.
    In the following scenario
    julia> observed(syss)
        2-element Vector{Equation}:
        sys₊y(tv) ~ sys₊x(tv)
        y(tv) ~ sys₊x(tv)
    sys₊y(t) is bound to the outer y(t) through the variable sys₊x(t) and should thus return is_bound(sys₊y(t)) = true.
    When asking is_bound(sys₊y(t)), we know that we are looking through observed equations and can thus ask
    if var is bound, if it is, then sys₊y(t) is also bound. This can lead to an infinite recursion, so we maintain a stack of variables we have previously asked about to be able to break cycles
    =#
    u in ibv.stack && return false # Cycle detected
    for vars in ibv.eqs_vars
        u in vars || continue
        for var in vars
            var === u && continue
            same_or_inner_namespace(u, var) || return true
        end
    end
    for vars in ibv.obs_vars
        u in vars || continue
        for var in vars
            var === u && continue
            same_or_inner_namespace(u, var) || return true
            push!(ibv.stack, u)
            isbound = ibv(var)
            pop!(ibv.stack)
            # The variable we are comparing to can not come from an inner namespace,
            # binding only counts outwards
            isbound && !inner_namespace(u, var) && return true
        end
    end
    return false
end

"""
    is_bound(sys, u)

Determine whether input/output variable `u` is "bound" within the system, i.e., if it's to be considered internal to `sys`.
A variable/signal is considered bound if it appears in an equation together with variables from other subsystems.
The typical usecase for this function is to determine whether the input to an IO component is connected to another component,
or if it remains an external input that the user has to supply before simulating the system.

See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
function is_bound(sys, u)
    return IsBoundValidator(sys)(unwrap(u))
end

"""
    same_or_inner_namespace(u, var)

Determine whether `var` is in the same namespace as `u`, or a namespace internal to the namespace of `u`.
Example: `sys.u ~ sys.inner.u` will bind `sys.inner.u`, but `sys.u` remains an unbound, external signal. The namespaced signal `sys.inner.u` lives in a namespace internal to `sys`.
"""
function same_or_inner_namespace(u, var)
    nu = get_namespace(u)
    nv = get_namespace(var)
    return nu == nv ||           # namespaces are the same
        startswith(nv, nu) || # or nv starts with nu, i.e., nv is an inner namespace to nu
        occursin(NAMESPACE_SEPARATOR, string(getname(var))) &&
        !occursin(NAMESPACE_SEPARATOR, string(getname(u))) # or u is top level but var is internal
end

function inner_namespace(u, var)
    nu = get_namespace(u)
    nv = get_namespace(var)
    nu == nv && return false
    return startswith(nv, nu) || # or nv starts with nu, i.e., nv is an inner namespace to nu
        occursin(NAMESPACE_SEPARATOR, string(getname(var))) &&
        !occursin(NAMESPACE_SEPARATOR, string(getname(u))) # or u is top level but var is internal
end

"""
    get_namespace(x)

Return the namespace of a variable as a string. If the variable is not namespaced, the string is empty.
"""
function get_namespace(x)
    sname = string(getname(x))
    parts = split(sname, NAMESPACE_SEPARATOR)
    if length(parts) == 1
        return ""
    end
    return join(parts[1:(end - 1)], NAMESPACE_SEPARATOR)
end

"""
    has_var(eq, x)

Determine whether an equation or expression contains variable `x`.
"""
function has_var(eq::Equation, x)
    return has_var(eq.rhs, x) || has_var(eq.lhs, x)
end

has_var(ex, x) = x ∈ Set(get_variables(ex))

# Build control function

"""
    (f_oop, f_ip), x_sym, p_sym, io_sys = generate_control_function(
            sys::System,
            inputs                   = unbound_inputs(sys),
            disturbance_inputs       = disturbances(sys);
            known_disturbance_inputs = nothing,
            implicit_dae             = false,
            simplify                 = false,
            split                    = true,
        )

For a system `sys` with inputs (as determined by [`unbound_inputs`](@ref) or user specified), generate functions with additional input argument `u`

The returned functions are the out-of-place (`f_oop`) and in-place (`f_ip`) forms:
```
f_oop : (x,u,p,t)      -> rhs         # basic form
f_oop : (x,u,p,t,w)    -> rhs         # with known_disturbance_inputs
f_ip  : (xout,x,u,p,t) -> nothing     # basic form
f_ip  : (xout,x,u,p,t,w) -> nothing   # with known_disturbance_inputs
```

The return values also include the chosen state-realization (the remaining unknowns) `x_sym` and parameters, in the order they appear as arguments to `f`.

# Disturbance Handling

- `disturbance_inputs`: Unknown disturbance inputs. The generated dynamics will preserve any state and dynamics associated with these disturbances, but the disturbance inputs themselves will not be included as function arguments. This is useful for state observers that estimate unmeasured disturbances.

- `known_disturbance_inputs`: Known disturbance inputs. The generated dynamics will preserve state and dynamics, and the disturbance inputs will be added as an additional input argument `w` to the generated function: `(x,u,p,t,w)->rhs`.

# Example

```julia
using ModelingToolkitBase: generate_control_function, varmap_to_vars, defaults
f, x_sym, ps = generate_control_function(sys, expression=Val{false}, simplify=false)
p = varmap_to_vars(defaults(sys), ps)
x = varmap_to_vars(defaults(sys), x_sym)
t = 0
f[1](x, inputs, p, t)
```
"""
function generate_control_function(
        sys::AbstractSystem, inputs = unbound_inputs(sys),
        disturbance_inputs = disturbances(sys);
        known_disturbance_inputs = nothing,
        disturbance_argument = false,
        implicit_dae = false,
        simplify = false,
        eval_expression = false,
        eval_module = @__MODULE__,
        split = true,
        kwargs...
    )
    isempty(inputs) && @warn("No unbound inputs were found in system.")

    # Handle backward compatibility for disturbance_argument
    if disturbance_argument
        Base.depwarn(
            "The `disturbance_argument` keyword argument is deprecated. Use `known_disturbance_inputs` instead. " *
                "For `disturbance_argument=true`, pass `known_disturbance_inputs=disturbance_inputs, disturbance_inputs=nothing`. " *
                "For `disturbance_argument=false`, use `disturbance_inputs` as before.",
            :generate_control_function
        )
        if known_disturbance_inputs !== nothing
            error("Cannot specify both `disturbance_argument=true` and `known_disturbance_inputs`")
        end
        known_disturbance_inputs = disturbance_inputs
        disturbance_inputs = nothing
    end

    # Collect all disturbance inputs for mtkcompile
    all_disturbances = vcat(
        disturbance_inputs === nothing ? [] : disturbance_inputs,
        known_disturbance_inputs === nothing ? [] : known_disturbance_inputs
    )

    if !isscheduled(sys)
        sys = mtkcompile(sys; inputs, disturbance_inputs = all_disturbances, split)
    end

    # Add all disturbances to inputs for the purposes of io processing
    if !isempty(all_disturbances)
        inputs = [inputs; all_disturbances]
    end
    inputs = vec(unwrap_vars(inputs))
    dvs = unknowns(sys)
    ps::Vector{SymbolicT} = parameters(sys; initial_parameters = true)
    ps = setdiff(ps, inputs)

    # Remove unknown disturbances from inputs (we don't want them as actual inputs to the dynamics)
    if disturbance_inputs !== nothing
        inputs = setdiff(inputs, disturbance_inputs)
    end

    inputs = map(value, inputs)

    # Prepare disturbance arrays for substitution and function arguments
    unknown_disturbances = disturbance_inputs === nothing ? [] : unwrap.(disturbance_inputs)
    known_disturbances = known_disturbance_inputs === nothing ? [] : unwrap.(known_disturbance_inputs)

    eqs = [eq for eq in full_equations(sys)]

    # Set unknown disturbance inputs to zero (we just want to keep the disturbance state)
    if !isempty(unknown_disturbances)
        subs = Dict(unknown_disturbances .=> 0)
        eqs = [eq.lhs ~ substitute(eq.rhs, subs) for eq in eqs]
    end
    check_operator_variables(eqs, Differential)
    # substitute x(t) by just x
    rhss = implicit_dae ? [_iszero(eq.lhs) ? eq.rhs : eq.rhs - eq.lhs for eq in eqs] :
        [eq.rhs for eq in eqs]

    # TODO: add an optional check on the ordering of observed equations
    p = reorder_parameters(sys, ps)
    t = get_iv(sys)

    # Construct args with known disturbances if provided
    if !isempty(known_disturbances)
        args = (dvs, inputs, p..., t, known_disturbances)
    else
        args = (dvs, inputs, p..., t)
    end
    if implicit_dae
        ddvs = map(Differential(get_iv(sys)), dvs)
        args = (ddvs, args...)
    end
    f = build_function_wrapper(
        sys, rhss, args...; p_start = 3 + implicit_dae,
        p_end = length(p) + 2 + implicit_dae, kwargs...
    )
    f = eval_or_rgf.(f; eval_expression, eval_module)
    f = GeneratedFunctionWrapper{
        (
            3 + implicit_dae, length(args) - length(p) + 1, is_split(sys),
        ),
    }(f...)
    # Return parameters excluding both control inputs and all disturbances
    ps = setdiff(parameters(sys), inputs, all_disturbances)
    return (; f = (f, f), dvs, ps, io_sys = sys)
end
