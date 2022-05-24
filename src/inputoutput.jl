using Symbolics: get_variables
"""
    inputs(sys)

Return all variables that mare marked as inputs. See also [`unbound_inputs`](@ref)
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref)
"""
inputs(sys) = filter(isinput, states(sys))

"""
    outputs(sys)

Return all variables that mare marked as outputs. See also [`unbound_outputs`](@ref)
See also [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
function outputs(sys)
    o = observed(sys)
    rhss = [eq.rhs for eq in o]
    lhss = [eq.lhs for eq in o]
    unique([filter(isoutput, states(sys))
            filter(x -> x isa Term && isoutput(x), rhss) # observed can return equations with complicated expressions, we are only looking for single Terms
            filter(x -> x isa Term && isoutput(x), lhss)])
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

"""
    is_bound(sys, u)

Determine whether or not input/output variable `u` is "bound" within the system, i.e., if it's to be considered internal to `sys`.
A variable/signal is considered bound if it appears in an equation together with variables from other subsystems.
The typical usecase for this function is to determine whether the input to an IO component is connected to another component,
or if it remains an external input that the user has to supply before simulating the system. 

See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
function is_bound(sys, u, stack = [])
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
    u ∈ Set(stack) && return false # Cycle detected
    eqs = equations(sys)
    eqs = filter(eq -> has_var(eq, u), eqs) # Only look at equations that contain u
    # isout = isoutput(u)
    for eq in eqs
        vars = [get_variables(eq.rhs); get_variables(eq.lhs)]
        for var in vars
            var === u && continue
            if !same_or_inner_namespace(u, var)
                return true
            end
        end
    end
    # Look through observed equations as well
    oeqs = observed(sys)
    oeqs = filter(eq -> has_var(eq, u), oeqs) # Only look at equations that contain u
    for eq in oeqs
        vars = [get_variables(eq.rhs); get_variables(eq.lhs)]
        for var in vars
            var === u && continue
            if !same_or_inner_namespace(u, var)
                return true
            end
            if is_bound(sys, var, [stack; u]) && !inner_namespace(u, var) # The variable we are comparing to can not come from an inner namespace, binding only counts outwards
                return true
            end
        end
    end
    false
end

"""
    same_or_inner_namespace(u, var)

Determine whether or not `var` is in the same namespace as `u`, or a namespace internal to the namespace of `u`.
Example: `sys.u ~ sys.inner.u` will bind `sys.inner.u`, but `sys.u` remains an unbound, external signal. The namepsaced signal `sys.inner.u` lives in a namspace internal to `sys`.
"""
function same_or_inner_namespace(u, var)
    nu = get_namespace(u)
    nv = get_namespace(var)
    nu == nv ||           # namespaces are the same
        startswith(nv, nu) || # or nv starts with nu, i.e., nv is an inner namepsace to nu
        occursin('₊', var) && !occursin('₊', u) # or u is top level but var is internal
end

function inner_namespace(u, var)
    nu = get_namespace(u)
    nv = get_namespace(var)
    nu == nv && return false
    startswith(nv, nu) || # or nv starts with nu, i.e., nv is an inner namepsace to nu
        occursin('₊', var) && !occursin('₊', u) # or u is top level but var is internal
end

"""
    get_namespace(x)

Return the namespace of a variable as a string. If the variable is not namespaced, the string is empty.
"""
function get_namespace(x)
    sname = string(x)
    parts = split(sname, '₊')
    if length(parts) == 1
        return ""
    end
    join(parts[1:(end - 1)], '₊')
end

"""
    has_var(eq, x)

Determine whether or not an equation or expression contains variable `x`.
"""
function has_var(eq::Equation, x)
    has_var(eq.rhs, x) || has_var(eq.lhs, x)
end

has_var(ex, x) = x ∈ Set(get_variables(ex))

# Build control function

"""
    (f_oop, f_ip), dvs, p = generate_control_function(sys::AbstractODESystem, dvs = states(sys), ps = parameters(sys); implicit_dae = false, ddvs = if implicit_dae

For a system `sys` that has unbound inputs (as determined by [`unbound_inputs`](@ref)), generate a function with additional input argument `in`
```
f_oop : (u,in,p,t)      -> rhs
f_ip  : (uout,u,in,p,t) -> nothing
```
The return values also include the remaining states and parameters, in the order they appear as arguments to `f`.

# Example
```
using ModelingToolkit: generate_control_function, varmap_to_vars, defaults
f, dvs, ps = generate_control_function(sys, expression=Val{false}, simplify=true)
p = varmap_to_vars(defaults(sys), ps)
x = varmap_to_vars(defaults(sys), dvs)
t = 0
f[1](x, inputs, p, t)
```
"""
function generate_control_function(sys::AbstractODESystem;
                                   implicit_dae = false,
                                   has_difference = false,
                                   simplify = true,
                                   kwargs...)
    ctrls = unbound_inputs(sys)
    if isempty(ctrls)
        error("No unbound inputs were found in system.")
    end

    # One can either connect unbound inputs to new parameters and allow structural_simplify, but then the unbound inputs appear as states :( .
    # One can also just remove them from the states and parameters for the purposes of code generation, but then structural_simplify fails :(
    # To have the best of both worlds, all unbound inputs must be converted to `@parameters` in which case structural_simplify handles them correctly :)
    sys = toparam(sys, ctrls)

    if simplify
        sys = structural_simplify(sys)
    end

    dvs = states(sys)
    ps = parameters(sys)

    dvs = setdiff(dvs, ctrls)
    ps = setdiff(ps, ctrls)
    inputs = map(x -> time_varying_as_func(value(x), sys), ctrls)

    eqs = [eq for eq in equations(sys) if !isdifferenceeq(eq)]
    check_operator_variables(eqs, Differential)
    # substitute x(t) by just x
    rhss = implicit_dae ? [_iszero(eq.lhs) ? eq.rhs : eq.rhs - eq.lhs for eq in eqs] :
           [eq.rhs for eq in eqs]

    # TODO: add an optional check on the ordering of observed equations
    u = map(x -> time_varying_as_func(value(x), sys), dvs)
    p = map(x -> time_varying_as_func(value(x), sys), ps)
    t = get_iv(sys)

    # pre = has_difference ? (ex -> ex) : get_postprocess_fbody(sys)

    args = (u, inputs, p, t)
    if implicit_dae
        ddvs = map(Differential(get_iv(sys)), dvs)
        args = (ddvs, args...)
    end
    pre, sol_states = get_substitutions_and_solved_states(sys)
    f = build_function(rhss, args...; postprocess_fbody = pre, states = sol_states,
                       kwargs...)
    f, dvs, ps
end

"""
    toparam(sys, ctrls::AbstractVector)

Transform all instances of `@varibales` in `ctrls` appearing as states and in equations of `sys` with similarly named `@parameters`. This allows [`structural_simplify`](@ref)(sys) in the presence unbound inputs.
"""
function toparam(sys, ctrls::AbstractVector)
    eqs = equations(sys)
    subs = Dict(ctrls .=> toparam.(ctrls))
    eqs = map(eqs) do eq
        substitute(eq.lhs, subs) ~ substitute(eq.rhs, subs)
    end
    ODESystem(eqs, name = nameof(sys))
end
