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
    unique([
        filter(isoutput, states(sys))
        filter(x -> x isa Term && isoutput(x), rhss) # observed can return equations with complicated expressions, we are only looking for single Terms
        filter(x -> x isa Term && isoutput(x), lhss)
    ])
end

"""
    bound_inputs(sys)

Return inputs that are bound within the system, i.e., internal inputs
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
bound_inputs(sys) = filter(x->is_bound(sys, x), inputs(sys))

"""
    unbound_inputs(sys)

Return inputs that are not bound within the system, i.e., external inputs
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
unbound_inputs(sys) = filter(x->!is_bound(sys, x), inputs(sys))

"""
    bound_outputs(sys)

Return outputs that are bound within the system, i.e., internal outputs
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
bound_outputs(sys) = filter(x->is_bound(sys, x), outputs(sys))

"""
    unbound_outputs(sys)

Return outputs that are not bound within the system, i.e., external outputs
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
unbound_outputs(sys) = filter(x->!is_bound(sys, x), outputs(sys))

"""
    is_bound(sys, u)

Determine whether or not input/output variable `u` is "bound" within the system, i.e., if it's to be considered internal to `sys`.
A variable/signal is considered bound if it appears in an equation together with variables from other subsystems.
The typical usecase for this function is to determine whether the input to an IO component is connected to another component,
or if it remains an external input that the user has to supply before simulating the system. 

See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
function is_bound(sys, u, stack=[])
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
    eqs = filter(eq->has_var(eq, u), eqs) # Only look at equations that contain u
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
    oeqs = filter(eq->has_var(eq, u), oeqs) # Only look at equations that contain u
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
    join(parts[1:end-1], '₊')
end

"""
    has_var(eq, x)

Determine whether or not an equation or expression contains variable `x`.
"""
function has_var(eq::Equation, x)
    has_var(eq.rhs, x) || has_var(eq.lhs, x)
end

has_var(ex, x) = x ∈ Set(get_variables(ex))

