using Symbolics: get_variables
"""
    inputs(sys)

Return all variables that mare marked as inputs. See also [`unbound_inputs`](@ref)
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref)
"""
inputs(sys) = [filter(isinput, unknowns(sys)); filter(isinput, parameters(sys))]

"""
    outputs(sys)

Return all variables that mare marked as outputs. See also [`unbound_outputs`](@ref)
See also [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
function outputs(sys)
    o = observed(sys)
    rhss = [eq.rhs for eq in o]
    lhss = [eq.lhs for eq in o]
    unique([filter(isoutput, unknowns(sys))
            filter(isoutput, parameters(sys))
            filter(x -> iscall(x) && isoutput(x), rhss) # observed can return equations with complicated expressions, we are only looking for single Terms
            filter(x -> iscall(x) && isoutput(x), lhss)])
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

Determine whether input/output variable `u` is "bound" within the system, i.e., if it's to be considered internal to `sys`.
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

Determine whether `var` is in the same namespace as `u`, or a namespace internal to the namespace of `u`.
Example: `sys.u ~ sys.inner.u` will bind `sys.inner.u`, but `sys.u` remains an unbound, external signal. The namespaced signal `sys.inner.u` lives in a namespace internal to `sys`.
"""
function same_or_inner_namespace(u, var)
    nu = get_namespace(u)
    nv = get_namespace(var)
    nu == nv ||           # namespaces are the same
        startswith(nv, nu) || # or nv starts with nu, i.e., nv is an inner namespace to nu
        occursin(NAMESPACE_SEPARATOR, string(getname(var))) &&
            !occursin(NAMESPACE_SEPARATOR, string(getname(u))) # or u is top level but var is internal
end

function inner_namespace(u, var)
    nu = get_namespace(u)
    nv = get_namespace(var)
    nu == nv && return false
    startswith(nv, nu) || # or nv starts with nu, i.e., nv is an inner namespace to nu
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
    join(parts[1:(end - 1)], NAMESPACE_SEPARATOR)
end

"""
    has_var(eq, x)

Determine whether an equation or expression contains variable `x`.
"""
function has_var(eq::Equation, x)
    has_var(eq.rhs, x) || has_var(eq.lhs, x)
end

has_var(ex, x) = x ∈ Set(get_variables(ex))

# Build control function

"""
    (f_oop, f_ip), x_sym, p_sym, io_sys = generate_control_function(
            sys::AbstractODESystem,
            inputs             = unbound_inputs(sys),
            disturbance_inputs = nothing;
            implicit_dae       = false,
            simplify           = false,
        )

For a system `sys` with inputs (as determined by [`unbound_inputs`](@ref) or user specified), generate a function with additional input argument `in`

```
f_oop : (x,u,p,t)      -> rhs
f_ip  : (xout,x,u,p,t) -> nothing
```

The return values also include the chosen state-realization (the remaining unknowns) `x_sym` and parameters, in the order they appear as arguments to `f`.

If `disturbance_inputs` is an array of variables, the generated dynamics function will preserve any state and dynamics associated with disturbance inputs, but the disturbance inputs themselves will (by default) not be included as inputs to the generated function. The use case for this is to generate dynamics for state observers that estimate the influence of unmeasured disturbances, and thus require unknown variables for the disturbance model, but without disturbance inputs since the disturbances are not available for measurement. To add an input argument corresponding to the disturbance inputs, either include the disturbance inputs among the control inputs, or set `disturbance_argument=true`, in which case an additional input argument `w` is added to the generated function `(x,u,p,t,w)->rhs`.

!!! note "Un-simplified system"
    This function expects `sys` to be un-simplified, i.e., `structural_simplify` or `@mtkbuild` should not be called on the system before passing it into this function. `generate_control_function` calls a special version of `structural_simplify` internally.

# Example

```
using ModelingToolkit: generate_control_function, varmap_to_vars, defaults
f, x_sym, ps = generate_control_function(sys, expression=Val{false}, simplify=false)
p = varmap_to_vars(defaults(sys), ps)
x = varmap_to_vars(defaults(sys), x_sym)
t = 0
f[1](x, inputs, p, t)
```
"""
function generate_control_function(sys::AbstractODESystem, inputs = unbound_inputs(sys),
        disturbance_inputs = disturbances(sys);
        disturbance_argument = false,
        implicit_dae = false,
        simplify = false,
        eval_expression = false,
        eval_module = @__MODULE__,
        kwargs...)
    isempty(inputs) && @warn("No unbound inputs were found in system.")

    if disturbance_inputs !== nothing
        # add to inputs for the purposes of io processing
        inputs = [inputs; disturbance_inputs]
    end

    sys, _ = io_preprocessing(sys, inputs, []; simplify, kwargs...)

    dvs = unknowns(sys)
    ps = parameters(sys)
    ps = setdiff(ps, inputs)
    if disturbance_inputs !== nothing
        # remove from inputs since we do not want them as actual inputs to the dynamics
        inputs = setdiff(inputs, disturbance_inputs)
        # ps = [ps; disturbance_inputs]
    end
    inputs = map(x -> time_varying_as_func(value(x), sys), inputs)
    disturbance_inputs = unwrap.(disturbance_inputs)

    eqs = [eq for eq in full_equations(sys)]
    eqs = map(subs_constants, eqs)
    if disturbance_inputs !== nothing && !disturbance_argument
        # Set all disturbance *inputs* to zero (we just want to keep the disturbance state)
        subs = Dict(disturbance_inputs .=> 0)
        eqs = [eq.lhs ~ substitute(eq.rhs, subs) for eq in eqs]
    end
    check_operator_variables(eqs, Differential)
    # substitute x(t) by just x
    rhss = implicit_dae ? [_iszero(eq.lhs) ? eq.rhs : eq.rhs - eq.lhs for eq in eqs] :
           [eq.rhs for eq in eqs]

    # TODO: add an optional check on the ordering of observed equations
    u = map(x -> time_varying_as_func(value(x), sys), dvs)
    p = map(x -> time_varying_as_func(value(x), sys), ps)
    p = reorder_parameters(sys, p)
    t = get_iv(sys)

    # pre = has_difference ? (ex -> ex) : get_postprocess_fbody(sys)
    if disturbance_argument
        args = (u, inputs, p..., t, disturbance_inputs)
    else
        args = (u, inputs, p..., t)
    end
    if implicit_dae
        ddvs = map(Differential(get_iv(sys)), dvs)
        args = (ddvs, args...)
    end
    f = build_function_wrapper(sys, rhss, args...; p_start = 3 + implicit_dae,
        p_end = length(p) + 2 + implicit_dae)
    f = eval_or_rgf.(f; eval_expression, eval_module)
    (; f, dvs, ps, io_sys = sys)
end

function inputs_to_parameters!(state::TransformationState, io)
    check_bound = io === nothing
    @unpack structure, fullvars, sys = state
    @unpack var_to_diff, graph, solvable_graph = structure
    @assert solvable_graph === nothing

    inputs = BitSet()
    var_reidx = zeros(Int, length(fullvars))
    ninputs = 0
    nvar = 0
    new_parameters = []
    input_to_parameters = Dict()
    new_fullvars = []
    for (i, v) in enumerate(fullvars)
        if isinput(v) && !(check_bound && is_bound(sys, v))
            if var_to_diff[i] !== nothing
                error("Input $(fullvars[i]) is differentiated!")
            end
            push!(inputs, i)
            ninputs += 1
            var_reidx[i] = -1
            p = toparam(v)
            push!(new_parameters, p)
            input_to_parameters[v] = p
        else
            nvar += 1
            var_reidx[i] = nvar
            push!(new_fullvars, v)
        end
    end
    ninputs == 0 && return (state, 1:0)

    nvars = ndsts(graph) - ninputs
    new_graph = BipartiteGraph(nsrcs(graph), nvars, Val(false))

    for ie in 1:nsrcs(graph)
        for iv in 𝑠neighbors(graph, ie)
            iv = var_reidx[iv]
            iv > 0 || continue
            add_edge!(new_graph, ie, iv)
        end
    end

    new_var_to_diff = DiffGraph(nvars, true)
    for (i, v) in enumerate(var_to_diff)
        new_i = var_reidx[i]
        (new_i < 1 || v === nothing) && continue
        new_v = var_reidx[v]
        @assert new_v > 0
        new_var_to_diff[new_i] = new_v
    end
    @set! structure.var_to_diff = complete(new_var_to_diff)
    @set! structure.graph = complete(new_graph)

    @set! sys.eqs = isempty(input_to_parameters) ? equations(sys) :
                    fast_substitute(equations(sys), input_to_parameters)
    @set! sys.unknowns = setdiff(unknowns(sys), keys(input_to_parameters))
    ps = parameters(sys)

    if io !== nothing
        inputs, = io
        # Change order of new parameters to correspond to user-provided order in argument `inputs`
        d = Dict{Any, Int}()
        for (i, inp) in enumerate(new_parameters)
            d[inp] = i
        end
        permutation = [d[i] for i in inputs]
        new_parameters = new_parameters[permutation]
    end

    @set! sys.ps = [ps; new_parameters]

    @set! state.sys = sys
    @set! state.fullvars = new_fullvars
    @set! state.structure = structure
    base_params = length(ps)
    return state, (base_params + 1):(base_params + length(new_parameters)) # (1:length(new_parameters)) .+ base_params
end

"""
    DisturbanceModel{M}

The structure represents a model of a disturbance, along with the input variable that is affected by the disturbance. See [`add_input_disturbance`](@ref) for additional details and an example.

# Fields:

  - `input`: The variable affected by the disturbance.
  - `model::M`: A model of the disturbance. This is typically an `ODESystem`, but type that implements [`ModelingToolkit.get_disturbance_system`](@ref)`(dist::DisturbanceModel) -> ::ODESystem` is supported.
"""
struct DisturbanceModel{M}
    input::Any
    model::M
    name::Symbol
end
DisturbanceModel(input, model; name) = DisturbanceModel(input, model, name)

# Point of overloading for libraries, e.g., to be able to support disturbance models from ControlSystemsBase
function get_disturbance_system(dist::DisturbanceModel{<:ODESystem})
    dist.model
end

"""
    (f_oop, f_ip), augmented_sys, dvs, p = add_input_disturbance(sys, dist::DisturbanceModel, inputs = nothing)

Add a model of an unmeasured disturbance to `sys`. The disturbance model is an instance of [`DisturbanceModel`](@ref).

The generated dynamics functions `(f_oop, f_ip)` will preserve any state and dynamics associated with disturbance inputs, but the disturbance inputs themselves will not be included as inputs to the generated function. The use case for this is to generate dynamics for state observers that estimate the influence of unmeasured disturbances, and thus require state variables for the disturbance model, but without disturbance inputs since the disturbances are not available for measurement.

`dvs` will be the states of the simplified augmented system, consisting of the states of `sys` as well as the states of the disturbance model.

For MIMO systems, all inputs to the system has to be specified in the argument `inputs`

# Example

The example below builds a double-mass model and adds an integrating disturbance to the input

```julia
using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
t = ModelingToolkitStandardLibrary.Blocks.t

# Parameters
m1 = 1
m2 = 1
k = 1000 # Spring stiffness
c = 10   # Damping coefficient

@named inertia1 = Inertia(; J = m1)
@named inertia2 = Inertia(; J = m2)
@named spring = Spring(; c = k)
@named damper = Damper(; d = c)
@named torque = Torque(; use_support = false)

eqs = [connect(torque.flange, inertia1.flange_a)
       connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
       connect(inertia2.flange_a, spring.flange_b, damper.flange_b)]
model = ODESystem(eqs, t; systems = [torque, inertia1, inertia2, spring, damper],
                  name = :model)
model = complete(model)
model_outputs = [model.inertia1.w, model.inertia2.w, model.inertia1.phi, model.inertia2.phi]

# Disturbance model
@named dmodel = Blocks.StateSpace([0.0], [1.0], [1.0], [0.0]) # An integrating disturbance
@named dist = ModelingToolkit.DisturbanceModel(model.torque.tau.u, dmodel)
(f_oop, f_ip), augmented_sys, dvs, p = ModelingToolkit.add_input_disturbance(model, dist)
```

`f_oop` will have an extra state corresponding to the integrator in the disturbance model. This state will not be affected by any input, but will affect the dynamics from where it enters, in this case it will affect additively from `model.torque.tau.u`.
"""
function add_input_disturbance(sys, dist::DisturbanceModel, inputs = nothing; kwargs...)
    t = get_iv(sys)
    @variables d(t)=0 [disturbance = true]
    @variables u(t)=0 [input = true] # New system input
    dsys = get_disturbance_system(dist)

    if inputs === nothing
        all_inputs = [u]
    else
        i = findfirst(isequal(dist.input), inputs)
        if i === nothing
            throw(ArgumentError("Input $(dist.input) indicated in the disturbance model was not found among inputs specified to add_input_disturbance"))
        end
        all_inputs = convert(Vector{Any}, copy(inputs))
        all_inputs[i] = u # The input where the disturbance acts is no longer an input, the new input is u
    end

    eqs = [dsys.input.u[1] ~ d
           dist.input ~ u + dsys.output.u[1]]
    augmented_sys = ODESystem(eqs, t, systems = [dsys], name = gensym(:outer))
    augmented_sys = extend(augmented_sys, sys)

    (f_oop, f_ip), dvs, p, io_sys = generate_control_function(augmented_sys, all_inputs,
        [d]; kwargs...)
    (f_oop, f_ip), augmented_sys, dvs, p, io_sys
end
