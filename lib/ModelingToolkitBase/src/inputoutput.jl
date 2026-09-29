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
    default_codegen_inputs(sys)

The inputs to use by default for input-aware code generation.

For scheduled (compiled) systems this is `inputs(sys)`: the inputs declared to
`mtkcompile` are the contract the simplified system was built around. The
[`unbound_inputs`](@ref) heuristic cannot be used there — after flattening and
simplification an effective input appears in equations together with variables
from other namespaces and is therefore classified as bound, so `unbound_inputs`
is empty for compiled hierarchical systems.

For unscheduled systems this is [`unbound_inputs`](@ref), which inspects the
connection structure of the hierarchy to find external inputs.
"""
default_codegen_inputs(sys) = isscheduled(sys) ? inputs(sys) : unbound_inputs(sys)

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
    generate_control_function(sys, inputs = default_codegen_inputs(sys),
        disturbance_inputs = disturbances(sys); kwargs...) -> (; f, dvs, ps, io_sys)

Generate the dynamics of an input-output system as callable functions of its state,
inputs, parameters, and independent variable.

# Arguments

- `sys::AbstractSystem`: The system to generate dynamics for. An unscheduled system is
  compiled with `mtkcompile`; a scheduled system is used as given.
- `inputs`: Symbolic variables that form the generated input argument `u`. By default,
  declared inputs are used for scheduled systems and external inputs for unscheduled systems.
- `disturbance_inputs`: Unknown disturbance inputs. Their state and dynamics are retained,
  but their values are set to zero and are not function arguments.

# Keywords

- `known_disturbance_inputs = nothing`: Disturbance inputs supplied as a final generated
  argument `w`; they are removed from the parameter arguments.
- `implicit_dae::Bool = false`: Generate residual dynamics for an implicit DAE.
- `simplify::Bool = false`: Forwarded to `mtkcompile` when `sys` is unscheduled.
- `split::Bool = true`: Forwarded to `mtkcompile` to select split-system generation.
- `eval_expression::Bool = false`: Evaluate generated code in `eval_module` instead of
  returning a runtime-generated function.
- `eval_module::Module = @__MODULE__`: Module used when `eval_expression = true`.
- `disturbance_argument = false`: Deprecated compatibility option. Use
  `known_disturbance_inputs` instead.
- `kwargs...`: Forwarded to `Symbolics.CodegenFunctionOptions`.

# Returns

A named tuple with:

- `f`: A pair `(f_oop, f_iip)` of generated out-of-place and in-place dynamics wrappers.
  The basic call signatures are `f_oop(x, u, p..., t)` and
  `f_iip(dx, x, u, p..., t)`. With known disturbances, both have a final `w` argument.
- `dvs`: The selected state variables, ordered as the `x` argument of `f`.
- `ps`: The selected parameter variables, ordered as the parameter arguments of `f`.
- `io_sys`: The scheduled system used to generate `f`.

# Example

```julia
using ModelingToolkitBase
import ModelingToolkitBase: t_nounits as t, D_nounits as D

@variables x(t) u(t)
@parameters k
@named sys = System([D(x) ~ -k * (x + u)], t)

(; f, dvs, ps, io_sys) = generate_control_function(sys, [u]; simplify = true)
p = [2.0]
f[1]([1.0], [3.0], p, 0.0) # [-8.0]
```
"""
function generate_control_function(
        sys::AbstractSystem, inputs, disturbance_inputs, opts::GeneratedFunctionOptions;
        known_disturbance_inputs = nothing,
        disturbance_argument = false,
        implicit_dae = false,
        simplify = false,
        split = true
    )
    (; eval_expression, eval_module) = opts
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
        sys = mtkcompile(sys; inputs, disturbance_inputs = all_disturbances, split, simplify)
    end

    # Add all disturbances to inputs for the purposes of io processing
    if !isempty(all_disturbances)
        inputs = [inputs; all_disturbances]
    end
    inputs = vec(unwrap_vars(inputs))
    dvs = unknowns(sys)
    ps::Vector{SymbolicT} = parameters(sys; initial_parameters = true)
    inputs_set = Set{SymbolicT}(inputs)
    for v in inputs
        arr, isarr = split_indexed_var(v)
        isarr || continue
        push!(inputs_set, arr)
    end
    ps = setdiff(ps, inputs_set)

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
        sys, rhss, collect(Any, args), BuildFunctionWrapperOptions(;
            u_arg = 1 + Int(implicit_dae), p_start = 3 + implicit_dae,
            p_end = length(p) + 2 + implicit_dae,
            codegen_function_options = opts.codegen
        )
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
