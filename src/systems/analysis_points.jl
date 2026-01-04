"""
    $(TYPEDSIGNATURES)

Given a list of analysis points, break the connection for each and set the output to zero.
"""
function handle_loop_openings(sys::AbstractSystem, aps)
    for ap in canonicalize_ap(sys, aps)
        sys, (d_v,) = apply_transformation(Break(ap, true, true, true), sys)
        guesses = copy(get_guesses(sys))
        guesses[d_v] = if symbolic_type(d_v) == ArraySymbolic()
            fill(NaN, size(d_v))
        else
            NaN
        end
        @set! sys.guesses = guesses
    end
    return sys
end

const DOC_LOOP_OPENINGS = """
- `loop_openings`: A list of analysis points whose connections should be removed and
  the outputs set to the input as a part of the linear analysis.
"""

const DOC_SYS_MODIFIER = """
- `system_modifier`: A function taking the transformed system and applying any
  additional transformations, returning the modified system. The modified system
  is passed to `linearization_function`. 
"""

"""
    $(TYPEDSIGNATURES)

Utility function for linear analyses that apply a transformation `transform`, which
returns the added variables `(du, u)`, to each of the analysis points in `aps` and then
calls `linearization_function` with all the `du`s as inputs and `u`s as outputs. Returns
the linearization function and modified, simplified system.

# Keyword arguments

$DOC_LOOP_OPENINGS
$DOC_SYS_MODIFIER

All other keyword arguments are forwarded to `linearization_function`.
"""
function get_linear_analysis_function(
        sys::AbstractSystem, transform, aps; system_modifier = identity, loop_openings = [], kwargs...
    )
    dus = []
    us = []
    sys = handle_loop_openings(sys, loop_openings)
    aps = canonicalize_ap(sys, aps)
    for ap in aps
        sys, (du, u) = apply_transformation(transform(ap), sys)
        push!(dus, du)
        push!(us, u)
    end
    return linearization_function(system_modifier(sys), dus, us; kwargs...)
end
"""
    $(TYPEDSIGNATURES)

Return the sensitivity function for the analysis point(s) `aps`, and the modified system
simplified with the appropriate inputs and outputs.

# Keyword Arguments

$DOC_LOOP_OPENINGS
$DOC_SYS_MODIFIER

All other keyword arguments are forwarded to `linearization_function`.
"""
function get_sensitivity_function(sys::AbstractSystem, aps; kwargs...)
    return get_linear_analysis_function(sys, SensitivityTransform, aps; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Return the complementary sensitivity function for the analysis point(s) `aps`, and the
modified system simplified with the appropriate inputs and outputs.

# Keyword Arguments

$DOC_LOOP_OPENINGS
$DOC_SYS_MODIFIER

All other keyword arguments are forwarded to `linearization_function`.
"""
function get_comp_sensitivity_function(sys::AbstractSystem, aps; kwargs...)
    return get_linear_analysis_function(sys, ComplementarySensitivityTransform, aps; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Return the loop-transfer function for the analysis point(s) `aps`, and the modified
system simplified with the appropriate inputs and outputs.

# Keyword Arguments

$DOC_LOOP_OPENINGS
$DOC_SYS_MODIFIER

All other keyword arguments are forwarded to `linearization_function`.
"""
function get_looptransfer_function(sys::AbstractSystem, aps; kwargs...)
    return get_linear_analysis_function(sys, LoopTransferTransform, aps; kwargs...)
end

for f in [:get_sensitivity, :get_comp_sensitivity, :get_looptransfer]
    utility_fun = Symbol(f, :_function)
    @eval function $f(
            sys, ap, args...; loop_openings = [], system_modifier = identity,
            allow_input_derivatives = true, kwargs...
        )
        lin_fun,
            ssys = $(utility_fun)(
            sys, ap, args...; loop_openings, system_modifier, kwargs...
        )
        mats, extras = ModelingToolkit.linearize(ssys, lin_fun; allow_input_derivatives)
        return mats, ssys, extras
    end
end

"""
    sys, input_vars, output_vars = $(TYPEDSIGNATURES)

Apply analysis-point transformations to prepare a system for linearization.

Returns
- `sys`: The transformed system.
- `input_vars`: A vector of input variables corresponding to the input analysis points.
- `output_vars`: A vector of output variables corresponding to the output analysis points.
"""
function linearization_ap_transform(
        sys,
        inputs::Union{Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint}},
        outputs, loop_openings
    )
    loop_openings = Set(map(nameof, canonicalize_ap(sys, loop_openings)))
    inputs = canonicalize_ap(sys, inputs)
    outputs = canonicalize_ap(sys, outputs)
    input_vars = []
    for input in inputs
        if nameof(input) in loop_openings
            delete!(loop_openings, nameof(input))
            sys, (input_var,) = apply_transformation(Break(input, true, true), sys)
        else
            sys, (input_var,) = apply_transformation(PerturbOutput(input), sys)
        end
        push!(input_vars, input_var)
    end
    output_vars = []
    for output in outputs
        if output isa AnalysisPoint
            sys, (output_var,) = apply_transformation(AddVariable(output), sys)
            sys, (input_var,) = apply_transformation(GetInput(output), sys)
            @set! sys.eqs = [get_eqs(sys); output_var ~ input_var]
        else
            output_var = output
        end
        push!(output_vars, output_var)
    end
    sys = handle_loop_openings(sys, map(AnalysisPoint, collect(loop_openings)))
    return sys, input_vars, output_vars
end

function linearization_function(
        sys::AbstractSystem,
        inputs::Union{Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint}},
        outputs; loop_openings = [], system_modifier = identity, kwargs...
    )
    sys, input_vars,
        output_vars = linearization_ap_transform(
        sys, inputs, outputs, loop_openings
    )
    return linearization_function(system_modifier(sys), input_vars, output_vars; kwargs...)
end

@doc """
    get_sensitivity(sys, ap::AnalysisPoint; kwargs)
    get_sensitivity(sys, ap_name::Symbol; kwargs)

Compute the sensitivity function in analysis point `ap`. The sensitivity function is obtained by introducing an infinitesimal perturbation `d` at the input of `ap`, linearizing the system and computing the transfer function between `d` and the output of `ap`.

# Arguments:

  - `kwargs`: Are sent to `ModelingToolkit.linearize`

See also [`get_comp_sensitivity`](@ref), [`get_looptransfer`](@ref).
""" get_sensitivity

@doc """
    get_comp_sensitivity(sys, ap::AnalysisPoint; kwargs)
    get_comp_sensitivity(sys, ap_name::Symbol; kwargs)

Compute the complementary sensitivity function in analysis point `ap`. The complementary sensitivity function is obtained by introducing an infinitesimal perturbation `d` at the output of `ap`, linearizing the system and computing the transfer function between `d` and the input of `ap`.

# Arguments:

  - `kwargs`: Are sent to `ModelingToolkit.linearize`

See also [`get_sensitivity`](@ref), [`get_looptransfer`](@ref).
""" get_comp_sensitivity

@doc """
    get_looptransfer(sys, ap::AnalysisPoint; kwargs)
    get_looptransfer(sys, ap_name::Symbol; kwargs)

Compute the (linearized) loop-transfer function in analysis point `ap`, from `ap.out` to `ap.in`.

!!! info "Negative feedback"

    Feedback loops often use negative feedback, and the computed loop-transfer function will in this case have the negative feedback included. Standard analysis tools often assume a loop-transfer function without the negative gain built in, and the result of this function may thus need negation before use.

# Arguments:

  - `kwargs`: Are sent to `ModelingToolkit.linearize`

See also [`get_sensitivity`](@ref), [`get_comp_sensitivity`](@ref), [`open_loop`](@ref).
""" get_looptransfer
