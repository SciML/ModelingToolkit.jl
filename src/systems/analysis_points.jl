"""
    $(TYPEDEF)
    $(TYPEDSIGNATURES)

Create an AnalysisPoint for linear analysis. Analysis points can be created by calling

```
connect(in, :ap_name, out...)
```

Where `in` is the input to the connection, and `out...` are the outputs. In the context of
ModelingToolkitStandardLibrary.jl, `in` is a `RealOutput` connector and `out...` are all
`RealInput` connectors. All involved connectors are required to either have an unknown named
`u` or a single unknown, all of which should have the same size.
"""
struct AnalysisPoint
    input::Any
    name::Symbol
    outputs::Union{Nothing, Vector{Any}}
end

AnalysisPoint() = AnalysisPoint(nothing, Symbol(), nothing)
AnalysisPoint(name::Symbol) = AnalysisPoint(nothing, name, nothing)

Base.nameof(ap::AnalysisPoint) = ap.name

Base.show(io::IO, ap::AnalysisPoint) = show(io, MIME"text/plain"(), ap)
function Base.show(io::IO, ::MIME"text/plain", ap::AnalysisPoint)
    if ap.input === nothing
        print(io, "0")
        return
    end
    if get(io, :compact, false)
        print(io,
            "AnalysisPoint($(ap_var(ap.input)), $(ap_var.(ap.outputs)); name=$(ap.name))")
    else
        print(io, "AnalysisPoint(")
        printstyled(io, ap.name, color = :cyan)
        if ap.input !== nothing && ap.outputs !== nothing
            print(io, " from ")
            printstyled(io, ap_var(ap.input), color = :green)
            print(io, " to ")
            if length(ap.outputs) == 1
                printstyled(io, ap_var(ap.outputs[1]), color = :blue)
            else
                printstyled(io, "[", join(ap_var.(ap.outputs), ", "), "]", color = :blue)
            end
        end
        print(io, ")")
    end
end

"""
    $(TYPEDSIGNATURES)

Convert an `AnalysisPoint` to a standard connection.
"""
function to_connection(ap::AnalysisPoint)
    return connect(ap.input, ap.outputs...)
end

"""
    $(TYPEDSIGNATURES)

Namespace an `AnalysisPoint` by namespacing the involved systems and the name of the point.
"""
function renamespace(sys, ap::AnalysisPoint)
    return AnalysisPoint(
        ap.input === nothing ? nothing : renamespace(sys, ap.input),
        renamespace(sys, ap.name),
        ap.outputs === nothing ? nothing : map(Base.Fix1(renamespace, sys), ap.outputs)
    )
end

# create analysis points via `connect`
function Symbolics.connect(in, ap::AnalysisPoint, outs...)
    return AnalysisPoint() ~ AnalysisPoint(in, ap.name, collect(outs))
end

"""
    $(TYPEDSIGNATURES)

Create an `AnalysisPoint` connection connecting `in` to `outs...`.
"""
function Symbolics.connect(in::AbstractSystem, name::Symbol, out, outs...)
    return AnalysisPoint() ~ AnalysisPoint(in, name, [out; collect(outs)])
end

"""
    $(TYPEDSIGNATURES)

Return all the namespaces in `name`. Namespaces should be separated by `.` or
`$NAMESPACE_SEPARATOR`.
"""
namespace_hierarchy(name::Symbol) = map(
    Symbol, split(string(name), ('.', NAMESPACE_SEPARATOR)))

"""
    $(TYPEDSIGNATURES)

Remove all `AnalysisPoint`s in `sys` and any of its subsystems, replacing them by equivalent connections.
"""
function remove_analysis_points(sys::AbstractSystem)
    eqs = map(get_eqs(sys)) do eq
        eq.lhs isa AnalysisPoint ? to_connection(eq.rhs) : eq
    end
    @set! sys.eqs = eqs
    @set! sys.systems = map(remove_analysis_points, get_systems(sys))

    return sys
end

"""
    $(TYPEDSIGNATURES)

Given a system involved in an `AnalysisPoint`, get the variable to be used in the
connection. This is the variable named `u` if present, and otherwise the only
variable in the system. If the system does not have a variable named `u` and
contains multiple variables, throw an error.
"""
function ap_var(sys)
    if hasproperty(sys, :u)
        return sys.u
    end
    x = unknowns(sys)
    length(x) == 1 && return renamespace(sys, x[1])
    error("Could not determine the analysis-point variable in system $(nameof(sys)). To use an analysis point, apply it to a connection between causal blocks which have a variable named `u` or a single unknown of the same size.")
end

"""
    $(TYPEDEF)

The supertype of all transformations that can be applied to an `AnalysisPoint`. All
concrete subtypes must implement `apply_transformation`.
"""
abstract type AnalysisPointTransformation end

"""
    apply_transformation(tf::AnalysisPointTransformation, sys::AbstractSystem)

Apply the given analysis point transformation `tf` to the system `sys`. Throw an error if
any analysis points referred to in `tf` are not present in `sys`. Return a tuple
containing the modified system as the first element, and a tuple of the additional
variables added by the transformation as the second element.
"""
function apply_transformation end

"""
    $(TYPEDSIGNATURES)

Given a namespaced subsystem `target` of root system `root`, return a modified copy of
`root` with `target` modified according to `fn` alongside any extra variables added
by `fn`.

`fn` is a function which takes the instance of `target` present in the hierarchy of
`root`, and returns a 2-tuple consisting of the modified version of `target` and a tuple
of the extra variables added.
"""
modify_nested_subsystem(fn, root::AbstractSystem, target::AbstractSystem) = modify_nested_subsystem(
    fn, root, nameof(target))
"""
    $(TYPEDSIGNATURES)

Apply the modification to the system containing the namespaced analysis point `target`.
"""
modify_nested_subsystem(fn, root::AbstractSystem, target::AnalysisPoint) = modify_nested_subsystem(
    fn, root, @view namespace_hierarchy(nameof(target))[1:(end - 1)])
"""
    $(TYPEDSIGNATURES)

Apply the modification to the nested subsystem of `root` whose namespaced name matches
the provided name `target`. The namespace separator in `target` should be `.` or
`$NAMESPACE_SEPARATOR`. The `target` may include `nameof(root)` as the first namespace.
"""
modify_nested_subsystem(fn, root::AbstractSystem, target::Symbol) = modify_nested_subsystem(
    fn, root, namespace_hierarchy(target))

"""
    $(TYPEDSIGNATURES)

Apply the modification to the nested subsystem of `root` where the name of the subsystem at
each level in the hierarchy is given by elements of `hierarchy`. For example, if
`hierarchy = [:a, :b, :c]`, the system being searched for will be `root.a.b.c`. Note that
the hierarchy may include the name of the root system, in which the first element will be
ignored. For example, `hierarchy = [:root, :a, :b, :c]` also searches for `root.a.b.c`.
An empty `hierarchy` will apply the modification to `root`.
"""
function modify_nested_subsystem(
        fn, root::AbstractSystem, hierarchy::AbstractVector{Symbol})
    # no hierarchy, so just apply to the root
    if isempty(hierarchy)
        return fn(root)
    end
    # ignore the name of the root
    if nameof(root) == hierarchy[1]
        hierarchy = @view hierarchy[2:end]
    end

    # recursive helper function which does the searching and modification
    function _helper(sys::AbstractSystem, i::Int)
        # we reached past the end, so everything matched and
        # `sys` is the system to modify.
        if i > length(hierarchy)
            sys, vars = fn(sys)
        else
            cur = hierarchy[i]
            idx = findfirst(subsys -> nameof(subsys) == cur, get_systems(sys))
            idx === nothing &&
                error("System $(join([nameof(root); hierarchy[1:i-1]], '.')) does not have a subsystem named $cur.")

            newsys, vars = _helper(get_systems(sys)[idx], i + 1)
            @set! sys.systems[idx] = newsys
        end
        if i != 1
            vars = ntuple(Val(length(vars))) do i
                renamespace(sys, vars[i])
            end
        end
        return sys, vars
    end

    return _helper(root, 1)
end

"""
    $(TYPEDSIGNATURES)

Given a system `sys` and analysis point `ap`, return the index in `get_eqs(sys)`
containing an equation which has as it's RHS an analysis point with name `nameof(ap)`.
"""
analysis_point_index(sys::AbstractSystem, ap::AnalysisPoint) = analysis_point_index(
    sys, nameof(ap))
"""
    $(TYPEDSIGNATURES)

Search for the analysis point with the given `name` in `get_eqs(sys)`.
"""
function analysis_point_index(sys::AbstractSystem, name::Symbol)
    name = namespace_hierarchy(name)[end]
    findfirst(get_eqs(sys)) do eq
        eq.lhs isa AnalysisPoint && nameof(eq.rhs) == name
    end
end

"""
    $(TYPEDSIGNATURES)

Create a new variable of the same `symtype` and size as `var`, using `name` as the base
name for the new variable. `iv` denotes the independent variable of the system. Prefix
`d_` to the name of the new variable if `perturb == true`. Return the new symbolic
variable and the appropriate zero value for it.
"""
function get_analysis_variable(var, name, iv; perturb = true)
    var = unwrap(var)
    if perturb
        name = Symbol(:d_, name)
    end
    if Symbolics.isarraysymbolic(var)
        T = Array{eltype(symtype(var)), ndims(var)}
        pvar = unwrap(only(@variables $name(iv)::T))
        pvar = setmetadata(pvar, Symbolics.ArrayShapeCtx, Symbolics.shape(var))
        default = zeros(symtype(var), size(var))
    else
        T = symtype(var)
        pvar = unwrap(only(@variables $name(iv)::T))
        default = zero(T)
    end
    return pvar, default
end

#### PRIMITIVE TRANSFORMATIONS

"""
    $(TYPEDEF)

A transformation which breaks the connection referred to by `ap`. If `add_input == true`,
it will add a new input variable which connects to the outputs of the analysis point.
`apply_transformation` returns the new input variable (if added) as the auxiliary
information. The new input variable will have the name `Symbol(:d_, nameof(ap))`.

Note that this transformation will remove `ap`, causing any subsequent transformations
referring to it to fail.

The added variable, if present, will have a default of zero, of the appropriate type and
size.

## Fields

$(TYPEDFIELDS)
"""
struct Break <: AnalysisPointTransformation
    """
    The analysis point to break.
    """
    ap::AnalysisPoint
    """
    Whether to add a new input variable connected to all the outputs of `ap`.
    """
    add_input::Bool
end

"""
    $(TYPEDSIGNATURES)

`Break` the given analysis point `ap` without adding an input.
"""
Break(ap::AnalysisPoint) = Break(ap, false)

function apply_transformation(tf::Break, sys::AbstractSystem)
    modify_nested_subsystem(sys, tf.ap) do breaksys
        # get analysis point
        ap_idx = analysis_point_index(breaksys, tf.ap)
        ap_idx === nothing &&
            error("Analysis point $(nameof(tf.ap)) not found in system $(nameof(sys)).")
        breaksys_eqs = copy(get_eqs(breaksys))
        @set! breaksys.eqs = breaksys_eqs

        ap = breaksys_eqs[ap_idx].rhs
        deleteat!(breaksys_eqs, ap_idx)

        tf.add_input || return sys, ()

        ap_ivar = ap_var(ap.input)
        new_var, new_def = get_analysis_variable(ap_ivar, nameof(ap), get_iv(sys))
        for outsys in ap.outputs
            push!(breaksys_eqs, ap_var(outsys) ~ new_var)
        end
        defs = copy(get_defaults(breaksys))
        defs[new_var] = new_def
        @set! breaksys.defaults = defs
        unks = copy(get_unknowns(breaksys))
        push!(unks, new_var)
        @set! breaksys.unknowns = unks

        return breaksys, (new_var,)
    end
end

"""
    $(TYPEDEF)

A transformation which returns the variable corresponding to the input of the analysis
point.

`apply_transformation` returns the variable as auxiliary information.

## Fields

$(TYPEDFIELDS)
"""
struct GetInput <: AnalysisPointTransformation
    """
    The analysis point to add the output to.
    """
    ap::AnalysisPoint
end

function apply_transformation(tf::GetInput, sys::AbstractSystem)
    modify_nested_subsystem(sys, tf.ap) do ap_sys
        # get the analysis point
        ap_idx = analysis_point_index(ap_sys, tf.ap)
        ap_idx === nothing &&
            error("Analysis point $(nameof(tf.ap)) not found in system $(nameof(sys)).")
        # get the anlysis point
        ap_sys_eqs = copy(get_eqs(ap_sys))
        ap = ap_sys_eqs[ap_idx].rhs

        # input variable
        ap_ivar = ap_var(ap.input)
        return ap_sys, (ap_ivar,)
    end
end

"""
    $(TYPEDEF)

A transformation that creates a new input variable which is added to the input of
the analysis point before connecting to the outputs. The new variable will have the name
`Symbol(:d_, nameof(ap))`.

If `with_output == true`, also creates an additional new variable which has the value
provided to the outputs after the above modification. This new variable has the same name
as the analysis point.

`apply_transformation` returns a 1-tuple of the perturbation variable if
`with_output == false` and a 2-tuple of the perturbation variable and output variable if
`with_output == true`.

Removes the analysis point `ap`, so any subsequent transformations requiring it will fail.

The added variable(s) will have a default of zero, of the appropriate type and size.

## Fields

$(TYPEDFIELDS)
"""
struct PerturbOutput <: AnalysisPointTransformation
    """
    The analysis point to modify
    """
    ap::AnalysisPoint
    """
    Whether to add an additional output variable.
    """
    with_output::Bool
end

"""
    $(TYPEDSIGNATURES)

Add an input without an additional output variable.
"""
PerturbOutput(ap::AnalysisPoint) = PerturbOutput(ap, false)

function apply_transformation(tf::PerturbOutput, sys::AbstractSystem)
    modify_nested_subsystem(sys, tf.ap) do ap_sys
        # get analysis point
        ap_idx = analysis_point_index(ap_sys, tf.ap)
        ap_idx === nothing &&
            error("Analysis point $(nameof(tf.ap)) not found in system $(nameof(sys)).")
        # modified quations
        ap_sys_eqs = copy(get_eqs(ap_sys))
        @set! ap_sys.eqs = ap_sys_eqs
        ap = ap_sys_eqs[ap_idx].rhs
        # remove analysis point
        deleteat!(ap_sys_eqs, ap_idx)

        # add equations involving new variable
        ap_ivar = ap_var(ap.input)
        new_var, new_def = get_analysis_variable(ap_ivar, nameof(ap), get_iv(sys))
        for outsys in ap.outputs
            push!(ap_sys_eqs, ap_var(outsys) ~ ap_ivar + new_var)
        end
        # add variable
        unks = copy(get_unknowns(ap_sys))
        push!(unks, new_var)
        @set! ap_sys.unknowns = unks
        # add default
        defs = copy(get_defaults(ap_sys))
        defs[new_var] = new_def
        @set! ap_sys.defaults = defs

        tf.with_output || return ap_sys, (new_var,)

        # add output variable, equation, default
        out_var, out_def = get_analysis_variable(
            ap_ivar, nameof(ap), get_iv(sys); perturb = false)
        defs[out_var] = out_def
        push!(ap_sys_eqs, out_var ~ ap_ivar + new_var)
        push!(unks, out_var)

        return ap_sys, (new_var, out_var)
    end
end

"""
    $(TYPEDSIGNATURES)

A transformation which adds a variable named `name` to the system containing the analysis
point `ap`. The added variable has the same type and size as the input of the analysis
point.
"""
struct AddVariable <: AnalysisPointTransformation
    ap::AnalysisPoint
    name::Symbol
end

AddVariable(ap::AnalysisPoint) = AddVariable(ap, nameof(ap))

function apply_transformation(tf::AddVariable, sys::AbstractSystem)
    modify_nested_subsystem(sys, tf.ap) do ap_sys
        # get analysis point
        ap_idx = analysis_point_index(ap_sys, tf.ap)
        ap_idx === nothing &&
            error("Analysis point $(nameof(tf.ap)) not found in system $(nameof(sys)).")
        ap_sys_eqs = copy(get_eqs(ap_sys))
        ap = ap_sys_eqs[ap_idx].rhs

        # add equations involving new variable
        ap_ivar = ap_var(ap.input)
        new_var, new_def = get_analysis_variable(
            ap_ivar, tf.name, get_iv(sys); perturb = false)
        # add variable
        unks = copy(get_unknowns(ap_sys))
        push!(unks, new_var)
        @set! ap_sys.unknowns = unks
        return ap_sys, (new_var,)
    end
end

#### DERIVED TRANSFORMATIONS

"""
    $(TYPEDSIGNATURES)

A transformation enable calculating the sensitivity function about the analysis point `ap`.
`apply_transformation` returns a 2-tuple `du, u` as auxiliary information.

Removes the analysis point `ap`, so any subsequent transformations requiring it will fail.

The added variables will have a default of zero, of the appropriate type and size.
"""
SensitivityTransform(ap::AnalysisPoint) = PerturbOutput(ap, true)

"""
    $(TYPEDEF)

A transformation to enable calculating the complementary sensitivity function about the
analysis point `ap`. `apply_transformation` returns a 2-tuple `du, u` as auxiliary
information.

Removes the analysis point `ap`, so any subsequent transformations requiring it will fail.

The added variables will have a default of zero, of the appropriate type and size.
"""
struct ComplementarySensitivityTransform <: AnalysisPointTransformation
    ap::AnalysisPoint
end

function apply_transformation(cst::ComplementarySensitivityTransform, sys::AbstractSystem)
    sys, (u,) = apply_transformation(GetInput(cst.ap), sys)
    sys, (du,) = apply_transformation(AddVariable(cst.ap, Symbol(:comp_sens_du)), sys)
    sys, (_du,) = apply_transformation(PerturbOutput(cst.ap), sys)

    # `PerturbOutput` adds the equation `input + _du ~ output`
    # but comp sensitivity wants `output + du ~ input`. Thus, `du ~ -_du`.
    eqs = copy(get_eqs(sys))
    @set! sys.eqs = eqs
    push!(eqs, du ~ -_du)

    defs = copy(get_defaults(sys))
    @set! sys.defaults = defs
    defs[du] = -_du
    return sys, (du, u)
end

struct LoopTransferTransform <: AnalysisPointTransformation
    ap::AnalysisPoint
end

function apply_transformation(tf::LoopTransferTransform, sys::AbstractSystem)
    sys, (u,) = apply_transformation(GetInput(tf.ap), sys)
    sys, (du,) = apply_transformation(Break(tf.ap, true), sys)
    return sys, (du, u)
end

### TODO: Move these

canonicalize_ap(ap::Symbol) = [AnalysisPoint(ap)]
canonicalize_ap(ap::AnalysisPoint) = [ap]
canonicalize_ap(ap) = [ap]
function canonicalize_ap(aps::Vector)
    mapreduce(canonicalize_ap, vcat, aps; init = [])
end

function handle_loop_openings(sys::AbstractSystem, aps)
    for ap in canonicalize_ap(aps)
        sys, (outvar,) = apply_transformation(Break(ap, true), sys)
        if Symbolics.isarraysymbolic(outvar)
            push!(get_eqs(sys), outvar ~ zeros(size(outvar)))
        else
            push!(get_eqs(sys), outvar ~ 0)
        end
    end
    return sys
end

function get_sensitivity_function(
        sys::AbstractSystem, aps; system_modifier = identity, loop_openings = [], kwargs...)
    sys = handle_loop_openings(sys, loop_openings)
    aps = canonicalize_ap(aps)
    dus = []
    us = []
    for ap in aps
        sys, (du, u) = apply_transformation(SensitivityTransform(ap), sys)
        push!(dus, du)
        push!(us, u)
    end
    linearization_function(system_modifier(sys), dus, us; kwargs...)
end

function get_comp_sensitivity_function(
        sys::AbstractSystem, aps; system_modifier = identity, loop_openings = [], kwargs...)
    sys = handle_loop_openings(sys, loop_openings)
    aps = canonicalize_ap(aps)
    dus = []
    us = []
    for ap in aps
        sys, (du, u) = apply_transformation(ComplementarySensitivityTransform(ap), sys)
        push!(dus, du)
        push!(us, u)
    end
    linearization_function(system_modifier(sys), dus, us; kwargs...)
end

function get_looptransfer_function(
        sys, aps; system_modifier = identity, loop_openings = [], kwargs...)
    sys = handle_loop_openings(sys, loop_openings)
    aps = canonicalize_ap(aps)
    dus = []
    us = []
    for ap in aps
        sys, (du, u) = apply_transformation(LoopTransferTransform(ap), sys)
        push!(dus, du)
        push!(us, u)
    end
    linearization_function(system_modifier(sys), dus, us; kwargs...)
end

for f in [:get_sensitivity, :get_comp_sensitivity, :get_looptransfer]
    @eval function $f(sys, ap, args...; kwargs...)
        lin_fun, ssys = $(Symbol(f, :_function))(sys, ap, args...; kwargs...)
        ModelingToolkit.linearize(ssys, lin_fun; kwargs...), ssys
    end
end

function open_loop(sys, ap::Union{Symbol, AnalysisPoint}; kwargs...)
    if ap isa Symbol
        ap = AnalysisPoint(ap)
    end
    tf = LoopTransferTransform(ap)
    return apply_transformation(tf, sys)
end

function linearization_function(sys::AbstractSystem,
        inputs::Union{Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint}},
        outputs; loop_openings = [], system_modifier = identity, kwargs...)
    sys = handle_loop_openings(sys, loop_openings)

    inputs = canonicalize_ap(inputs)
    outputs = canonicalize_ap(outputs)

    input_vars = []
    for input in inputs
        sys, (input_var,) = apply_transformation(PerturbOutput(input), sys)
        push!(input_vars, input_var)
    end
    output_vars = []
    for output in outputs
        if output isa AnalysisPoint
            sys, (output_var,) = apply_transformation(GetInput(output), sys)
            push!(output_vars, output_var)
        else
            push!(output_vars, output)
        end
    end

    return linearization_function(system_modifier(sys), input_vars, output_vars; kwargs...)
end
