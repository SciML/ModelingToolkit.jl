"""
    $(TYPEDEF)
    AnalysisPoint(input, name::Symbol, outputs::Vector)

Create an AnalysisPoint for linear analysis. Analysis points can be created by calling

```
connect(out, :ap_name, in...)
```

Where `out` is the output being connected to the inputs `in...`. All involved
connectors (input and outputs) are required to either have an unknown named
`u` or a single unknown, all of which should have the same size.

See also [`get_sensitivity`](@ref), [`get_comp_sensitivity`](@ref), [`get_looptransfer`](@ref), [`open_loop`](@ref)

# Fields

$(TYPEDFIELDS)

# Example

```julia
using ModelingToolkitBase
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitBase: t_nounits as t

@named P = FirstOrder(k = 1, T = 1)
@named C = Gain(; k = -1)
t = ModelingToolkitBase.get_iv(P)

eqs = [connect(P.output, C.input)
       connect(C.output, :plant_input, P.input)]
sys = System(eqs, t, systems = [P, C], name = :feedback_system)

matrices_S, _ = get_sensitivity(sys, :plant_input) # Compute the matrices of a state-space representation of the (input) sensitivity function.
matrices_T, _ = get_comp_sensitivity(sys, :plant_input)
```

Continued linear analysis and design can be performed using ControlSystemsBase.jl.
Create `ControlSystemsBase.StateSpace` objects using

```julia
using ControlSystemsBase, Plots
S = ss(matrices_S...)
T = ss(matrices_T...)
bodeplot([S, T], lab = ["S" "T"])
```

The sensitivity functions obtained this way should be equivalent to the ones obtained with the code below

```julia
using ControlSystemsBase
P = tf(1.0, [1, 1])
C = 1                      # Negative feedback assumed in ControlSystems
S = sensitivity(P, C)      # or feedback(1, P*C)
T = comp_sensitivity(P, C) # or feedback(P*C)
```
"""
struct AnalysisPoint
    """
    The input to the connection. In the context of ModelingToolkitStandardLibrary.jl,
    this is a `RealOutput` connector.
    """
    input::Any
    """
    The name of the analysis point.
    """
    name::Symbol
    """
    The outputs of the connection. In the context of ModelingToolkitStandardLibrary.jl,
    these are all `RealInput` connectors.
    """
    outputs::Union{Nothing, Vector{System}, Vector{SymbolicT}}

    function AnalysisPoint(input, name::Symbol, outputs; verbose = true)
        # input to analysis point should be an output variable
        if verbose && input !== nothing
            var = ap_var(input)
            isoutput(var) || ap_warning(1, name, true)
        end
        # outputs of analysis points should be input variables
        if verbose && outputs !== nothing
            for (i, output) in enumerate(outputs)
                var = ap_var(output)
                isinput(var) || ap_warning(2 + i, name, false)
            end
        end

        return new(input, name, outputs)
    end
end

function ap_warning(arg::Int, name::Symbol, should_be_output)
    causality = should_be_output ? "output" : "input"
    return @warn """
    The $(arg)-th argument to analysis point $(name) was not a $causality. This is supported in \
    order to handle inverse models, but may not be what you intended.

    If you are building a forward mode (causal), you may want to swap this argument with \
    one on the opposite side of the name of the analysis point provided to `connect`. \
    Learn more about the causality of analysis points in the docstring for `AnalysisPoint`. \
    Silence this message using `connect(out, :name, in...; warn = false)`.
    """
end

AnalysisPoint() = AnalysisPoint(nothing, Symbol(), nothing)
"""
    $(TYPEDSIGNATURES)

Create an `AnalysisPoint` with the given name, with no input or outputs specified.
"""
AnalysisPoint(name::Symbol) = AnalysisPoint(nothing, name, nothing)

Base.nameof(ap::AnalysisPoint) = ap.name

Base.show(io::IO, ap::AnalysisPoint) = show(io, MIME"text/plain"(), ap)
function Base.show(io::IO, ::MIME"text/plain", ap::AnalysisPoint)
    Symbolics.warn_load_latexify()
    if ap.input === nothing
        print(io, "0")
        return
    end
    return if get(io, :compact, false)
        print(
            io,
            "AnalysisPoint($(ap_var(ap.input)), $(ap_var.(ap.outputs)); name=$(ap.name))"
        )
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

Symbolics.hide_lhs(::AnalysisPoint) = true

"""
    $(TYPEDSIGNATURES)

Convert an `AnalysisPoint` to a standard connection.
"""
function to_connection(ap::AnalysisPoint)
    if ap.input isa System
        vs = System[ap.input]
        append!(vs, ap.outputs::Vector{System})
        return Connection() ~ Connection(vs)
    elseif ap.input isa SymbolicT
        vs = SymbolicT[ap.input]
        append!(vs, ap.outputs::Vector{SymbolicT})
        return Connection() ~ Connection(vs)
    else
        error("Unreachable!")
    end
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
function connect(in, ap::AnalysisPoint, outs...; verbose = true)
    return AnalysisPoint() ~ AnalysisPoint(unwrap(in), ap.name, collect(unwrap.(outs)); verbose)
end

"""
    connect(output_connector, ap_name::Symbol, input_connector; verbose = true)
    connect(output_connector, ap::AnalysisPoint, input_connector; verbose = true)

Connect `output_connector` and `input_connector` with an [`AnalysisPoint`](@ref) inbetween.
The incoming connection `output_connector` is expected to be an output connector (for
example, `ModelingToolkitStandardLibrary.Blocks.RealOutput`), and vice versa.

*PLEASE NOTE*: The connection is assumed to be *causal*, meaning that

```julia
@named P = FirstOrder(k = 1, T = 1)
@named C = Gain(; k = -1)
connect(C.output, :plant_input, P.input)
```

is correct, whereas

```julia
connect(P.input, :plant_input, C.output)
```

typically is not (unless the model is an inverse model).

# Arguments

- `output_connector`: An output connector
- `input_connector`: An input connector
- `ap`: An explicitly created [`AnalysisPoint`](@ref)
- `ap_name`: If a name is given, an [`AnalysisPoint`](@ref) with the given name will be
  created automatically.

# Keyword arguments

- `verbose`: Warn if an input is connected to an output (reverse causality). Silence this
  warning if you are analyzing an inverse model.
"""
function connect(in::AbstractSystem, name::Symbol, out, outs...; verbose = true)
    return AnalysisPoint() ~ AnalysisPoint(in, name, System[out; collect(outs)]; verbose)
end

function connect(
        in::ConnectableSymbolicT, name::Symbol, out::ConnectableSymbolicT,
        outs::ConnectableSymbolicT...; verbose = true
    )
    allvars = SymbolicT[]
    push!(allvars, unwrap(in))
    push!(allvars, unwrap(out))
    for var in outs
        push!(allvars, unwrap(var))
    end
    validate_causal_variables_connection(allvars)
    return AnalysisPoint() ~ AnalysisPoint(
        allvars[1], name, allvars[2:end]; verbose
    )
end

"""
    $(TYPEDSIGNATURES)

Return all the namespaces in `name`. Namespaces should be separated by `.` or
`$NAMESPACE_SEPARATOR`.
"""
function namespace_hierarchy(name::Symbol)
    return map(
        Symbol, split(string(name), ('.', NAMESPACE_SEPARATOR))
    )
end

"""
    $(TYPEDSIGNATURES)

Remove all `AnalysisPoint`s in `sys` and any of its subsystems, replacing them by equivalent connections.
"""
function remove_analysis_points(sys::AbstractSystem)
    eqs = Equation[]
    for eq in get_eqs(sys)
        if unwrap_const(eq.lhs) isa AnalysisPoint
            push!(eqs, to_connection(unwrap_const(eq.rhs)::AnalysisPoint))
        else
            push!(eqs, eq)
        end
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
function ap_var(sys::AbstractSystem)
    if hasproperty(sys, :u)
        return unwrap(sys.u)
    end
    x = unknowns(sys)
    length(x) == 1 && return renamespace(sys, x[1])
    error("Could not determine the analysis-point variable in system $(nameof(sys)). To use an analysis point, apply it to a connection between causal blocks which have a variable named `u` or a single unknown of the same size.")
end

"""
    $(TYPEDSIGNATURES)

For an `AnalysisPoint` involving causal variables. Simply return the variable.
"""
function ap_var(var::ConnectableSymbolicT)
    return unwrap(var)::SymbolicT
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
function modify_nested_subsystem(fn, root::AbstractSystem, target::AbstractSystem)
    return modify_nested_subsystem(
        fn, root, nameof(target)
    )
end
"""
    $(TYPEDSIGNATURES)

Apply the modification to the system containing the namespaced analysis point `target`.
"""
function modify_nested_subsystem(fn, root::AbstractSystem, target::AnalysisPoint)
    return modify_nested_subsystem(
        fn, root, @view namespace_hierarchy(nameof(target))[1:(end - 1)]
    )
end
"""
    $(TYPEDSIGNATURES)

Apply the modification to the nested subsystem of `root` whose namespaced name matches
the provided name `target`. The namespace separator in `target` should be `.` or
`$NAMESPACE_SEPARATOR`. The `target` may include `nameof(root)` as the first namespace.
"""
function modify_nested_subsystem(fn, root::AbstractSystem, target::Symbol)
    return modify_nested_subsystem(
        fn, root, namespace_hierarchy(target)
    )
end

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
        fn, root::AbstractSystem, hierarchy::AbstractVector{Symbol}
    )
    # no hierarchy, so just apply to the root
    if isempty(hierarchy)
        return fn(root)
    end
    # ignore the name of the root
    if nameof(root) != hierarchy[1]
        error(
            """
            Invalid analysis point name `$(join(hierarchy, NAMESPACE_SEPARATOR))`. The name
            must include the name of the root system `$(nameof(root))`. This typically happens
            when  using an analysis point obtained by calling `getproperty` on a system marked
            as `complete` to linearize a system that is not marked as `complete`.
            """
        )
    end
    hierarchy = @view hierarchy[2:end]

    # recursive helper function which does the searching and modification
    function _helper(sys::AbstractSystem, i::Int)
        if i > length(hierarchy)
            # we reached past the end, so everything matched and
            # `sys` is the system to modify.
            sys, vars = fn(sys)
        else
            # find the subsystem with the given name and error otherwise
            cur = hierarchy[i]
            idx = findfirst(subsys -> nameof(subsys) == cur, get_systems(sys))
            idx === nothing &&
                error("System $(join([nameof(root); hierarchy[1:(i - 1)]], '.')) does not have a subsystem named $cur.")

            # recurse into new subsystem
            newsys, vars = _helper(get_systems(sys)[idx], i + 1)
            # update this system with modified subsystem
            @set! sys.systems[idx] = newsys
        end
        # only namespace variables from inner systems
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
function analysis_point_index(sys::AbstractSystem, ap::AnalysisPoint)
    return analysis_point_index(
        sys, nameof(ap)
    )
end
"""
    $(TYPEDSIGNATURES)

Search for the analysis point with the given `name` in `get_eqs(sys)`.
"""
function analysis_point_index(sys::AbstractSystem, name::Symbol)
    name = namespace_hierarchy(name)[end]
    return findfirst(get_eqs(sys)) do eq
        value(eq.lhs) isa AnalysisPoint && nameof(value(eq.rhs)::AnalysisPoint) == name
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
    if symbolic_type(var) == ArraySymbolic()
        T = eltype(symtype(var))
        pvar = unwrap(only(@variables $name(iv)[SU.shape(var)...]::T))
        default = zeros(eltype(symtype(var)), size(var))
    else
        T = symtype(var)
        pvar = unwrap(only(@variables $name(iv)::T))
        default = zero(T)
    end
    return pvar, default
end

function with_analysis_point_ignored(sys::AbstractSystem, ap::AnalysisPoint)
    has_ignored_connections(sys) || return sys
    ignored = get_ignored_connections(sys)
    if ignored === nothing
        ignored = Connection[]
    else
        ignored = copy(ignored)
    end
    if ap.outputs === nothing
        error("Empty analysis point")
    end

    push!(ignored, Connection([unwrap(ap.input); unwrap.(ap.outputs)]))

    return @set sys.ignored_connections = ignored
end

#### PRIMITIVE TRANSFORMATIONS

const DOC_WILL_REMOVE_AP = """
Note that this transformation will remove `ap`, causing any subsequent transformations \
referring to it to fail.\
"""

const DOC_ADDED_VARIABLE = """
The added variable(s) will have a default of zero, of the appropriate type and size.\
"""

"""
    $(TYPEDEF)

A transformation which breaks the connection referred to by `ap`.
`apply_transformation` returns the vector of output variables as auxiliary information.
`outputs_to_params` controls whether the output variables are turned into parameters.

$DOC_WILL_REMOVE_AP

## Fields

$(TYPEDFIELDS)
"""
struct Break <: AnalysisPointTransformation
    """
    The analysis point to break.
    """
    ap::AnalysisPoint
    """
    Whether to turn outputs of the analysis point into parameters.
    """
    outputs_to_params::Bool
end

Break(ap::AnalysisPoint) = Break(ap, false)

function apply_transformation(tf::Break, sys::AbstractSystem)
    return modify_nested_subsystem(sys, tf.ap) do breaksys
        # get analysis point
        ap_idx = analysis_point_index(breaksys, tf.ap)
        ap_idx === nothing &&
            error("Analysis point $(nameof(tf.ap)) not found in system $(nameof(sys)).")
        breaksys_eqs = copy(get_eqs(breaksys))
        @set! breaksys.eqs = breaksys_eqs

        ap = value(breaksys_eqs[ap_idx].rhs)::AnalysisPoint
        deleteat!(breaksys_eqs, ap_idx)

        breaksys = with_analysis_point_ignored(breaksys, ap)

        out_vars = SymbolicT[]
        outs = ap.outputs::Union{Vector{System}, Vector{SymbolicT}}
        if outs isa Vector{System}
            for out in outs
                if tf.outputs_to_params
                    breaksys, (out_var,) = modify_nested_subsystem(breaksys, renamespace(breaksys, out)) do outsys
                        new_dvs = copy(get_unknowns(outsys))
                        @set! outsys.unknowns = new_dvs
                        new_ps = copy(get_ps(outsys))
                        @set! outsys.ps = new_ps

                        out_var = ap_var(toggle_namespacing(outsys, false))
                        deleteat!(new_dvs, findfirst(isequal(out_var), new_dvs)::Int)
                        push!(new_ps, out_var)

                        return outsys, (out_var,)
                    end
                else
                    out_var = ap_var(out)
                end
                push!(out_vars, out_var)
            end
        else
            for out in outs
                out_var = ap_var(out)
                if tf.outputs_to_params
                    new_dvs = copy(get_unknowns(breaksys))
                    @set! breaksys.unknowns = new_dvs
                    new_ps = copy(get_ps(breaksys))
                    @set! breaksys.ps = new_ps
                    deleteat!(new_dvs, findfirst(isequal(out_var), new_dvs)::Int)
                    push!(new_ps, out_var)
                end
                push!(out_vars, ap_var(out))
            end
        end
        return breaksys, (out_vars,)
    end
end

"""
    $(TYPEDEF)

A transformation which returns the variable corresponding to the input of the analysis
point. Does not modify the system.

## Fields

$(TYPEDFIELDS)
"""
struct GetInput <: AnalysisPointTransformation
    """
    The analysis point to get the input of.
    """
    ap::AnalysisPoint
end

function apply_transformation(tf::GetInput, sys::AbstractSystem)
    return modify_nested_subsystem(sys, tf.ap) do ap_sys
        # get the analysis point
        ap_idx = analysis_point_index(ap_sys, tf.ap)
        ap_idx === nothing &&
            error("Analysis point $(nameof(tf.ap)) not found in system $(nameof(sys)).")
        # get the analysis point
        ap_sys_eqs = get_eqs(ap_sys)
        ap = value(ap_sys_eqs[ap_idx].rhs)::AnalysisPoint

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
as the analysis point and will be the second variable in the tuple of new variables returned
from `apply_transformation`.

$DOC_WILL_REMOVE_AP

$DOC_ADDED_VARIABLE

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
    return modify_nested_subsystem(sys, tf.ap) do ap_sys
        # get analysis point
        ap_idx = analysis_point_index(ap_sys, tf.ap)
        ap_idx === nothing &&
            error("Analysis point $(nameof(tf.ap)) not found in system $(nameof(sys)).")
        # modified equations
        ap_sys_eqs = copy(get_eqs(ap_sys))
        @set! ap_sys.eqs = ap_sys_eqs
        ap = value(ap_sys_eqs[ap_idx].rhs)::AnalysisPoint
        # remove analysis point
        deleteat!(ap_sys_eqs, ap_idx)
        ap_sys = with_analysis_point_ignored(ap_sys, ap)

        # add equations involving new variable
        ap_ivar = ap_var(ap.input)
        new_var, new_def = get_analysis_variable(ap_ivar, nameof(ap), get_iv(sys))
        for outsys in ap.outputs
            push!(ap_sys_eqs, ap_var(outsys) ~ ap_ivar + wrap(new_var))
        end
        # add variable
        unks = copy(get_unknowns(ap_sys))
        push!(unks, new_var)
        @set! ap_sys.unknowns = unks
        # add default
        defs = copy(get_initial_conditions(ap_sys))
        defs[new_var] = new_def
        @set! ap_sys.initial_conditions = defs

        tf.with_output || return ap_sys, (new_var,)

        # add output variable, equation, default
        out_var,
            out_def = get_analysis_variable(
            ap_ivar, nameof(ap), get_iv(sys); perturb = false
        )
        push!(ap_sys_eqs, out_var ~ ap_ivar + wrap(new_var))
        push!(unks, out_var)

        return ap_sys, (new_var, out_var)
    end
end

"""
    $(TYPEDEF)

A transformation which adds a variable named `name` to the system containing the analysis
point `ap`. $DOC_ADDED_VARIABLE

# Fields

$(TYPEDFIELDS)
"""
struct AddVariable <: AnalysisPointTransformation
    """
    The analysis point in the system to modify, and whose input should be used as the
    template for the new variable.
    """
    ap::AnalysisPoint
    """
    The name of the added variable.
    """
    name::Symbol
end

"""
    $(TYPEDSIGNATURES)

Add a new variable to the system containing analysis point `ap` with the same name as the
analysis point.
"""
AddVariable(ap::AnalysisPoint) = AddVariable(ap, nameof(ap))

function apply_transformation(tf::AddVariable, sys::AbstractSystem)
    return modify_nested_subsystem(sys, tf.ap) do ap_sys
        # get analysis point
        ap_idx = analysis_point_index(ap_sys, tf.ap)
        ap_idx === nothing &&
            error("Analysis point $(nameof(tf.ap)) not found in system $(nameof(sys)).")
        ap_sys_eqs = get_eqs(ap_sys)
        ap = value(ap_sys_eqs[ap_idx].rhs)::AnalysisPoint

        # add equations involving new variable
        ap_ivar = ap_var(ap.input)
        new_var,
            new_def = get_analysis_variable(
            ap_ivar, tf.name, get_iv(sys); perturb = false
        )
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
The returned added variables are `(du, u)` where `du` is the perturbation added to the
input, and `u` is the output after perturbation.

$DOC_WILL_REMOVE_AP

$DOC_ADDED_VARIABLE
"""
SensitivityTransform(ap::AnalysisPoint) = PerturbOutput(ap, true)

"""
    $(TYPEDEF)

A transformation to enable calculating the complementary sensitivity function about the
analysis point `ap`. The returned added variables are `(du, u)` where `du` is the
perturbation added to the outputs and `u` is the input to the analysis point.

$DOC_WILL_REMOVE_AP

$DOC_ADDED_VARIABLE

# Fields

$(TYPEDFIELDS)
"""
struct ComplementarySensitivityTransform <: AnalysisPointTransformation
    """
    The analysis point to modify.
    """
    ap::AnalysisPoint
end

function apply_transformation(cst::ComplementarySensitivityTransform, sys::AbstractSystem)
    sys, (u,) = apply_transformation(GetInput(cst.ap), sys)
    sys,
        (du,) = apply_transformation(
        AddVariable(
            cst.ap, Symbol(namespace_hierarchy(nameof(cst.ap))[end], :_comp_sens_du)
        ),
        sys
    )
    sys, (_du,) = apply_transformation(PerturbOutput(cst.ap), sys)

    # `PerturbOutput` adds the equation `input + _du ~ output`
    # but comp sensitivity wants `output + du ~ input`. Thus, `du ~ -_du`.
    eqs = copy(get_eqs(sys))
    @set! sys.eqs = eqs
    push!(eqs, du ~ -wrap(_du))

    defs = copy(get_initial_conditions(sys))
    @set! sys.initial_conditions = defs
    defs[du] = -wrap(_du)
    return sys, (du, u)
end

"""
    $(TYPEDEF)

A transformation to enable calculating the loop transfer function about the analysis point
`ap`. The returned added variables are `(du, u)` where `du` feeds into the outputs of `ap`
and `u` is the input of `ap`.

$DOC_WILL_REMOVE_AP

$DOC_ADDED_VARIABLE

# Fields

$(TYPEDFIELDS)
"""
struct LoopTransferTransform <: AnalysisPointTransformation
    """
    The analysis point to modify.
    """
    ap::AnalysisPoint
end

function apply_transformation(tf::LoopTransferTransform, sys::AbstractSystem)
    sys, (u,) = apply_transformation(GetInput(tf.ap), sys)
    sys, (dus,) = apply_transformation(Break(tf.ap), sys)
    return sys, (dus, u)
end

"""
    $(TYPEDSIGNATURES)

A utility function to get the "canonical" form of a list of analysis points. Always returns
a list of values. Any value that cannot be turned into an `AnalysisPoint` (i.e. isn't
already an `AnalysisPoint` or `Symbol`) is simply wrapped in an array. `Symbol` names of
`AnalysisPoint`s are namespaced with `sys`.
"""
canonicalize_ap(sys::AbstractSystem, ap::Symbol) = [AnalysisPoint(renamespace(sys, ap))]
function canonicalize_ap(sys::AbstractSystem, ap::AnalysisPoint)
    if does_namespacing(sys)
        return [ap]
    else
        return [renamespace(sys, ap)]
    end
end
canonicalize_ap(sys::AbstractSystem, ap) = [ap]
function canonicalize_ap(sys::AbstractSystem, aps::Vector)
    return mapreduce(Base.Fix1(canonicalize_ap, sys), vcat, aps; init = [])
end

"""
    $(TYPEDSIGNATURES)

Apply `LoopTransferTransform` to the analysis point `ap` and return the
result of `apply_transformation`.

# Keyword Arguments

- `system_modifier`: a function which takes the modified system and returns a new system
  with any required further modifications performed.
"""
function open_loop(sys, ap::Union{Symbol, AnalysisPoint}; system_modifier = identity)
    ap = only(canonicalize_ap(sys, ap))
    tf = LoopTransferTransform(ap)
    sys, vars = apply_transformation(tf, sys)
    return system_modifier(sys), vars
end

"""
    generate_control_function(sys::ModelingToolkitBase.AbstractSystem, input_ap_name::Union{Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint}}, dist_ap_name::Union{Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint}}; system_modifier = identity, known_disturbance_inputs = nothing, kwargs)

When called with analysis points as input arguments, we assume that all analysis points corresponds to connections that should be opened (broken). The use case for this is to get rid of input signal blocks, such as `Step` or `Sine`, since these are useful for simulation but are not needed when using the plant model in a controller or state estimator.

The `known_disturbance_inputs` keyword accepts analysis points (or symbols) for disturbances that should be treated as known (i.e., provided as an additional function argument `w`). The loops at these analysis points are opened and the resulting variables are forwarded to the base `generate_control_function`.
"""
function generate_control_function(
        sys::ModelingToolkitBase.AbstractSystem, input_ap_name::Union{
            Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint},
        },
        dist_ap_name::Union{
            Nothing, Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint},
        } = nothing;
        system_modifier = identity,
        known_disturbance_inputs::Union{
            Nothing, Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint},
        } = nothing,
        kwargs...
    )
    input_ap_name = canonicalize_ap(sys, input_ap_name)
    u = []
    for input_ap in input_ap_name
        sys, (du, _) = open_loop(sys, input_ap)
        push!(u, du)
    end

    # Open loops for known disturbance APs and collect their variables.
    # Use append! (not push!) because open_loop returns a Vector{SymbolicT}
    # from Break, and the base function expects a flat vector of symbolic vars
    known_dist_vars = nothing
    if known_disturbance_inputs !== nothing
        known_disturbance_inputs = canonicalize_ap(sys, known_disturbance_inputs)
        known_dist_vars = []
        for kd_ap in known_disturbance_inputs
            sys, (du, _) = open_loop(sys, kd_ap)
            append!(known_dist_vars, du)
        end
    end

    # Open loops for unknown disturbance APs (positional dist_ap_name)
    d = nothing
    if dist_ap_name !== nothing
        dist_ap_name = canonicalize_ap(sys, dist_ap_name)
        d = []
        for dist_ap in dist_ap_name
            sys, (du, _) = open_loop(sys, dist_ap)
            push!(d, du)
        end
    end

    # Always pass disturbance_inputs explicitly to prevent the base function
    # from defaulting to disturbances(sys), which would overlap with variables
    # we already opened and get substituted to zero.
    return ModelingToolkitBase.generate_control_function(
        system_modifier(sys), u, d;
        known_disturbance_inputs = known_dist_vars, kwargs...)
end
