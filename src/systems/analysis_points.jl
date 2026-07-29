"""
    $(TYPEDSIGNATURES)

Given a list of analysis points, break the connection for each. The output of each broken
analysis point is turned into a parameter of the system. Returns the modified system and a
vector of the variables that were turned into parameters (the "loop-opening parameters"),
whose operating-point values must be supplied by the user when linearizing.
"""
function handle_loop_openings(sys::AbstractSystem, aps)
    loop_opening_params = SymbolicT[]
    for ap in canonicalize_ap(sys, aps)
        sys, (d_vs,) = apply_transformation(Break(ap, true), sys)
        append!(loop_opening_params, d_vs)
    end
    return sys, loop_opening_params
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
    dus = SymbolicT[]
    us = SymbolicT[]
    sys, loop_opening_params = handle_loop_openings(sys, loop_openings)
    aps = canonicalize_ap(sys, aps)
    for ap in aps
        sys, (du, u) = apply_transformation(transform(ap), sys)
        du = du::Union{SymbolicT, Vector{SymbolicT}}
        u = u::Union{SymbolicT, Vector{SymbolicT}}
        if du isa SymbolicT
            push!(dus, du)
        else
            append!(dus, du)
        end
        if u isa SymbolicT
            push!(us, u)
        else
            append!(us, u)
        end
    end
    return linearization_function(system_modifier(sys), dus, us; loop_opening_params, kwargs...)
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
            allow_input_derivatives = true, op = Dict{SymbolicT, SymbolicT}(), kwargs...
        )
        lin_fun,
            ssys = $(utility_fun)(
            sys, ap, args...; loop_openings, system_modifier, op, kwargs...
        )
        mats, extras = ModelingToolkit.linearize(ssys, lin_fun; op, allow_input_derivatives)
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
- `loop_opening_params`: A vector of the variables that were turned into parameters by the
  `loop_openings` (their operating-point values must be supplied when linearizing).
"""
function linearization_ap_transform(
        sys,
        inputs::Union{Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint}},
        outputs, loop_openings
    )
    loop_openings = Set(map(nameof, canonicalize_ap(sys, loop_openings)))
    inputs = canonicalize_ap(sys, inputs)
    outputs = canonicalize_ap(sys, outputs)
    input_vars = SymbolicT[]
    for input in inputs
        if nameof(input) in loop_openings
            delete!(loop_openings, nameof(input))
            sys, (input_vs,) = apply_transformation(Break(input), sys)
            append!(input_vars, input_vs)
        else
            sys, (input_var,) = apply_transformation(PerturbOutput(input), sys)
            push!(input_vars, input_var)
        end
    end
    output_vars = SymbolicT[]
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
    sys, loop_opening_params = handle_loop_openings(sys, map(AnalysisPoint, collect(loop_openings)))
    return sys, input_vars, output_vars, loop_opening_params
end

function linearization_function(
        sys::AbstractSystem,
        inputs::Union{Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint}},
        outputs; loop_openings = [], system_modifier = identity, kwargs...
    )
    sys, input_vars, output_vars,
        loop_opening_params = linearization_ap_transform(
        sys, inputs, outputs, loop_openings
    )
    return linearization_function(
        system_modifier(sys), input_vars, output_vars; loop_opening_params, kwargs...
    )
end

"""
    sys, input_vars, output_vars = isolate_subsystem(sys, input_aps, output_aps)

Isolate the part of the unsimplified, hierarchical system `sys` bounded by the input
analysis points `input_aps` and the output analysis points `output_aps`. The returned
`sys` contains only the subsystems between the boundary analysis points at every level
of the hierarchy; all upstream and downstream components, and all equations involving
them, are removed.

Boundary analysis points may reside at any level of the hierarchy and in different
branches of the subsystem tree.

Returns:
- `sys`: The system with only the isolated subsystems and their internal connections.
- `input_vars`: Variables at the inside face of each input analysis point.
- `output_vars`: Variables at the inside face of each output analysis point.
"""
function isolate_subsystem(
        sys::AbstractSystem,
        input_aps::Union{Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint}},
        output_aps::Union{Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint}}
    )
    # Run clock inference first to be able to port the results to the isolated subsystem
    ts = TearingState(expand_connections(sys))
    ci = MTKTearing.ClockInference(ts)
    MTKTearing.infer_clocks!(ci)
    all_clock_subs = Dict{SymbolicT, SymbolicT}()
    if length(ci.all_clocks) > 1
        sizehint!(all_clock_subs, length(ci.var_domain))
        for (var, clk) in zip(ts.fullvars, ci.var_domain)
            # Clock operators either define a clock or depend on the clock of their inputs.
            # Neither case needs to be propagated, and this simplifies downstream
            # implementation.
            Moshi.Match.@match var begin
                BSImpl.Term(; f) && if f isa Operator end => continue
                _ => nothing
            end
            all_clock_subs[var] = setmetadata(var, VariableTimeDomain, clk)
            arr, isidx = split_indexed_var(var)
            if isidx
                get!(() -> setmetadata(arr, VariableTimeDomain, clk), all_clock_subs, var)
            end
        end
    end

    input_aps = canonicalize_ap(sys, input_aps)
    output_aps = canonicalize_ap(sys, output_aps)

    # Step 1: index every named subsystem in the hierarchy.
    # paths[i] is the path of names from sys to that subsystem, e.g. [:inner, :P].
    paths = Vector{Symbol}[]
    path_to_idx = Dict{Vector{Symbol}, Int}()
    function _collect_subsystems!(cur, path)
        for s in get_systems(cur)
            p = [path; nameof(s)]
            push!(paths, p)
            path_to_idx[p] = length(paths)
            _collect_subsystems!(s, p)
        end
        return
    end
    _collect_subsystems!(sys, Symbol[])

    g = SimpleGraph(length(paths))

    input_ap_names = Set{Symbol}(nameof(ap) for ap in input_aps)
    output_ap_names = Set{Symbol}(nameof(ap) for ap in output_aps)
    boundary_ap_names = union(input_ap_names, output_ap_names)
    sys_root_name = nameof(sys)

    # Reconstruct the full canonical AP name (as produced by canonicalize_ap) from the
    # local name seen in an equation at depth `parent_path` within `sys`.
    function _full_ap_name(parent_path, local_name)
        return Symbol(join(string.([sys_root_name; parent_path; local_name]), NAMESPACE_SEPARATOR))
    end

    # Map a connector (port System or signal SymbolicT) at a given parent path to the
    # component path by stripping the last namespace segment (port or variable name).
    function _conn_to_path(conn, parent_path)
        segs = conn isa AbstractSystem ? namespace_hierarchy(nameof(conn)) :
            namespace_hierarchy(getname(unwrap(conn)))
        length(segs) <= 1 && return nothing
        return [parent_path; segs[1:(end - 1)]]
    end

    function _try_add_edge!(p1, p2)
        i1 = get(path_to_idx, p1, nothing)
        i2 = get(path_to_idx, p2, nothing)
        (i1 === nothing || i2 === nothing || i1 == i2) && return
        return add_edge!(g, i1, i2)
    end

    # Step 2: walk every equation at every level of the hierarchy.
    # Every connection (including boundary APs) adds edges between the connected
    # components. Boundary AP edges are additionally recorded in `cut_edges` so they can
    # be removed after the full graph is built. Recording rather than skipping is
    # essential: a boundary connection may be duplicated as a plain `connect` equation,
    # which would re-add the edge. Since the graph is simple, the duplicate `add_edge!` is
    # a no-op and a single `rem_edge!` afterwards guarantees the edge stays cut.
    inside_seeds = Int[]
    input_vars = SymbolicT[]
    output_vars = SymbolicT[]
    cut_edges = Tuple{Int, Int}[]

    function _build_graph!(cur, parent_path)
        for eq in get_eqs(cur)
            lhs_val = value(eq.lhs)
            rhs_val = value(eq.rhs)
            if lhs_val isa AnalysisPoint
                ap_data = rhs_val::AnalysisPoint
                fname = _full_ap_name(parent_path, nameof(ap_data))
                in_conn = ap_data.input
                out_conns = something(ap_data.outputs, [])
                is_boundary = fname in boundary_ap_names
                if is_boundary && fname in input_ap_names
                    for c in out_conns
                        push!(input_vars, ap_var(c))
                        p = _conn_to_path(c, parent_path)
                        if p !== nothing
                            idx = get(path_to_idx, p, nothing)
                            idx !== nothing && push!(inside_seeds, idx)
                        end
                    end
                end
                if is_boundary && fname in output_ap_names
                    if in_conn !== nothing
                        push!(output_vars, ap_var(in_conn))
                        p = _conn_to_path(in_conn, parent_path)
                        if p !== nothing
                            idx = get(path_to_idx, p, nothing)
                            idx !== nothing && push!(inside_seeds, idx)
                        end
                    end
                end
                in_path = in_conn !== nothing ? _conn_to_path(in_conn, parent_path) : nothing
                in_idx = in_path === nothing ? nothing : get(path_to_idx, in_path, nothing)
                for c in out_conns
                    out_path = _conn_to_path(c, parent_path)
                    out_path === nothing && continue
                    in_path === nothing && continue
                    _try_add_edge!(in_path, out_path)
                    if is_boundary && in_idx !== nothing
                        out_idx = get(path_to_idx, out_path, nothing)
                        (out_idx === nothing || out_idx == in_idx) && continue
                        push!(cut_edges, (in_idx, out_idx))
                    end
                end
            elseif rhs_val isa Connection
                conn_list = get_systems(rhs_val)
                conn_list === nothing && continue
                cps = [
                    p for c in conn_list
                        for p in (_conn_to_path(c, parent_path),) if p !== nothing
                ]
                for i in eachindex(cps), j in (i + 1):length(cps)
                    _try_add_edge!(cps[i], cps[j])
                end
            end
        end
        for s in get_systems(cur)
            _build_graph!(s, [parent_path; nameof(s)])
        end
        return
    end
    _build_graph!(sys, Symbol[])

    # Remove the boundary-connection edges now that the full graph (including any
    # duplicate plain connections) has been built.
    for (i1, i2) in cut_edges
        rem_edge!(g, i1, i2)
    end

    # Step 3: BFS from inside seeds to find all reachable (inside) components.
    inside = falses(length(paths))
    queue = unique!(copy(inside_seeds))
    for v in queue
        inside[v] = true
    end
    qi = 1
    while qi <= length(queue)
        v = queue[qi]
        qi += 1
        for w in Graphs.neighbors(g, v)
            if !inside[w]
                inside[w] = true
                push!(queue, w)
            end
        end
    end

    # Build a prefix set so that has_inside(path) is an O(1) lookup:
    # a path is "has inside" if any inside component lives under it.
    inside_prefixes = Set{Vector{Symbol}}()
    for i in eachindex(paths)
        inside[i] || continue
        p = paths[i]
        for k in eachindex(p)
            push!(inside_prefixes, p[1:k])
        end
    end
    has_inside(path) = path in inside_prefixes

    # Step 4: rebuild the system hierarchy keeping only inside components and equations.
    # Systems that are directly inside (their path is in the `inside` set) are returned
    # as-is — their internal structure belongs to the isolated subsystem.
    # Systems that are containers (in inside_prefixes only because they wrap inside
    # components) have their equations filtered and their own variables/parameters cleared.
    function _get_next_clock_substitutions(cur_subs::Dict{SymbolicT, SymbolicT}, name::Symbol)
        new_subs = empty(cur_subs)
        sizehint!(new_subs, length(cur_subs))
        for (k, v) in cur_subs
            hierarchy = namespace_hierarchy(getname(k))
            length(hierarchy) == 1 && continue
            first(hierarchy) == name || continue
            new_name = Symbol(join(Iterators.drop(hierarchy, 1), NAMESPACE_SEPARATOR_SYMBOL))
            new_k = rename(k, new_name)
            new_v = rename(v, new_name)
            new_subs[new_k] = new_v
        end

        return new_subs
    end

    function _recursively_clock_subsystems(cur, clock_subs)
        isempty(clock_subs) && return cur
        eqs = copy(get_eqs(cur))
        clock_subber = SU.IRSubstituter{false}(get_irstructure(sys), clock_subs)
        map!(clock_subber, eqs, eqs)
        syss = copy(get_systems(cur))
        map!(s -> _recursively_clock_subsystems(s, _get_next_clock_substitutions(clock_subs, nameof(s))), syss, syss)
        return ConstructionBase.setproperties(cur; eqs = eqs, systems = syss)
    end

    function _reconstruct!(cur, parent_path, clock_subs)
        idx = get(path_to_idx, parent_path, nothing)
        if idx !== nothing && inside[idx]
            # Directly-inside component — preserve it entirely.
            return _recursively_clock_subsystems(cur, clock_subs)
        end

        # Container: filter subsystems and equations, clear own vars/params so that
        # nothing from outside the isolated region leaks into the result.
        new_systems = [
            _reconstruct!(s, [parent_path; nameof(s)], _get_next_clock_substitutions(clock_subs, nameof(s)))
                for s in get_systems(cur) if has_inside([parent_path; nameof(s)])
        ]
        new_eqs = filter(get_eqs(cur)) do eq
            lhs_val = value(eq.lhs)
            rhs_val = value(eq.rhs)
            if lhs_val isa AnalysisPoint
                ap_data = rhs_val::AnalysisPoint
                _full_ap_name(parent_path, nameof(ap_data)) in boundary_ap_names && return false
                in_conn = ap_data.input
                out_conns = something(ap_data.outputs, [])
                all_conns = in_conn === nothing ? out_conns : [in_conn; out_conns]
                all_paths = [
                    p for c in all_conns
                        for p in (_conn_to_path(c, parent_path),) if p !== nothing
                ]
                return all(has_inside, all_paths)
            elseif rhs_val isa Connection
                conn_list = get_systems(rhs_val)
                conn_list === nothing && return true
                cps = [
                    p for c in conn_list
                        for p in (_conn_to_path(c, parent_path),) if p !== nothing
                ]
                return all(has_inside, cps)
            else
                return false
            end
        end
        if !isempty(clock_subs)
            clock_subber = SU.IRSubstituter{false}(get_irstructure(sys), clock_subs)
            map!(clock_subber, new_eqs, new_eqs)
        end
        # Build a fresh container so its own variables, parameters, observed equations,
        # initial conditions and non-connection equations are stripped — only the filtered
        # connections and inside subsystems are retained.
        return System(new_eqs, get_iv(cur), SymbolicT[], SymbolicT[]; name = nameof(cur), systems = new_systems)
    end

    return _reconstruct!(sys, Symbol[], all_clock_subs), input_vars, output_vars
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
