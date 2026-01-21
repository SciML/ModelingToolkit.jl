"""
    $(TYPEDSIGNATURES)

Validate the rules for replacement of subcomponents as defined in `substitute_component`.
"""
function validate_replacement_rule(
        rule::Pair{T, T}; namespace = []
    ) where {T <: AbstractSystem}
    lhs, rhs = rule

    iscomplete(lhs) && throw(ArgumentError("LHS of replacement rule cannot be completed."))
    iscomplete(rhs) && throw(ArgumentError("RHS of replacement rule cannot be completed."))

    rhs_h = namespace_hierarchy(nameof(rhs))
    if length(rhs_h) != 1
        throw(ArgumentError("RHS of replacement rule must not be namespaced."))
    end
    rhs_h[1] == namespace_hierarchy(nameof(lhs))[end] ||
        throw(ArgumentError("LHS and RHS must have the same name."))

    if !isequal(get_iv(lhs), get_iv(rhs))
        throw(ArgumentError("LHS and RHS of replacement rule must have the same independent variable."))
    end

    lhs_u = get_unknowns(lhs)
    rhs_u = Dict(get_unknowns(rhs) .=> nothing)
    for u in lhs_u
        if !haskey(rhs_u, u)
            if isempty(namespace)
                throw(ArgumentError("RHS of replacement rule does not contain unknown $u."))
            else
                throw(ArgumentError("Subsystem $(join([namespace; nameof(lhs)], NAMESPACE_SEPARATOR)) of RHS does not contain unknown $u."))
            end
        end
        ru = getkey(rhs_u, u, nothing)
        name = join(
            [namespace; nameof(lhs); (hasname(u) ? getname(u) : Symbol(u))],
            NAMESPACE_SEPARATOR
        )
        l_connect = something(getconnect(u), Equality)
        r_connect = something(getconnect(ru), Equality)
        if l_connect != r_connect
            throw(ArgumentError("Variable $(name) should have connection metadata $(l_connect),"))
        end

        l_input = isinput(u)
        r_input = isinput(ru)
        if l_input != r_input
            throw(ArgumentError("Variable $name has differing causality. Marked as `input = $l_input` in LHS and `input = $r_input` in RHS."))
        end
        l_output = isoutput(u)
        r_output = isoutput(ru)
        if l_output != r_output
            throw(ArgumentError("Variable $name has differing causality. Marked as `output = $l_output` in LHS and `output = $r_output` in RHS."))
        end
    end

    lhs_p = get_ps(lhs)
    rhs_p = Set(get_ps(rhs))
    for p in lhs_p
        if !(p in rhs_p)
            if isempty(namespace)
                throw(ArgumentError("RHS of replacement rule does not contain parameter $p"))
            else
                throw(ArgumentError("Subsystem $(join([namespace; nameof(lhs)], NAMESPACE_SEPARATOR)) of RHS does not contain parameter $p."))
            end
        end
    end

    lhs_s = get_systems(lhs)
    rhs_s = Dict(nameof(s) => s for s in get_systems(rhs))

    for s in lhs_s
        if haskey(rhs_s, nameof(s))
            rs = rhs_s[nameof(s)]
            if isconnector(s)
                name = join([namespace; nameof(lhs); nameof(s)], NAMESPACE_SEPARATOR)
                if !isconnector(rs)
                    throw(ArgumentError("Subsystem $name of RHS is not a connector."))
                end
                if (lct = get_connector_type(s)) !== (rct = get_connector_type(rs))
                    throw(ArgumentError("Subsystem $name of RHS has connection type $rct but LHS has $lct."))
                end
            end
            validate_replacement_rule(s => rs; namespace = [namespace; nameof(rhs)])
            continue
        end
        name1 = join([namespace; nameof(lhs)], NAMESPACE_SEPARATOR)
        throw(ArgumentError("$name1 of replacement rule does not contain subsystem $(nameof(s))."))
    end
    return
end

"""
    $(TYPEDSIGNATURES)

Chain `getproperty` calls on `root` in the order given in `hierarchy`.

# Keyword Arguments

- `skip_namespace_first`: Whether to avoid namespacing in the first `getproperty` call.
"""
function recursive_getproperty(
        root::AbstractSystem, hierarchy::Vector{Symbol}; skip_namespace_first = true
    )
    cur = root
    for (i, name) in enumerate(hierarchy)
        cur = getproperty(cur, name; namespace = i > 1 || !skip_namespace_first)
    end
    return unwrap(cur)
end

"""
    $(TYPEDSIGNATURES)

Recursively descend through `sys`, finding all connection equations and re-creating them
using the names of the involved variables/systems and finding the required variables/
systems in the hierarchy.
"""
function recreate_connections(sys::AbstractSystem)
    eqs = map(get_eqs(sys)) do eq
        eq.lhs isa Union{Connection, AnalysisPoint} || return eq
        if eq.lhs isa Connection
            oldargs = get_systems(eq.rhs)
        else
            ap::AnalysisPoint = eq.rhs
            oldargs = [ap.input; ap.outputs]
        end
        newargs = map(get_systems(eq.rhs)::Union{Vector{System}, Vector{SymbolicT}}) do arg
            name = arg isa AbstractSystem ? nameof(arg) : getname(arg)
            hierarchy = namespace_hierarchy(name)
            newarg = recursive_getproperty(sys, hierarchy)
            return newarg
        end
        if eq.lhs isa Connection
            return eq.lhs ~ Connection(newargs)
        else
            return eq.lhs ~ AnalysisPoint(newargs[1], eq.rhs.name, newargs[2:end])
        end
    end
    @set! sys.eqs = eqs
    @set! sys.systems = map(recreate_connections, get_systems(sys))
    return sys
end

"""
    $(TYPEDSIGNATURES)

Given a hierarchical system `sys` and a rule `lhs => rhs`, replace the subsystem `lhs` in
`sys` by `rhs`. The `lhs` must be the namespaced version of a subsystem of `sys` (e.g.
obtained via `sys.inner.component`). The `rhs` must be valid as per the following
conditions:

1. `rhs` must not be namespaced.
2. The name of `rhs` must be the same as the unnamespaced name of `lhs`.
3. Neither one of `lhs` or `rhs` can be marked as complete.
4. Both `lhs` and `rhs` must share the same independent variable.
5. `rhs` must contain at least all of the unknowns and parameters present in
   `lhs`.
6. Corresponding unknowns in `rhs` must share the same connection and causality
   (input/output) metadata as their counterparts in `lhs`.
7. For each subsystem of `lhs`, there must be an identically named subsystem of `rhs`.
   These two corresponding subsystems must satisfy conditions 3, 4, 5, 6, 7. If the
   subsystem of `lhs` is a connector, the corresponding subsystem of `rhs` must also
   be a connector of the same type.

`sys` also cannot be marked as complete.
"""
function substitute_component(sys::T, rule::Pair{T, T}) where {T <: AbstractSystem}
    iscomplete(sys) &&
        throw(ArgumentError("Cannot replace subsystems of completed systems"))

    validate_replacement_rule(rule)

    lhs, rhs = rule
    hierarchy = namespace_hierarchy(nameof(lhs))

    newsys, _ = modify_nested_subsystem(sys, hierarchy) do inner
        return rhs, ()
    end
    return recreate_connections(newsys)
end
