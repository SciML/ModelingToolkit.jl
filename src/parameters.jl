import SymbolicUtils: symtype, term, hasmetadata, issym

"""
    @enum VariableType

The type of the declared variable, used for automatic identification of
variables/parameters/brownians/etc. by the `System` constructor.
"""
@enum VariableType VARIABLE PARAMETER BROWNIAN

"""
    $TYPEDEF

The symbolic metadata key for storing the `VariableType`.
"""
struct MTKVariableTypeCtx end

getvariabletype(x, def = VARIABLE) = getmetadata(unwrap(x), MTKVariableTypeCtx, def)

"""
    $TYPEDEF

Check if the variable contains the metadata identifying it as a parameter.
"""
function isparameter(x)
    x = unwrap(x)

    if x isa Symbolic && (varT = getvariabletype(x, nothing)) !== nothing
        return varT === PARAMETER
        #TODO: Delete this branch
    elseif x isa Symbolic && Symbolics.getparent(x, false) !== false
        p = Symbolics.getparent(x)
        isparameter(p) ||
            (hasmetadata(p, Symbolics.VariableSource) &&
             getmetadata(p, Symbolics.VariableSource)[1] == :parameters)
    elseif iscall(x) && operation(x) isa Symbolic
        varT === PARAMETER || isparameter(operation(x))
    elseif iscall(x) && operation(x) == (getindex)
        isparameter(arguments(x)[1])
    elseif x isa Symbolic
        varT === PARAMETER
    else
        false
    end
end

function iscalledparameter(x)
    x = unwrap(x)
    return isparameter(getmetadata(x, CallWithParent, nothing))
end

function getcalledparameter(x)
    x = unwrap(x)
    # `parent` is a `CallWithMetadata` with the correct metadata,
    # but no namespacing. `operation(x)` has the correct namespacing,
    # but is not a `CallWithMetadata` and doesn't have any metadata.
    # This approach combines both.
    parent = getmetadata(x, CallWithParent)
    return CallWithMetadata(operation(x), metadata(parent))
end

"""
    toparam(s)

Maps the variable to a parameter.
"""
function toparam(s)
    if s isa Symbolics.Arr
        Symbolics.wrap(toparam(Symbolics.unwrap(s)))
    elseif s isa AbstractArray
        map(toparam, s)
    else
        setmetadata(s, MTKVariableTypeCtx, PARAMETER)
    end
end
toparam(s::Num) = wrap(toparam(value(s)))

"""
    tovar(s)

Maps the variable to an unknown.
"""
tovar(s::Symbolic) = setmetadata(s, MTKVariableTypeCtx, VARIABLE)
tovar(s::Union{Num, Symbolics.Arr}) = wrap(tovar(unwrap(s)))

"""
$(SIGNATURES)

Define one or more known parameters.

See also [`@independent_variables`](@ref), [`@variables`](@ref) and [`@constants`](@ref).
"""
macro parameters(xs...)
    Symbolics._parse_vars(:parameters,
        Real,
        xs,
        toparam) |> esc
end

function find_types(array)
    by = let set = Dict{Any, Int}(), counter = Ref(0)
        x -> begin
            # t = typeof(x)

            get!(set, typeof(x)) do
                # if t == Float64
                #     1
                # else
                counter[] += 1
                # end
            end
        end
    end
    return by.(array)
end

function split_parameters_by_type(ps)
    if ps === SciMLBase.NullParameters()
        return Float64[], [] #use Float64 to avoid Any type warning
    else
        by = let set = Dict{Any, Int}(), counter = Ref(0)
            x -> begin
                get!(set, typeof(x)) do
                    counter[] += 1
                end
            end
        end
        idxs = by.(ps)
        split_idxs = [Int[]]
        for (i, idx) in enumerate(idxs)
            if idx > length(split_idxs)
                push!(split_idxs, Int[])
            end
            push!(split_idxs[idx], i)
        end
        tighten_types = x -> identity.(x)
        split_ps = tighten_types.(Base.Fix1(getindex, ps).(split_idxs))

        if ps isa StaticArray
            parrs = map(x -> SArray{Tuple{size(x)...}}(x), split_ps)
            split_ps = SArray{Tuple{size(parrs)...}}(parrs)
        end
        if length(split_ps) == 1  #Tuple not needed, only 1 type
            return split_ps[1], split_idxs
        else
            return (split_ps...,), split_idxs
        end
    end
end

"""
    $(TYPEDSIGNATURES)

Change the tunable parameters of a system to a new set of tunables.

The new tunable parameters must be a subset of the current tunables as discovered by [`tunable_parameters`](@ref).
The remaining parameters will be set as constants in the system.
"""
function subset_tunables(sys, new_tunables)
    if !iscomplete(sys)
        throw(ArgumentError("System must be `complete` before changing tunables."))
    end
    if !is_split(sys)
        throw(ArgumentError("Tunable parameters can only be changed for split systems."))
    end

    cur_tunables = tunable_parameters(sys, parameters(sys))
    diff_params = setdiff(cur_tunables, new_tunables)

    if !isempty(setdiff(new_tunables, cur_tunables))
        throw(ArgumentError("""New tunable parameters must be a subset of the current tunable parameters. Found tunable parameters not in the system: $(setdiff(new_tunables, cur_tunables)).
        Note that array parameters can only be set as tunable or non-tunable, not partially tunable. They should be specified in the un-scalarized form.
        """))
    end
    cur_ps = copy(get_ps(sys))
    const_ps = toconstant.(diff_params)

    for (idx, p) in enumerate(cur_ps)
        for (d, c) in zip(diff_params, const_ps)
            if isequal(p, d)
                setindex!(cur_ps, c, idx)
                break
            end
        end
    end
    @set! sys.ps = cur_ps
    complete(sys)
end
