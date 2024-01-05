import SymbolicUtils: symtype, term, hasmetadata, issym
@enum VariableType VARIABLE PARAMETER BROWNIAN
struct MTKVariableTypeCtx end

getvariabletype(x, def = VARIABLE) = getmetadata(unwrap(x), MTKVariableTypeCtx, def)

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
    elseif istree(x) && operation(x) isa Symbolic
        varT === PARAMETER || isparameter(operation(x))
    elseif istree(x) && operation(x) == (getindex)
        isparameter(arguments(x)[1])
    elseif x isa Symbolic
        varT === PARAMETER
    else
        false
    end
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

Maps the variable to a state.
"""
tovar(s::Symbolic) = setmetadata(s, MTKVariableTypeCtx, VARIABLE)
tovar(s::Num) = Num(tovar(value(s)))

"""
$(SIGNATURES)

Define one or more known parameters.
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
