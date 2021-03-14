struct VariableUnit end
struct VariableConnectType end
Symbolics.option_to_metadata_type(::Val{:unit}) = VariableUnit
Symbolics.option_to_metadata_type(::Val{:connect}) = VariableConnectType

"""
$(SIGNATURES)
Takes a list of pairs of `variables=>values` and an ordered list of variables
and creates the array of values in the correct order with default values when
applicable.
"""
function varmap_to_vars(varmap::Dict, varlist; defaults=Dict())
    varmap = merge(defaults, varmap) # prefers the `varmap`
    varmap = Dict(value(k)=>value(varmap[k]) for k in keys(varmap))
    # resolve symbolic parameter expressions
    for (p, v) in pairs(varmap)
        varmap[p] = fixpoint_sub(v, varmap)
    end
    T′ = eltype(values(varmap))
    T = Base.isconcretetype(T′) ? T′ : Base.promote_typeof(values(varmap)...)
    out = Vector{T}(undef, length(varlist))
    missingvars = setdiff(varlist, keys(varmap))
    isempty(missingvars) || throw(ArgumentError("$missingvars are missing from the variable map."))

    for (i, var) in enumerate(varlist)
        out[i] = varmap[var]
    end
    out
end

function varmap_to_vars(varmap::Union{AbstractArray,Tuple},varlist; kw...)
    if eltype(varmap) <: Pair
        out = varmap_to_vars(Dict(varmap), varlist; kw...)
        if varmap isa Tuple
            (out..., )
        else
            # Note that `varmap` might be longer than `varlist`
            construct_state(varmap, out)
        end
    else
        varmap
    end
end
varmap_to_vars(varmap::DiffEqBase.NullParameters,varlist; kw...) = varmap
varmap_to_vars(varmap::Nothing,varlist; kw...) = varmap

construct_state(x::StaticArray, y) = StaticArrays.similar_type(x, eltype(y), StaticArrays.Size(size(y)...))(y)
construct_state(x::Array, y) = y
