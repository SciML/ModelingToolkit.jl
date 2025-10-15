import SymbolicUtils: symtype, term, hasmetadata, issym
""""""
@enum VariableType VARIABLE PARAMETER BROWNIAN
""""""
struct MTKVariableTypeCtx end
getvariabletype(x, def = VARIABLE) = safe_getmetadata(MTKVariableTypeCtx, unwrap(x), def)::Union{typeof(def), VariableType}
""""""
isparameter(x::Union{Num, Symbolics.Arr, Symbolics.CallAndWrap}) = isparameter(unwrap(x))
function isparameter(x::SymbolicT)
    varT = getvariabletype(x, nothing)
    return varT === PARAMETER
end
isparameter(x) = false
function iscalledparameter(x)
    x = unwrap(x)
    return SymbolicUtils.is_called_function_symbolic(x) && isparameter(operation(x))
end
function getcalledparameter(x)
    x = unwrap(x)
    @assert iscalledparameter(x)
    return operation(x)
end
""""""
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
""""""
tovar(s::SymbolicT) = setmetadata(s, MTKVariableTypeCtx, VARIABLE)
tovar(s::Union{Num, Symbolics.Arr}) = wrap(tovar(unwrap(s)))
""""""
macro parameters(xs...)
    Symbolics.parse_vars(:parameters,
        Real,
        xs,
        toparam)
end
function find_types(array)
    by = let set = Dict{Any, Int}(), counter = Ref(0)
        x -> begin
            get!(set, typeof(x)) do
                counter[] += 1
            end
        end
    end
    return by.(array)
end
function split_parameters_by_type(ps)
    if ps === SciMLBase.NullParameters()
        return Float64[], []
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
        if length(split_ps) == 1
            return split_ps[1], split_idxs
        else
            return (split_ps...,), split_idxs
        end
    end
end
""""""
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
