import SymbolicUtils: symtype, term, hasmetadata, issym

"""
    @enum VariableType

The type of the declared variable, used for automatic identification of
variables/parameters/brownians/etc. by the `System` constructor.
"""
@enum VariableType VARIABLE PARAMETER BROWNIAN POISSONIAN

"""
    $TYPEDEF

The symbolic metadata key for storing the `VariableType`.
"""
struct MTKVariableTypeCtx end

getvariabletype(x, def = VARIABLE) = safe_getmetadata(MTKVariableTypeCtx, unwrap(x), def)::Union{typeof(def), VariableType}

"""
    $TYPEDEF

The symbolic metadata key for storing the rate expression of a poissonian variable.
"""
struct PoissonianRateCtx end

"""
    getpoissonianrate(x)

Get the rate expression associated with poissonian variable `x`. Returns `nothing`
if `x` is not a poissonian or has no rate.
"""
getpoissonianrate(x::Union{Num, Symbolics.Arr}) = getpoissonianrate(Symbolics.unwrap(x))
function getpoissonianrate(x)
    return Symbolics.getmetadata(unwrap(x), PoissonianRateCtx, nothing)
end

"""
    $TYPEDEF

Check if the variable contains the metadata identifying it as a parameter.
"""
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

"""
    toparam(s)

Maps the variable to a parameter.
"""
function toparam(s)
    return if s isa Symbolics.Arr
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
tovar(s::SymbolicT) = setmetadata(s, MTKVariableTypeCtx, VARIABLE)
tovar(s::Union{Num, Symbolics.Arr}) = wrap(tovar(unwrap(s)))

function toparam_validate(s::SymbolicT)
    if iscall(s)
        error(
            lazy"""
            `@parameters` cannot create time-dependent parameters. Encountered $s. Use \
            `@discretes` for this purpose.
            """
        )
    end
    return toparam(s)
end
function toparam_validate(s::Union{Num, Symbolics.Arr, Symbolics.CallAndWrap})
    return typeof(s)(toparam_validate(unwrap(s)))
end

"""
$(SIGNATURES)

Define one or more known parameters. A parameter is a non-time-dependent quantity
in the model.

See also [`@independent_variables`](@ref),
[`@variables`](https://docs.sciml.ai/Symbolics/stable/manual/variables/#Symbolics.@variables)
and [`@constants`](@ref).
"""
macro parameters(xs...)
    return Symbolics.parse_vars(
        :parameters,
        Real,
        xs,
        toparam_validate
    )
end

"""
    $(TYPEDSIGNATURES)

Create a parameter `name` of the (non-`Real`) type `T` that is declared to have the array
shape `dims`, one `Int` or `AbstractUnitRange` per dimension. `@parameters K::T` leaves the
shape unknown, so `K * u` takes the shape of `u`, and `@parameters K[1:m, 1:n]::T` builds an
array of `T`. This is the way to declare a rectangular opaque operator, for example a
`SparseMatrixCSC` or a user-defined struct used through `@register_array_symbolic`:

```julia
Kop = shaped_parameter(:Kop, FEMOp, 4, 6)   # maps R^6 -> R^4
@register_array_symbolic apply_op(op::FEMOp, x::AbstractVector) begin
    size = (size(op, 1),); eltype = Float64; ndims = 1
end
```

`T` must be a concrete type. The parameter is not tunable unless `tunable = true`.
"""
function shaped_parameter(
        name::Symbol, T::Type, dims::Union{Int, AbstractUnitRange{Int}}...;
        tunable::Bool = false
    )
    isconcretetype(T) || throw(ArgumentError("`shaped_parameter` requires a concrete type, got `$T`."))
    isempty(dims) && throw(ArgumentError("`shaped_parameter` requires at least one dimension."))
    shape = SymbolicUtils.ShapeVecT(map(d -> d isa Int ? (1:d) : UnitRange{Int}(d), dims))
    sym = SymbolicUtils.Sym{SymReal}(name; type = T, shape)
    sym = setmetadata(sym, VariableTunable, tunable)
    return wrap(toparam(sym))
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
        throw(
            ArgumentError(
                """New tunable parameters must be a subset of the current tunable parameters. Found tunable parameters not in the system: $(setdiff(new_tunables, cur_tunables)).
                Note that array parameters can only be set as tunable or non-tunable, not partially tunable. They should be specified in the un-scalarized form.
                """
            )
        )
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
    return complete(sys)
end
