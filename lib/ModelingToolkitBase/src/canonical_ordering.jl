"""
    $TYPEDSIGNATURES

If `x` is a small integer (fits in `Int8`) return `Int(x)`. Otherwise, return
`typemin(Int)`. Used to build the integer index portion of [`canonical_sort_key`](@ref)
without allocating.
"""
function manual_dispatch_is_small_int(@nospecialize(x::Number))::Int
    if x isa Int
        return typemin(Int8) <= x <= typemax(Int8) ? x : typemin(Int)
    elseif x isa BigInt
        return typemin(Int8) <= x <= typemax(Int8) ? Int(x) : typemin(Int)
    elseif x isa Int32
        return typemin(Int8) <= x <= typemax(Int8) ? Int(x) : typemin(Int)
    elseif x isa Float64
        return isinteger(x) && typemin(Int8) <= x <= typemax(Int8) ? Int(x) : typemin(Int)
    elseif x isa Float32
        return isinteger(x) && typemin(Int8) <= x <= typemax(Int8) ? Int(x) : typemin(Int)
    elseif x isa BigFloat
        return isinteger(x) && typemin(Int8) <= x <= typemax(Int8) ? Int(x) : typemin(Int)
    elseif x isa Rational{Int}
        return isinteger(x) && typemin(Int8) <= x <= typemax(Int8) ? Int(x) : typemin(Int)
    elseif x isa Rational{Int32}
        return isinteger(x) && typemin(Int8) <= x <= typemax(Int8) ? Int(x) : typemin(Int)
    elseif x isa Rational{BigInt}
        return isinteger(x) && typemin(Int8) <= x <= typemax(Int8) ? Int(x) : typemin(Int)
    else
        return if isinteger(x)::Bool && (typemin(Int8) <= x <= typemax(Int8))::Bool
            Int(x)::Int
        else
            typemin(Int)
        end
    end
end

"""
    $TYPEDSIGNATURES

Total, allocation-light "canonical name" for any expression that may appear as a
variable/parameter: the name for named variables (`Sym` / call-variable / `getindex`),
the operation's name for operator/function applications, and a fixed sentinel symbol for
the remaining structural variants. Key collisions are acceptable — they only mean the
canonical tie-break falls back to the original order among the colliding entries.
"""
function canonical_name(x::SymbolicT)
    Moshi.Match.@match x begin
        BSImpl.Sym(; name) => name
        BSImpl.Term(; f, args) && if f === getindex end => canonical_name(args[1])
        BSImpl.Term(; f) => canonical_opname(f)
        BSImpl.AddMul(; variant) => variant === SU.AddMulVariant.ADD ? :+ : :*
        BSImpl.Div(;) => :/
        BSImpl.ArrayOp(;) => Symbol("#arrayop")
        BSImpl.Const(;) => Symbol("#const")
        _ => Symbol("#expr")
    end
end

function canonical_opname(@nospecialize(f))
    f isa SymbolicT && return canonical_name(f)
    f isa Function && return nameof(f)::Symbol
    return nameof(typeof(f))::Symbol
end

"""
    $TYPEDSIGNATURES

Structured, hash-independent canonical sort key for a symbolic variable/parameter: a tuple
`(name, indices, opsig)` where `name` is the [`canonical_name`](@ref) of the underlying
(array) variable, `indices` are the integer indices when `v` is a scalarized array element
(empty otherwise) and `opsig` encodes the operator chain wrapping the variable
(`Differential` → `1`; `Shift` → `2` followed by its step count; any other single-argument
operator → `3`). Comparing these tuples orders variables deterministically regardless of
hashing/`objectid`, declaration order or equation order, without stringifying symbolic
expressions. Use this (rather than `Dict`/`Set` iteration order or `hash`) wherever the
order of a collection of symbolic quantities reaches observable output.
"""
function canonical_sort_key(v::SymbolicT)
    # `opsig`/`idxs` are built as `Vector{Int}` (not growing tuples) so the key type is
    # concrete and inferrable after the loop. Vectors compare lexicographically, so the
    # tuple ordering is unchanged.
    opsig = Int[]
    x = v
    while true
        stripped = Moshi.Match.@match x begin
            BSImpl.Term(; f, args) && if f isa Differential end => begin
                push!(opsig, 1)
                args[1]
            end
            BSImpl.Term(; f, args) && if f isa Shift end => begin
                push!(opsig, 2, Int(f.steps))
                args[1]
            end
            BSImpl.Term(; f, args) && if f isa SU.Operator && length(args) == 1 end => begin
                push!(opsig, 3)
                args[1]
            end
            _ => nothing
        end
        stripped === nothing && break
        x = stripped
    end
    idxs = Int[]
    Moshi.Match.@match x begin
        BSImpl.Term(; f, args) && if f === getindex end => begin
            for i in Iterators.drop(args, 1)
                ival = SU.isconst(i) ? manual_dispatch_is_small_int(unwrap_const(i))::Int : 0
                push!(idxs, ival)
            end
            x = args[1]
        end
        _ => nothing
    end
    return (canonical_name(x), idxs, opsig)
end

"""
    $TYPEDSIGNATURES

Return a new vector containing the symbolic quantities `xs` in deterministic
[`canonical_sort_key`](@ref) order. Stable: ties retain their original relative order.
"""
canonical_sort(xs) = sort(collect(xs); by = canonical_sort_key, alg = Base.Sort.DEFAULT_STABLE)
