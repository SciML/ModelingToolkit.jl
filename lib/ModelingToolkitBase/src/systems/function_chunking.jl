"""
Number of non-trivial statements above which `eval_or_rgf` outlines the body of a generated
function into helper functions when `CompilerOptions.chunk_size` is left at its default.
"""
const DEFAULT_CHUNK_THRESHOLD = 20000

"""
Number of statements per helper function when `CompilerOptions.chunk_size` is left at its
default.
"""
const DEFAULT_CHUNK_SIZE = 500

"""
Largest ratio of values returned by helpers to statements they contain for which automatic
chunking goes ahead. Densely coupled bodies (e.g. 2D stencils sharing common subexpressions
across the whole grid) return so many values that chunking costs more than it saves.
"""
const DEFAULT_MAX_TRAFFIC = 0.15

# Heads that may appear inside the right-hand side of a statement that is safe to move into
# a helper function: none of them introduce bindings.
const _CHUNKABLE_HEADS = (:call, :ref, :., :curly)

_is_inbounds_call(ex) = ex isa Expr && ex.head === :macrocall && ex.args[1] === Symbol("@inbounds")

function _is_binding_free(ex)
    ex isa Expr || return true
    ex.head in _CHUNKABLE_HEADS || return false
    return all(_is_binding_free, ex.args)
end

# Every symbol mentioned in `ex`, including ones that are only names of globals or functions.
function _collect_symbols!(out::Vector{Symbol}, ex)
    if ex isa Symbol
        push!(out, ex)
    elseif ex isa Expr
        foreach(Base.Fix1(_collect_symbols!, out), ex.args)
    end
    return out
end

# Every symbol bound anywhere inside `ex`. Over-approximates, which is what is needed for
# statements that stay in the calling function.
function _collect_bound!(out::Vector{Symbol}, ex)
    ex isa Expr || return out
    if ex.head in (:(=), :local, :global, :(.=)) || ex.head === :(::)
        _collect_lhs!(out, ex.head === :(=) ? ex.args[1] : ex)
    end
    foreach(Base.Fix1(_collect_bound!, out), ex.args)
    return out
end
function _collect_lhs!(out, lhs)
    if lhs isa Symbol
        push!(out, lhs)
    elseif lhs isa Expr
        foreach(Base.Fix1(_collect_lhs!, out), lhs.args)
    end
    return out
end

# Symbols read by `ex`: like `_collect_symbols!`, but the target of a plain assignment is not
# a read.
function _collect_reads!(out::Vector{Symbol}, ex)
    if ex isa Symbol
        push!(out, ex)
    elseif ex isa Expr
        if ex.head === :(=) && ex.args[1] isa Symbol
            _collect_reads!(out, ex.args[2])
        else
            foreach(Base.Fix1(_collect_reads!, out), ex.args)
        end
    end
    return out
end

# Descend through single-statement blocks and `@inbounds` to the first block that holds
# several statements. Returns the statements, a function rebuilding the original nesting
# around a replacement list of statements, and whether an `@inbounds` was passed on the way
# (so that the statements run without bounds checks), or `nothing`.
function _find_statement_block(ex, inbounds = false)
    if ex isa Expr && ex.head === :block
        kids = filter(!Base.Fix2(isa, LineNumberNode), ex.args)
        length(kids) == 1 || return (kids, stmts -> Expr(:block, stmts...), inbounds)
        found = _find_statement_block(only(kids), inbounds)
        found === nothing && return nothing
        stmts, rebuild, inb = found
        return stmts, s -> Expr(:block, rebuild(s)), inb
    elseif _is_inbounds_call(ex) && ex.args[end] isa Expr && ex.args[end].head === :block
        found = _find_statement_block(ex.args[end], true)
        found === nothing && return nothing
        stmts, rebuild, inb = found
        return stmts, s -> Expr(ex.head, ex.args[1:(end - 1)]..., rebuild(s)), inb
    end
    return nothing
end

# Statements of a block, with `@inbounds` and nested `begin` blocks spliced in.
function _splice_block!(out, ex)
    if ex isa Expr && ex.head === :block
        for kid in ex.args
            kid isa LineNumberNode || _splice_block!(out, kid)
        end
    elseif _is_inbounds_call(ex) && ex.args[end] isa Expr && ex.args[end].head === :block
        _splice_block!(out, ex.args[end])
    else
        push!(out, ex)
    end
    return out
end

# The final `x = @inbounds begin ...; y = z end` written by array codegen is spliced into
# the statement list so that the stores can live next to the values they store.
function _splice_tail(stmts)
    isempty(stmts) && return stmts
    last_stmt = stmts[end]
    inner = last_stmt
    target = nothing
    if last_stmt isa Expr && last_stmt.head === :(=) && last_stmt.args[1] isa Symbol
        target = last_stmt.args[1]
        inner = last_stmt.args[2]
    end
    (inner isa Expr && (inner.head === :block || _is_inbounds_call(inner))) || return stmts
    spliced = _splice_block!(Any[], inner)
    isempty(spliced) && return stmts
    final = spliced[end]
    if target !== nothing
        final isa Expr && final.head === :(=) && final.args[1] isa Symbol || return stmts
        push!(spliced, Expr(:(=), target, final.args[1]))
    end
    return [stmts[1:(end - 1)]; spliced]
end

struct _ChunkStatement
    expr::Any
    defs::Vector{Symbol}
    uses::Vector{Symbol}
    movable::Bool
end

function _analyze_statement(st)
    inner = st
    if st isa Expr && st.head === :local && length(st.args) == 1 &&
            st.args[1] isa Expr && st.args[1].head === :(=)
        inner = st.args[1]
    end
    if inner isa Expr && inner.head === :(=) && _is_binding_free(inner.args[2])
        lhs = inner.args[1]
        if lhs isa Symbol
            return _ChunkStatement(st, [lhs], _collect_symbols!(Symbol[], inner.args[2]), true)
        elseif lhs isa Expr && lhs.head === :ref && _is_binding_free(lhs)
            return _ChunkStatement(st, Symbol[], _collect_symbols!(Symbol[], inner), true)
        end
    elseif inner isa Expr && inner.head === :call && _is_binding_free(inner)
        return _ChunkStatement(st, Symbol[], _collect_symbols!(Symbol[], inner), true)
    end
    return _ChunkStatement(st, _collect_bound!(Symbol[], st), _collect_symbols!(Symbol[], st), false)
end

# Right-hand sides that are duplicated into helpers instead of passed to them: literals and
# plain renamings of function arguments. Values read from arguments are never recomputed, so
# every read happens exactly where and when the original function performed it.
function _is_cheap_rhs(rhs, cheap::Set{Symbol}, fargs::Set{Symbol})
    rhs isa Number && return true
    rhs isa Symbol && return rhs in fargs || rhs in cheap
    return false
end

# Forms whose scoping or side effects the pass does not model: closures, inner functions,
# comprehensions and macros other than `@inbounds`.
function _has_unsafe_form(ex)
    ex isa Expr || return false
    ex.head in (:->, :function, :do, :comprehension, :generator, :opaque_closure, :macro) && return true
    ex.head === :macrocall && !_is_inbounds_call(ex) && return true
    ex.head === :(=) && ex.args[1] isa Expr && ex.args[1].head === :call && return true
    return any(_has_unsafe_form, ex.args)
end

# Functions that real generated code calls and that neither mutate their arguments nor
# anything else. Only these may be called by a chunked function.
const _PURE_NAMES = (
    :+, :-, :*, :/, :^, :\, :inv, :abs, :abs2, :sqrt, :cbrt, :exp, :exp2, :exp10, :expm1,
    :log, :log2, :log10, :log1p, :sin, :cos, :tan, :asin, :acos, :atan, :sinh, :cosh,
    :tanh, :asinh, :acosh, :atanh, :sec, :csc, :cot, :sign, :signbit, :floor, :ceil,
    :round, :trunc, :mod, :rem, :div, :fld, :cld, :min, :max, :ifelse, :hypot, :copysign,
    :flipsign, :muladd, :fma, :literal_pow, :Val, :typeof, :>=, :<=, :>, :<, :(==), :!=,
    :isless, :isequal, :isnan, :isinf, :isfinite, :!,
)

# The same functions as values, as they appear embedded in generated expressions, plus the
# array constructors of SymbolicUtils and the sparse-Jacobian store.
const _PURE_VALUES = let vals = Any[Val, typeof]
    for name in _PURE_NAMES
        name in (:Val, :typeof) && continue
        isdefined(Base, name) && push!(vals, getfield(Base, name))
    end
    push!(vals, Base.copy)
    for name in (:create_array, :array_literal)
        isdefined(SU.Code, name) && push!(vals, getfield(SU.Code, name))
        isdefined(SU, name) && push!(vals, getfield(SU, name))
    end
    vals
end

_is_sparse_store(f) = f isa Function && nameof(f) === :set_nzval && parentmodule(f) === Symbolics

# How a call is classified: `:pure`, `:store` (writes into an array) or `:bad`. A callee
# written as a name is only trusted if no parameter or local of the function has that name.
function _callee_kind(f, bound)
    if f isa Symbol
        return f in bound || !(f in _PURE_NAMES) ? :bad : :pure
    elseif f isa Expr && f.head === :. && length(f.args) == 2
        mod, name = f.args
        ok = mod in (:Base, :NaNMath) && !(mod in bound) && name isa QuoteNode &&
            name.value in _PURE_NAMES
        return ok ? :pure : :bad
    elseif f isa Expr && f.head === :curly && f.args[1] === :Val && !(:Val in bound)
        return all(a -> !(a isa Expr || a isa Symbol), f.args[2:end]) ? :pure : :bad
    elseif f isa GlobalRef
        ok = nameof(f.mod) in (:Base, :Core, :NaNMath) && f.name in _PURE_NAMES
        return ok ? :pure : :bad
    elseif _is_sparse_store(f)
        return :store
    end
    return any(v -> v === f, _PURE_VALUES) ? :pure : :bad
end

"""
Strict whitelist scan of the body: only the expression heads that occur in generated
functions (`=`, `local`, `block`, `ref`, `call`, broadcast calls and `@inbounds`), calls to the
functions above or `Symbolics.set_nzval`, and no rebinding of arguments. Anything else, such
as conditionals, tuples, type assertions or broadcast assignment, makes the function
unchunked.
"""
function _is_whitelisted(stmts, fargs, bound)
    # `top` is true while the expression is evaluated unconditionally as part of a top-level
    # statement; only there may a name be bound, so that every name a helper reads has
    # definitely been assigned.
    function visit(ex, top)
        ex isa Expr || return true
        h = ex.head
        if h === :ref
            return all(a -> visit(a, false), ex.args)
        elseif h === :block
            return all(a -> visit(a, top), ex.args)
        elseif h === :(=)
            lhs, rhs = ex.args
            if lhs isa Symbol
                return top && !(lhs in fargs) && visit(rhs, top)
            end
            return lhs isa Expr && lhs.head === :ref && visit(lhs, false) && visit(rhs, false)
        elseif h === :local
            # Only `local name = expr`: a bare declaration would not initialize the name.
            return top && length(ex.args) == 1 && ex.args[1] isa Expr &&
                ex.args[1].head === :(=) && ex.args[1].args[1] isa Symbol && visit(ex.args[1], top)
        elseif h === :call
            _callee_kind(ex.args[1], bound) === :bad && return false
            return all(a -> visit(a, false), @view ex.args[2:end])
        elseif h === :. && length(ex.args) == 2 && ex.args[2] isa Expr && ex.args[2].head === :tuple
            # broadcast call `f.(args...)`
            _callee_kind(ex.args[1], bound) === :bad && return false
            return all(a -> visit(a, false), ex.args[2].args)
        elseif _is_inbounds_call(ex)
            return visit(ex.args[end], top)
        end
        return false
    end
    return all(st -> visit(st, true), stmts)
end

"""
    chunk_function_expr(fex::Expr, chunk_size::Int, threshold::Int)

Split the straight-line body of the generated function expression `fex` into helper
functions of about `chunk_size` statements, returning the rewritten function expression.
The helpers are local anonymous functions that capture nothing, so they are lowered and
inferred separately from the calling function. Values that are cheap to recompute (literals and constant-index reads
of arguments) are recomputed in each helper instead of being passed to it. Returns `fex`
unchanged if it has fewer than `threshold` statements, an unrecognized shape, or if the
helpers would return more than `max_traffic` values per statement they contain.
"""
function chunk_function_expr(fex::Expr, chunk_size::Int, threshold::Int, max_traffic::Real = Inf)
    # Chunking is an optimization: any shape the analysis trips over leaves `fex` as is.
    try
        return _chunk_function_expr(fex, chunk_size, threshold, max_traffic)
    catch
        return fex
    end
end

function _chunk_function_expr(fex::Expr, chunk_size::Int, threshold::Int, max_traffic::Real)
    fex.head === :function && length(fex.args) == 2 || return fex
    sig = fex.args[1]
    sig isa Expr && sig.head === :tuple && all(a -> a isa Symbol, sig.args) || return fex
    fargs = Set{Symbol}(sig.args)
    found = _find_statement_block(fex.args[2])
    found === nothing && return fex
    stmts, rebuild, inbounds = found
    isempty(stmts) && return fex
    # The stores of array codegen sit in their own `@inbounds` block; only fold them into
    # the statement list when the whole body already runs inside `@inbounds`.
    inbounds && (stmts = _splice_tail(stmts))
    any(_has_unsafe_form, stmts) && return fex
    infos = map(_analyze_statement, stmts)
    n_defs = Dict{Symbol, Int}()
    for info in infos, d in info.defs
        n_defs[d] = get(n_defs, d, 0) + 1
    end
    # Rebinding an argument would need renaming across helpers.
    any(in(fargs), keys(n_defs)) && return fex
    _is_whitelisted(stmts, fargs, union(fargs, keys(n_defs))) || return fex
    # Every local that is read must have been bound by an earlier statement.
    first_def = Dict{Symbol, Int}()
    for (i, info) in enumerate(infos), d in info.defs
        get!(first_def, d, i)
    end
    for (i, st) in enumerate(stmts), r in _collect_reads!(Symbol[], st)
        haskey(first_def, r) && first_def[r] >= i && return fex
    end

    # A cheap statement stays in the calling function and is also duplicated into any helper
    # that needs it.
    cheap = Set{Symbol}()
    is_cheap = falses(length(infos))
    for (i, info) in enumerate(infos)
        info.movable && length(info.defs) == 1 || continue
        d = only(info.defs)
        n_defs[d] == 1 || continue
        inner = info.expr
        inner = inner.head === :local ? inner.args[1] : inner
        _is_cheap_rhs(inner.args[2], cheap, fargs) || continue
        push!(cheap, d)
        is_cheap[i] = true
    end
    # The last statement provides the return value and stays in the calling function.
    outlinable = [infos[i].movable && !is_cheap[i] && i < length(infos) for i in eachindex(infos)]
    count(outlinable) < threshold && return fex

    last_use = Dict{Symbol, Int}()
    for (i, info) in enumerate(infos), u in info.uses
        last_use[u] = i
    end
    cheap_stmt = Dict{Symbol, Int}(only(infos[i].defs) => i for i in eachindex(infos) if is_cheap[i])

    # Pieces are contiguous runs of outlinable statements, cut to about `chunk_size`.
    new_stmts = Any[]
    traffic = outlined = 0
    # Helper names must not collide with any name already used in the function.
    used = _collect_symbols!(Symbol[], fex)
    prefix = "__mtk_chunk_"
    while any(u -> startswith(string(u), prefix), used)
        prefix = "_" * prefix
    end
    i = 1
    n = length(infos)
    while i <= n
        if !outlinable[i]
            push!(new_stmts, stmts[i])
            i += 1
            continue
        end
        # Cheap statements interleaved with outlinable ones do not end a run.
        j = i
        k = i
        while k <= n && (outlinable[k] || is_cheap[k])
            outlinable[k] && (j = k)
            k += 1
        end
        run = [m for m in i:j if outlinable[m]]
        npieces = max(1, cld(length(run), chunk_size))
        piece_len = cld(length(run), npieces)
        # Cheap statements of the run stay in the caller, placed before the first call
        # that could need them; helpers recompute what they use.
        for m in i:j
            is_cheap[m] && push!(new_stmts, stmts[m])
        end
        for piece in Iterators.partition(run, piece_len)
            helper, call, n_out = _outline_piece(
                prefix, length(new_stmts), inbounds, infos, stmts, piece, cheap, cheap_stmt, fargs,
                last_use, first_def
            )
            push!(new_stmts, helper, call)
            traffic += n_out
            outlined += length(piece)
        end
        i = j + 1
    end
    traffic > max_traffic * outlined && return fex
    return Expr(:function, fex.args[1], rebuild(new_stmts))
end

function _outline_piece(prefix, id, inbounds, infos, stmts, piece, cheap, cheap_stmt, fargs, last_use, first_def)
    first_idx, last_idx = first(piece), last(piece)
    defined = Set{Symbol}()
    inputs = Symbol[]
    needed_cheap = Symbol[]
    function note_use!(u)
        u in defined && return
        if u in cheap
            u in needed_cheap || push!(needed_cheap, u)
        elseif u in fargs
            u in inputs || push!(inputs, u)
        elseif get(first_def, u, typemax(Int)) < first_idx
            u in inputs || push!(inputs, u)
        end
        return
    end
    for m in piece
        foreach(note_use!, infos[m].uses)
        union!(defined, infos[m].defs)
    end
    # Close `needed_cheap` over its own dependencies; they only depend on arguments and
    # other cheap values.
    queue = copy(needed_cheap)
    while !isempty(queue)
        c = pop!(queue)
        for u in infos[cheap_stmt[c]].uses
            if u in cheap
                u in needed_cheap || (push!(needed_cheap, u); push!(queue, u))
            elseif u in fargs
                u in inputs || push!(inputs, u)
            end
        end
    end
    # A helper cannot both take a variable as input and rebind it (`a = a + 1`); the caller
    # of this function treats the error as "leave the function unchunked".
    any(in(defined), inputs) && error("variable rebound from an earlier value")
    sort!(needed_cheap; by = c -> cheap_stmt[c])
    outputs = Symbol[]
    for m in piece, d in infos[m].defs
        get(last_use, d, 0) > last_idx && !(d in outputs) && push!(outputs, d)
    end

    body = Any[_as_local(stmts[cheap_stmt[c]]) for c in needed_cheap]
    append!(body, _as_local(stmts[m]) for m in piece)
    ret = length(outputs) == 1 ? only(outputs) : Expr(:tuple, outputs...)
    push!(body, ret)
    helper_body = Expr(:block, body...)
    inbounds && (helper_body = Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(0), helper_body))
    helper = Expr(:function, Expr(:tuple, inputs...), helper_body)
    name = Symbol(prefix, id)
    call = Expr(:call, name, inputs...)
    isempty(outputs) && return (Expr(:(=), name, helper), call, 0)
    lhs = length(outputs) == 1 ? only(outputs) : Expr(:tuple, outputs...)
    return (Expr(:(=), name, helper), Expr(:(=), lhs, call), length(outputs))
end

# A helper must bind its own variables: assigning a name that the calling function also
# binds would make the helper capture it, boxing the variable.
_as_local(st) = st isa Expr && st.head === :(=) && st.args[1] isa Symbol ? Expr(:local, st) : st
