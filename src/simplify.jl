function simplify_constants(t::Term, shorten_tree = true)
    while true
        t′ = _simplify_constants(t, shorten_tree)
        if is_branch(t′)
            t′ = map(x -> simplify_constants(x, shorten_tree), t′)
        end
        t == t′ && return t
        t = t′
    end
end

const AC_OPERATORS = (*, +)

function _simplify_constants(t::Term, shorten_tree = true)
    is_branch(t) || return t

    head = root(t)
    head === :call || return t

    fn, args = unpack(t)

    # Tree shrinking
    if shorten_tree && fn ∈ AC_OPERATORS
        # Flatten tree
        idxs = findall(x -> is_branch(x) && unpack(x)[1] === fn, args)
        if !isempty(idxs)
            keep_idxs = eachindex(args) .∉ (idxs,)
            new_args = Vector{Term}[unpack(args[i])[2] for i in idxs]
            push!(new_args, args[keep_idxs])
            return pack(fn, vcat(new_args...))
        end

        # Collapse constants
        idxs = findall(is_constant, args)
        if length(idxs) > 1
            other_idxs = eachindex(args) .∉ (idxs,)
            new_const = mapreduce(root, fn, args[idxs])
            new_args = push!(args[other_idxs], new_const)

            length(new_args) == 1 && return first(new_args)
            return pack(fn, new_args)
        end
    end

    if fn === (*)
        # If any variable is `0`, zero the whole thing
        any(_iszero, args) && return @term(0)

        # If any variable is `1`, remove that `1` unless they are all `1`,
        # in which case simplify to a single variable
        if any(_isone, args)
            args = filter(!_isone, args)

            isempty(args)     && return @term(1)
            length(args) == 1 && return first(args)
            return pack(fn, args)
        end

        return t
    end

    if fn === (+) && any(_iszero, args)
        # If there are `0`s in a big `+` expression, get rid of them
        args = filter(!_iszero, args)

        isempty(args)     && return @term(0)
        length(args) == 1 && return first(args)
        return pack(fn, args)
    end

    (fn, length(args)) === (identity, 1) && return args[1]

    (fn, length(args)) === (-, 1) && return @term(-1 * $(args[1]))

    return t
end

_iszero(t::Term) = !is_branch(t) && iszero(root(t))
_isone(t::Term) = !is_branch(t) && isone(root(t))

export simplify_constants
