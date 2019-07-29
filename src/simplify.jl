export simplify_constants


function simplify_constants(O::Operation, shorten_tree)
    while true
        O′ = _simplify_constants(O, shorten_tree)
        if is_operation(O′)
            O′ = Operation(O′.op, simplify_constants.(O′.args, shorten_tree))
        end
        isequal(O, O′) && return O
        O = O′
    end
end
simplify_constants(x, shorten_tree) = x
simplify_constants(x) = simplify_constants(x, true)


const AC_OPERATORS = (*, +)

function _simplify_constants(O::Operation, shorten_tree)
    # Tree shrinking
    if shorten_tree && O.op ∈ AC_OPERATORS
        # Flatten tree
        idxs = findall(x -> is_operation(x) && x.op === O.op, O.args)
        if !isempty(idxs)
            keep_idxs = eachindex(O.args) .∉ (idxs,)
            args = Vector{Expression}[O.args[i].args for i in idxs]
            push!(args, O.args[keep_idxs])
            return Operation(O.op, vcat(args...))
        end

        # Collapse constants
        idxs = findall(is_constant, O.args)
        if length(idxs) > 1
            other_idxs = eachindex(O.args) .∉ (idxs,)
            new_const = Constant(mapreduce(get, O.op, O.args[idxs]))
            args = push!(O.args[other_idxs], new_const)

            length(args) == 1 && return first(args)
            return Operation(O.op, args)
        end
    end

    if O.op === (*)
        # If any variable is `Constant(0)`, zero the whole thing
        any(iszero, O.args) && return Constant(0)

        # If any variable is `Constant(1)`, remove that `Constant(1)` unless
        # they are all `Constant(1)`, in which case simplify to a single variable
        if any(isone, O.args)
            args = filter(!isone, O.args)

            isempty(args)     && return Constant(1)
            length(args) == 1 && return first(args)
            return Operation(O.op, args)
        end

        return O
    end

    if O.op === (^) && length(O.args) == 2 && iszero(O.args[2])
        return Constant(1)
    end

    if O.op === (^) && length(O.args) == 2 && isone(O.args[2])
        return O.args[1]
    end

    if O.op === (+) && any(iszero, O.args)
        # If there are Constant(0)s in a big `+` expression, get rid of them
        args = filter(!iszero, O.args)

        isempty(args)     && return Constant(0)
        length(args) == 1 && return first(args)
        return Operation(O.op, args)
    end

    if (O.op === (-) || O.op === (+) || O.op === (*)) && all(is_constant, O.args) && !isempty(O.args)
        v = O.args[1].value
        for i in 2:length(O.args)
            v = O.op(v, O.args[i].value)
        end
        return Constant(v)
    end

    (O.op, length(O.args)) === (identity, 1) && return O.args[1]

    (O.op, length(O.args)) === (-, 1) && return Operation(*, Expression[-1, O.args[1]])

    return O
end
_simplify_constants(x, shorten_tree) = x
_simplify_constants(x) = _simplify_constants(x, true)
