function simplify_constants(O::Operation, shorten_tree = true)
    O_last = nothing
    _O = O
    while _O != O_last
        O_last = _O
        _O = _simplify_constants(_O,shorten_tree)
        if is_operation(_O)
            _O = Operation(_O.op,simplify_constants.(_O.args,shorten_tree))
        end
    end
    _O
end

const AC_OPERATORS = [*, +]

function _simplify_constants(O, shorten_tree = true)
    # Tree shrinking
    if shorten_tree && O.op ∈ AC_OPERATORS
        # Flatten tree
        idxs = findall(x -> is_operation(x) && x.op == O.op, O.args)
        if !isempty(idxs)
            keep_idxs = eachindex(O.args) .∉ Ref(idxs)
            args = Vector{Expression}[O.args[i].args for i in idxs]
            push!(args, O.args[keep_idxs])
            return Operation(O.op, vcat(args...))
        end

        # Collapse constants
        idxs = findall(is_constant, O.args)
        if length(idxs) > 1
            other_idxs = eachindex(O.args) .∉ (idxs,)
            new_var = Constant(mapreduce(get, O.op, O.args[idxs]))
            new_args = O.args[other_idxs]
            push!(new_args,new_var)

            return length(new_args) > 1 ? Operation(O.op, new_args) : first(new_args)
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

    if O.op === (+) && any(iszero, O.args)
        # If there are Constant(0)s in a big `+` expression, get rid of them
        args = filter(!iszero, O.args)

        isempty(args)     && return Constant(0)
        length(args) == 1 && return first(args)
        return Operation(O.op, args)
    end

    (O.op, length(O.args)) === (identity, 1) && return O.args[1]

    (O.op, length(O.args)) === (-, 1) && return Operation(*, Expression[-1, O.args[1]])

    return O
end
simplify_constants(x::Variable,y=false) = x
_simplify_constants(x::Variable,y=false) = x

export simplify_constants
