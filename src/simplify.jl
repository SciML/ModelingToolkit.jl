function simplify_constants(O::Operation,shorten_tree = true)
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

const TREE_SHRINK_OPS = [*, +]

function _simplify_constants(O,shorten_tree = true)
    # Tree shrinking
    if shorten_tree
        for cur_op in TREE_SHRINK_OPS
            if O.op == cur_op
                # Shrink tree
                if any(x -> is_operation(x) && x.op == cur_op, O.args)
                    idxs = findall(x -> is_operation(x) && x.op == cur_op, O.args)
                    keep_idxs = 1:length(O.args) .∉ (idxs,)
                    args = Vector{Expression}[O.args[i].args for i in idxs]
                    push!(args,O.args[keep_idxs])
                    return Operation(O.op, vcat(args...))
                end
                # Collapse constants
                idxs = findall(is_constant, O.args)
                if length(idxs) > 1
                    other_idxs = 1:length(O.args) .∉ (idxs,)
                    if cur_op == (*)
                        new_var = Constant(prod(get, O.args[idxs]))
                    elseif cur_op == (+)
                        new_var = Constant(sum(get, O.args[idxs]))
                    end
                    new_args = O.args[other_idxs]
                    push!(new_args,new_var)
                    if length(new_args) > 1
                        return Operation(O.op,new_args)
                    else
                        return new_args[1]
                    end
                end
            end
        end
    end

    if O.op == (*)
        # If any variable is `Constant(0)`, zero the whole thing
        # If any variable is `Constant(1)`, remove that `Constant(1)` unless
        # they are all `Constant(1)`, in which case simplify to a single variable
        if any(iszero, O.args)
            return Constant(0)
        elseif any(isone, O.args)
            idxs = findall(isone, O.args)
            _O = Operation(O.op,O.args[1:length(O.args) .∉ (idxs,)])
            if isempty(_O.args)
                return Constant(1)
            elseif length(_O.args) == 1
                return _O.args[1]
            else
                return _O
            end
        else
            return O
        end
    elseif O.op == (+) && any(iszero, O.args)
        # If there are Constant(0)s in a big `+` expression, get rid of them
        idxs = findall(iszero, O.args)
        _O = Operation(O.op,O.args[1:length(O.args) .∉ (idxs,)])
        if isempty(_O.args)
            return Constant(0)
        elseif length(_O.args) == 1
            return _O.args[1]
        else
            return O
        end
    elseif O.op == identity
        return O.args[1]
    elseif O.op == (-) && length(O.args) == 1
        return Operation(*,Expression[-1,O.args[1]])
    else
        return O
    end
    return O
end
simplify_constants(x::Variable,y=false) = x
_simplify_constants(x::Variable,y=false) = x

export simplify_constants
