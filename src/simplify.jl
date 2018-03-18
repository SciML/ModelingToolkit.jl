function simplify_constants(O::Operation)
    O_last = nothing
    _O = O
    while _O != O_last
        O_last = _O
        _O = _simplify_constants(_O)
        if typeof(_O) <: Operation
            _O = Operation(_O.op,simplify_constants.(_O.args))
        end
    end
    _O
end

function _simplify_constants(O)
    if Symbol(O.op) == :*
        # If any variable is `Constant(0)`, zero the whole thing
        # If any variable is `Constant(1)`, remove that `Constant(1)` unless
        # they are all `Constant(1)`, in which case simplify to a single variable
        if any(x->typeof(x)<:Variable && isequal(x,Constant(0)),O.args)
            return Constant(0)
        elseif any(x->typeof(x)<:Variable && isequal(x,Constant(1)),O.args)
            idxs = find(x->typeof(x)<:Variable && isequal(x,Constant(1)),O.args)
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
    elseif Symbol(O.op) == :+ && any(x->typeof(x)<:Variable && isequal(x,Constant(0)),O.args)
        # If there are Constant(0)s in a big `+` expression, get rid of them
        idxs = find(x->typeof(x)<:Variable && isequal(x,Constant(0)),O.args)
        _O = Operation(O.op,O.args[1:length(O.args) .∉ (idxs,)])
        if isempty(_O.args)
            return Constant(0)
        elseif length(_O.args) == 1
            return _O.args[1]
        else
            return O
        end
    else
        return O
    end
end
simplify_constants(x::Variable) = x
_simplify_constants(x::Variable) = x

export simplify_constants
