function simplify_constants(O::Operation)
    if Symbol(O.op) == :* && any(x->isequal(x,Constant(0)),O.args)
        return Constant(0)
    elseif Symbol(O.op) == :* && any(x->isequal(x,Constant(1)),O.args)
        idxs = find(x->isequal(x,Constant(1)),O.args)
        _O = Operation(O.op,O.args[1:length(O.args) .∉ (idxs,)])
        if isempty(_O.args)
            return Constant(1)
        else
            return _O
        end
    elseif Symbol(O.op) == :+ && any(x->isequal(x,Constant(0)),O.args)
        idxs = find(x->isequal(x,Constant(1)),O.args)
        _O = Operation(O.op,O.args[1:length(O.args) .∉ (idxs,)])
        if isempty(_O.args)
            return _O
        else
            return O
        end
    else
        return O
    end
end
simplify_constants(x::Variable) = x

export simplify_constants
