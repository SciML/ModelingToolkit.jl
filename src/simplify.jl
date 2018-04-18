function simplify_constants(O::Operation,shorten_tree = true)
    O_last = nothing
    _O = O
    while _O != O_last
        O_last = _O
        _O = _simplify_constants(_O,shorten_tree)
        if typeof(_O) <: Operation
            _O = Operation(_O.op,simplify_constants.(_O.args,shorten_tree))
        end
    end
    _O
end

const TREE_SHRINK_OPS = [:*,:+]

function _simplify_constants(O,shorten_tree = true)
    # Tree shrinking
    if shorten_tree
      for cur_op in TREE_SHRINK_OPS
        if Symbol(O.op) == cur_op
          # Shrink tree
          if any(x->typeof(x)<:Operation && Symbol(x.op) == cur_op ,O.args)
              idxs = find(x->typeof(x)<:Operation && Symbol(x.op) == cur_op,O.args)
              keep_idxs = 1:length(O.args) .∉ (idxs,)
              args = Vector{Expression}[O.args[i].args for i in idxs]
              push!(args,O.args[keep_idxs])
              return Operation(O.op,vcat(args...))
          # Collapse constants
          elseif length(find(x->typeof(x)<:Variable && x.subtype == :Constant ,O.args)) > 1
              idxs = find(x->typeof(x)<:Variable && x.subtype == :Constant ,O.args)
              other_idxs = 1:length(O.args) .∉ (idxs,)
              if cur_op == :*
                new_var = Constant(prod(x->x.value,O.args[idxs]))
              elseif cur_op == :+
                new_var = Constant(sum(x->x.value,O.args[idxs]))
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

    if Symbol(O.op) == :*
        # If any variable is `Constant(0)`, zero the whole thing
        # If any variable is `Constant(1)`, remove that `Constant(1)` unless
        # they are all `Constant(1)`, in which case simplify to a single variable
        if any(x->typeof(x)<:Variable && (isequal(x,Constant(0)) || isequal(x,Constant(-0))),O.args)
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
    elseif Symbol(O.op) == :+ && any(x->typeof(x)<:Variable &&
           (isequal(x,Constant(0)) || isequal(x,Constant(-0))),O.args)
        # If there are Constant(0)s in a big `+` expression, get rid of them
        idxs = find(x->typeof(x)<:Variable && (isequal(x,Constant(0)) || isequal(x,Constant(-0))),O.args)
        _O = Operation(O.op,O.args[1:length(O.args) .∉ (idxs,)])
        if isempty(_O.args)
            return Constant(0)
        elseif length(_O.args) == 1
            return _O.args[1]
        else
            return O
        end
    #=
    elseif O.op == identity
        return O.args[1]
    =#
    elseif Symbol(O.op) == :- && length(O.args) == 1
        return Operation(*,Expression[-1,O.args[1]])
    else
        return O
    end
    return O
end
simplify_constants(x::Variable,y=false) = x
_simplify_constants(x::Variable,y=false) = x

export simplify_constants
