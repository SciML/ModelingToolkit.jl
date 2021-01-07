export alias_elimination, flatten

function flatten(sys::ODESystem)
    if isempty(sys.systems)
        return sys
    else
        return ODESystem(equations(sys),
                         independent_variable(sys),
                         states(sys),
                         parameters(sys),
                         observed=observed(sys))
    end
end


using SymbolicUtils: Rewriters

function fixpoint_sub(x, dict)
    y = substitute(x, dict)
    while !isequal(x, y)
        y = x
        x = substitute(y, dict)
    end

    return x
end

function substitute_aliases(diffeqs, dict)
    lhss(diffeqs) .~ fixpoint_sub.(rhss(diffeqs), (dict,))
end

# Note that we reduce parameters, too
# i.e. `2param = 3` will be reduced away
isvar(s::Sym) = true
isvar(s::Term) = isvar(operation(s))
isvar(s::Any) = false

function get_α_x(αx)
    if isvar(αx)
        return 1, αx
    elseif αx isa Term && operation(αx) === (*)
        args = arguments(αx)
        nums = []
        syms = []
        for arg in args
            isvar(arg) ? push!(syms, arg) : push!(nums, arg)
        end

        if length(syms) == 1
            return prod(nums), syms[1]
        end
    else
        return nothing
    end
end

function alias_elimination(sys::ODESystem)
    eqs = vcat(equations(sys), observed(sys))
    subs = Pair[]
    diff_vars = filter(!isnothing, map(eqs) do eq
            if isdiffeq(eq)
                eq.lhs.args[1]
            else
                nothing
            end
        end) |> Set

    # only substitute when the variable is algebraic
    del = Int[]
    for (i, eq) in enumerate(eqs)
        isdiffeq(eq) && continue
        res_left = get_α_x(eq.lhs)
        if !isnothing(res_left) && !(res_left[2] in diff_vars)
            # `α x = rhs` => `x = rhs / α`
            α, x = res_left
            push!(subs, x => _isone(α) ? eq.rhs : eq.rhs / α)
            push!(del, i)
        else
            res_right = get_α_x(eq.rhs)
            if !isnothing(res_right) && !(res_right[2] in diff_vars)
                # `lhs = β y` => `y = lhs / β`
                β, y = res_right
                push!(subs, y => _isone(β) ? eq.lhs : β * eq.lhs)
                push!(del, i)
            end
        end
    end
    deleteat!(eqs, del)

    eqs′ = substitute_aliases(eqs, Dict(subs))
    alias_vars = first.(subs)

    newstates = setdiff(states(sys), alias_vars)
    ODESystem(eqs′, sys.iv, newstates, parameters(sys), observed=alias_vars .~ last.(subs))
end
