export alias_elimination

function flatten(sys::ODESystem)
    if isempty(sys.systems)
        return sys
    else
        return ODESystem(equations(sys),
                         independent_variable(sys),
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

isvar(s::Sym) = !isparameter(s)
isvar(s::Term) = isvar(s.op)
isvar(s::Any) = false

function filterexpr(f, s)
    vs = []
    Rewriters.Prewalk(Rewriters.Chain([@rule((~x::f) => push!(vs, ~x))]))(s)
    vs
end

function make_lhs_0(eq)
    if eq.lhs isa Number && iszero(eq.lhs)
        return eq
    else
        0 ~ eq.lhs - eq.rhs
    end
end

function alias_elimination(sys::ODESystem)
    eqs = vcat(equations(sys), observed(sys))

    # make all algebraic equations have 0 on LHS
    eqs = map(eqs) do eq
        if eq.lhs isa Term && eq.lhs.op isa Differential
            eq
        else
            make_lhs_0(eq)
        end
    end

    newstates = map(eqs) do eq
        if eq.lhs isa Term && eq.lhs.op isa Differential
            filterexpr(isvar, eq.lhs)
        else
            []
        end
    end |> Iterators.flatten |> collect |> unique


    all_vars = map(eqs) do eq
        filterexpr(isvar, eq.rhs)
    end |> Iterators.flatten |> collect |> unique

    alg_idxs = findall(x->!(x.lhs isa Term) && iszero(x.lhs), eqs)

    eliminate = setdiff(all_vars, newstates)

    outputs = solve_for(eqs[alg_idxs], eliminate)

    diffeqs = eqs[setdiff(1:length(eqs), alg_idxs)]

    diffeqs′ = substitute_aliases(diffeqs, Dict(eliminate .=> outputs))

    ODESystem(diffeqs′, sys.iv, newstates, parameters(sys), observed=eliminate .~ outputs)
end

function get_α_x(αx)
    if isvar(αx)
        return αx, 1
    elseif αx isa Term && operation(αx) === (*)
        args = arguments(αx)
        nums = filter(!isvar, args)
        syms = filter(isvar, args)

        if length(syms) == 1
            return syms[1], prod(nums)
        end
    else
        return nothing
    end
end

function alias_elimination2(sys)
    eqs = vcat(equations(sys), observed(sys))

    subs = Pair[]
    # Case 1: Right hand side is a constant
    ii = findall(eqs) do eq
        (eq.lhs isa Sym || (eq.lhs isa Term && !(eq.lhs.op isa Differential))) && !(eq.rhs isa Symbolic)
    end
    for eq in eqs[ii]
        substitution_dict[eq.lhs] = eq.rhs
        push!(subs, eq.lhs => eq.rhs)
    end
    deleteat!(eqs, ii) # remove them

    # Case 2: One side is a differentiated var, the other is an algebraic var
    #         substitute the algebraic var with the diff var
    diff_vars = findall(eqs) do eq
        if eq.lhs isa Term && eq.lhs.op isa Differential
            eq.lhs.args[1]
        else
            nothing
        end
    end

    for eq in eqs
        res_left = get_α_x(eq.lhs)
        if !isnothing(res)
            res_right = get_α_x(eq.rhs)
            β, y = res
            if y in diff_vars && !(x in diff_vars)
                multiple = β / α
                push!(subs, x => isone(multiple) ? y : multiple * y)
            elseif x in diff_vars && !(y in diff_vars)
                multiple = α / β
                push!(subs, y => isone(multiple) ? y : multiple * y)
            end
        end
    end

    # Case 3: Explicit substitutions
    for eq in eqs
        res_left = get_α_x(eq.lhs)
        if !isnothing(res)
            res_right = get_α_x(eq.rhs)
            β, y = res
            multiple =  β / α
            push!(subs, x => isone(multiple) ? x : multiple * x)
        end
    end

    diffeqs = filter(eq -> eq.lhs isa Term && eq.lhs.op isa Differential, eqs)
    diffeqs′ = substitute_aliases(diffeqs, Dict(subs))
    ODESystem(diffeqs′, sys.iv, newstates, parameters(sys), observed=first.(subs) .~ last.(subs))
end
