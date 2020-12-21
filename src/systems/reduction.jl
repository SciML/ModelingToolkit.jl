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

isvar(s::Sym; param=false) = param ? true : !isparameter(s)
isvar(s::Term; param=false) = isvar(s.op; param=param)
isvar(s::Any;param=false) = false

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
    if isvar(αx, param=true)
        return 1, αx
    elseif αx isa Term && operation(αx) === (*)
        args = arguments(αx)
        nums = filter(!isvar, args)
        syms = filter(isvar, args)

        if length(syms) == 1
            return prod(nums), syms[1]
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
        !(eq.rhs isa Symbolic)
    end
    for eq in eqs[ii]
        α,x = get_α_x(eq.lhs)
        push!(subs, x => isone(α) ? eq.rhs : eq.rhs / α)
    end
    deleteat!(eqs, ii) # remove them

    # Case 2: One side is a differentiated var, the other is an algebraic var
    #         substitute the algebraic var with the diff var
    diff_vars = filter(!isnothing, map(eqs) do eq
            if eq.lhs isa Term && eq.lhs.op isa Differential
                eq.lhs.args[1]
            else
                nothing
            end
        end) |> Set

    del = Int[]
    for (i, eq) in enumerate(eqs)
        res_left = get_α_x(eq.lhs)
        if !isnothing(res_left)
            α, x = res_left
            res_right = get_α_x(eq.rhs)
            if !isnothing(res_right)
                β, y = res_right
                if y in diff_vars && !(x in diff_vars)
                    multiple = β / α
                    push!(subs, x => isone(multiple) ? y : multiple * y)
                    push!(del, i)
                elseif x in diff_vars && !(y in diff_vars)
                    multiple = α / β
                    push!(subs, y => isone(multiple) ? x : multiple * x)
                    push!(del, i)
                end
            end
        end
    end
    deleteat!(eqs, del)

    # Case 3: Explicit substitutions
    del = Int[]
    for (i, eq) in enumerate(eqs)
        res_left = get_α_x(eq.lhs)
        if !isnothing(res_left)
            α, x = res_left
            res_right = get_α_x(eq.rhs)
            if !isnothing(res_right)
                β, y = res_right
                multiple =  β / α
                push!(subs, x => _isone(multiple) ? x : multiple * x)
                push!(del, i)
            end
        end
    end
    deleteat!(eqs, del)

    diffeqs = filter(eq -> eq.lhs isa Term && eq.lhs.op isa Differential, eqs)
    diffeqs′ = substitute_aliases(diffeqs, Dict(subs))

    newstates = map(diffeqs) do eq
        eq.lhs.args[1]
    end
    ODESystem(diffeqs′, sys.iv, newstates, parameters(sys), observed=first.(subs) .~ last.(subs))
end
