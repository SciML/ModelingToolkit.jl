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

