struct Reaction
    rate
    reactants
    products
end

struct ReactionSystem
    rxs
    iv
    dvs
    ps
end

function Base.convert(::Type{<:ODESystem},rs::ReactionSystem)
    D = Differential(rs.iv)
    eqs = [D(x) ~ 0 for x in rs.dvs]

    for rx in rs.rxs
        for reactant in rx.reactants
            i = findfirst(x->x.op == reactant.op,rs.dvs)
            eqs[i] = Equation(eqs[i].lhs,eqs[i].rhs - reactant)
        end

        for product in rx.products
            i = findfirst(x->x.op == product.op,rs.dvs)
            eqs[i] = Equation(eqs[i].lhs,+(eqs[i].rhs,rx.reactants...))
        end
    end

    ODESystem(eqs,rs.iv,rs.dvs,rs.ps)
end
