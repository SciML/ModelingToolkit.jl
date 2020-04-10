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

# TODO: Make it do the combinatorics stuff
reaction_expr(reactants) = *(reactants...)

function essemble_drift(rs)
    D = Differential(rs.iv)
    eqs = [D(x) ~ 0 for x in rs.dvs]

    for rx in rs.rxs
        for reactant in rx.reactants
            i = findfirst(x->x.op == reactant.op,rs.dvs)
            eqs[i] = Equation(eqs[i].lhs,eqs[i].rhs - rx.rate * reaction_expr(rx.reactants))
        end

        for product in rx.products
            i = findfirst(x->x.op == product.op,rs.dvs)
            eqs[i] = Equation(eqs[i].lhs,eqs[i].rhs + rx.rate * reaction_expr(rx.reactants))
        end
    end
    eqs
end

function essemble_diffusion(rs)
    eqs = Expression[Constant(0) for x in rs.dvs, y in rs.rxs]
    @show size(eqs)

    for (j,rx) in enumerate(rs.rxs)
        for reactant in rx.reactants
            i = findfirst(x->x.op == reactant.op,rs.dvs)
            eqs[i,j] -= sqrt(rx.rate) * reaction_expr(rx.reactants)
        end

        for product in rx.products
            i = findfirst(x->x.op == product.op,rs.dvs)
            eqs[i,j] += sqrt(rx.rate) * reaction_expr(rx.reactants)
        end
    end
    eqs
end

function Base.convert(::Type{<:ODESystem},rs::ReactionSystem)
    eqs = essemble_drift(rs)
    ODESystem(eqs,rs.iv,rs.dvs,rs.ps)
end

function Base.convert(::Type{<:SDESystem},rs::ReactionSystem)
    D = Differential(rs.iv)
    eqs = essemble_drift(rs)
    noiseeqs = essemble_diffusion(rs)
    SDESystem(eqs,noiseeqs,rs.iv,rs.dvs,rs.ps)
end
