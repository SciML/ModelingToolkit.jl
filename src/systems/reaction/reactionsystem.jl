struct Reaction
    rate
    reactants::Vector{Operation}
    products::Vector{Operation}
end

struct ReactionSystem <: AbstractSystem
    eqs::Vector{Reaction}
    iv::Variable
    states::Vector{Variable}
    ps::Vector{Variable}
    name::Symbol
    systems::Vector{ReactionSystem}
end

function ReactionSystem(eqs,iv,dvs,ps;
                        systems = ReactionSystem[],
                        name=gensym(:ReactionSystem))
    ReactionSystem(eqs,iv,convert.(Variable,dvs),convert.(Variable,ps),name,systems)
end

# TODO: Make it do the combinatorics stuff
reaction_expr(reactants) = *(reactants...)

function essemble_drift(rs)
    D = Differential(rs.iv())
    eqs = [D(x(rs.iv())) ~ 0 for x in rs.states]

    for rx in rs.eqs
        for reactant in rx.reactants
            i = findfirst(x->x == reactant.op,rs.states)
            eqs[i] = Equation(eqs[i].lhs,eqs[i].rhs - rx.rate * reaction_expr(rx.reactants))
        end

        for product in rx.products
            i = findfirst(x->x == product.op,rs.states)
            eqs[i] = Equation(eqs[i].lhs,eqs[i].rhs + rx.rate * reaction_expr(rx.reactants))
        end
    end
    eqs
end

function essemble_diffusion(rs)
    eqs = Expression[Constant(0) for x in rs.states, y in rs.eqs]

    for (j,rx) in enumerate(rs.eqs)
        for reactant in rx.reactants
            i = findfirst(x->x == reactant.op,rs.states)
            eqs[i,j] = -sqrt(rx.rate) * reaction_expr(rx.reactants)
        end

        for product in rx.products
            i = findfirst(x->x == product.op,rs.states)
            eqs[i,j] = sqrt(rx.rate) * reaction_expr(rx.reactants)
        end
    end
    eqs
end

function Base.convert(::Type{<:ODESystem},rs::ReactionSystem)
    eqs = essemble_drift(rs)
    ODESystem(eqs,rs.iv,rs.states,rs.ps,name=rs.name,
              systems=convert.(ODESystem,rs.systems))
end

function Base.convert(::Type{<:SDESystem},rs::ReactionSystem)
    eqs = essemble_drift(rs)
    noiseeqs = essemble_diffusion(rs)
    SDESystem(eqs,noiseeqs,rs.iv,rs.states,rs.ps,
              name=rs.name,systems=convert.(SDESystem,rs.systems))
end
