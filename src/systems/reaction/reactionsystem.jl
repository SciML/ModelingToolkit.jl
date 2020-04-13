struct Reaction{T <: Number}
    rate
    substrates::Vector{Operation}
    products::Vector{Operation}
    substoich::Vector{T}
    prodstoich::Vector{T}
end

struct ReactionSystem <: AbstractSystem
    eqs::Vector{Reaction}
    iv::Variable
    species::Vector{Variable}
    params::Vector{Variable}
    name::Symbol
    systems::Vector{ReactionSystem}
end

function ReactionSystem(eqs, iv, species, params;
                        systems = ReactionSystem[],
                        name = gensym(:ReactionSystem))

    ReactionSystem(eqs, iv, convert.(Variable,species), 
                   convert.(Variable,params), name, systems)
end

# Calculate the ODE rate law
function oderatelaw(rx)
    @unpack rate, substrates, substoich = rx    
    rl   = rate
    coef = one(eltype(substoich))
    for (i,stoich) in enumerate(substoich)
        coef *= factorial(stoich)        
        rl   *= (stoich != one(stoich)) ? substrates[i]^stoich : substrates[i]
    end
    (coef != one(coef)) && (rl /= coef)

    rl
end

function essemble_drift(rs)
    D   = Differential(rs.iv())
    eqs = [D(x(rs.iv())) ~ 0 for x in rs.species]
    species_to_idx = Dict(enumerate(rs.species))

    # note, this should really use the net stoichiometry to avoid 
    # adding and substracting the same term for A + X -> B + X
    for rx in rs.eqs        
        rl     = oderatelaw(rx.rate, rx.substrates, rx.substoich)
        stoich = rx.substoich
        for (sidx,substrate) in enumerate(rx.substrates)
            i      = species_to_idx[substrate.op]
            eqs[i] = Equation(eqs[i].lhs, eqs[i].rhs - stoich[sidx]*rl)
        end

        stoich = rx.prodstoich
        for (pidx,product) in enumerate(rx.products)
            i      = species_to_idx[product.op]
            eqs[i] = Equation(eqs[i].lhs,eqs[i].rhs + stoich[pidx]*rl)
        end
    end

    eqs
end

function essemble_diffusion(rs)
    eqs = Expression[Constant(0) for x in rs.species, y in rs.eqs]
    species_to_idx = Dict(enumerate(rs.species))

    for (j,rx) in enumerate(rs.eqs)
        rlsqrt = sqrt(oderatelaw(rx.rate, rx.substrates, rx.substoich))
        stoich = rx.substoich
        for (sidx,substrate) in enumerate(rx.substrates)
            i        = species_to_idx[substrate.op]
            eqs[i,j] = -stoich[sidx] * rlsqrt
        end

        for (pidx,product) in enumerate(rx.products)
            i        = species_to_idx[product.op]
            Δp       = stoich[pidx] * rlsqrt            
            eqs[i,j] = (eqs[i,j]==Constant(0)) ? Δp : (eqs[i,j] + Δp)
        end
    end
    eqs
end

function Base.convert(::Type{<:ODESystem},rs::ReactionSystem)
    eqs = essemble_drift(rs)
    ODESystem(eqs,rs.iv,rs.species,rs.params,name=rs.name,
              systems=convert.(ODESystem,rs.systems))
end

function Base.convert(::Type{<:SDESystem},rs::ReactionSystem)
    eqs = essemble_drift(rs)
    noiseeqs = essemble_diffusion(rs)
    SDESystem(eqs,noiseeqs,rs.iv,rs.species,rs.params,
              name=rs.name,systems=convert.(SDESystem,rs.systems))
end
