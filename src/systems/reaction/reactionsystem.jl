struct Reaction{S <: Variable, T <: Number}
    rate
    substrates::Vector{Operation}
    products::Vector{Operation}
    substoich::Vector{T}
    prodstoich::Vector{T}    
    netstoich::Vector{Pair{S,T}}
    only_use_rate::Bool
end

function Reaction(rate, subs, prods, substoich, prodstoich; 
                  netstoich=nothing, only_use_rate=false, kwargs...)

    ns = isnothing(netstoich) ? get_netstoich(subs, prods, substoich, prodstoich) : netstoich
    Reaction(rate, subs, prods, substoich, prodstoich, ns, only_use_rate)
end

# three argument constructor assumes stoichiometric coefs are one and integers
Reaction(rate, subs, prods; kwargs...) = Reaction(rate, subs, prods,
                                                  ones(Int,length(subs)),
                                                  ones(Int,length(prods)); 
                                                  kwargs...)

# calculates the net stoichiometry of a reaction as a vector of pairs (sub,substoich)
function get_netstoich(subs, prods, sstoich, pstoich)
    # stoichiometry as a Dictionary
    nsdict = Dict(sub.op => -sstoich[i] for (i,sub) in enumerate(subs))
    for (i,p) in enumerate(prods)
        coef = pstoich[i]
        prod = p.op
        @inbounds nsdict[prod] = haskey(nsdict, prod) ? nsdict[prod] + coef : coef
    end

    # stoichiometry as a vector
    ns = [el for el in nsdict if el[2] != zero(el[2])]

    ns
end


struct ReactionSystem <: AbstractSystem
    eqs::Vector{Reaction}
    iv::Variable
    states::Vector{Variable}
    ps::Vector{Variable}
    name::Symbol
    systems::Vector{ReactionSystem}
end

function ReactionSystem(eqs, iv, species, params; systems = ReactionSystem[],
                                                  name = gensym(:ReactionSystem))

    ReactionSystem(eqs, iv, convert.(Variable,species), convert.(Variable,params), 
                   name, systems)
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

function assemble_drift(rs)
    D   = Differential(rs.iv())
    eqs = [D(x(rs.iv())) ~ 0 for x in rs.states]
    species_to_idx = Dict((x => i for (i,x) in enumerate(rs.states)))

    for rx in rs.eqs        
        rl = oderatelaw(rx)
        for (spec,stoich) in rx.netstoich
            i = species_to_idx[spec]
            if iszero(eqs[i].rhs)
                Δspec  = (stoich == one(stoich)) ? rl : stoich * rl            
                eqs[i] = Equation(eqs[i].lhs, Δspec)
            else
                Δspec = (abs(stoich) == one(stoich)) ? rl : abs(stoich) * rl            
                if stoich > zero(stoich)
                    eqs[i] = Equation(eqs[i].lhs, eqs[i].rhs + Δspec)
                else
                    eqs[i] = Equation(eqs[i].lhs, eqs[i].rhs - Δspec)
                end
            end
        end
    end
    eqs
end

function assemble_diffusion(rs)
    eqs = Expression[Constant(0) for x in rs.states, y in rs.eqs]
    species_to_idx = Dict((x => i for (i,x) in enumerate(rs.states)))

    for (j,rx) in enumerate(rs.eqs)
        rlsqrt = sqrt(oderatelaw(rx))
        for (spec,stoich) in rx.netstoich
            i        = species_to_idx[spec]
            eqs[i,j] = (stoich == one(stoich)) ? rlsqrt : stoich * rlsqrt            
        end
    end
    eqs
end

function Base.convert(::Type{<:ODESystem},rs::ReactionSystem)
    eqs = assemble_drift(rs)
    ODESystem(eqs,rs.iv,rs.states,rs.ps,name=rs.name,
              systems=convert.(ODESystem,rs.systems))
end

function Base.convert(::Type{<:SDESystem},rs::ReactionSystem)
    eqs = assemble_drift(rs)
    noiseeqs = assemble_diffusion(rs)
    SDESystem(eqs,noiseeqs,rs.iv,rs.states,rs.ps,
              name=rs.name,systems=convert.(SDESystem,rs.systems))
end
