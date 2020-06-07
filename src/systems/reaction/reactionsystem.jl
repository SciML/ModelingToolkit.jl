
"""
$(TYPEDEF)

One chemical reaction.

# Fields
$(FIELDS)

# Examples

```
using ModelingToolkit
@parameters t k[1:20]
@variables A(t) B(t) C(t) D(t)
rxs = [Reaction(k[1], nothing, [A]),            # 0 -> A
       Reaction(k[2], [B], nothing),            # B -> 0
       Reaction(k[3],[A],[C]),                  # A -> C
       Reaction(k[4], [C], [A,B]),              # C -> A + B
       Reaction(k[5], [C], [A], [1], [2]),      # C -> A + A
       Reaction(k[6], [A,B], [C]),              # A + B -> C
       Reaction(k[7], [B], [A], [2], [1]),      # 2B -> A
       Reaction(k[8], [A,B], [A,C]),            # A + B -> A + C
       Reaction(k[9], [A,B], [C,D]),            # A + B -> C + D
       Reaction(k[10], [A], [C,D], [2], [1,1]), # 2A -> C + D
       Reaction(k[11], [A], [A,B], [2], [1,1]), # 2A -> A + B
       Reaction(k[12], [A,B,C], [C,D], [1,3,4], [2, 3]),          # A+3B+4C -> 2C + 3D
       Reaction(k[13], [A,B], nothing, [3,1], nothing),           # 3A+B -> 0
       Reaction(k[14], nothing, [A], nothing, [2]),               # 0 -> 2A
       Reaction(k[15]*A/(2+A), [A], nothing; only_use_rate=true), # A -> 0 with custom rate
       Reaction(k[16], [A], [B]; only_use_rate=true),             # A -> B with custom rate.
       Reaction(k[17]*A*exp(B), [C], [D], [2], [1]),              # 2C -> D with non constant rate.
       Reaction(k[18]*B, nothing, [B], nothing, [2]),             # 0 -> 2B with non constant rate.
       Reaction(k[19]*t, [A], [B]),                                # A -> B with non constant rate.
       Reaction(k[20]*t*A, [B,C], [D],[2,1],[2])                  # 2A +B -> 2C with non constant rate.
  ]
```

Notes:
- `nothing` can be used to indicate a reaction that has no reactants or no products.
  In this case the corresponding stoichiometry vector should also be set to `nothing`.
- The three-argument form assumes all reactant and product stoichiometric coefficients
  are one.
"""
struct Reaction{S <: Variable, T <: Number}
    """The rate function (excluding mass action terms)."""
    rate
    """Reaction substrates."""
    substrates::Vector{Operation}
    """Reaction products."""
    products::Vector{Operation}
    """The stoichiometric coefficients of the reactants."""
    substoich::Vector{T}
    """The stoichiometric coefficients of the products."""
    prodstoich::Vector{T}
    """The net stoichiometric coefficients of all species changed by the reaction."""
    netstoich::Vector{Pair{S,T}}
    """
    `false` (default) if `rate` should be multiplied by mass action terms to give the rate law.
    `true` if `rate` represents the full reaction rate law.
    """
    only_use_rate::Bool
end

function Reaction(rate, subs, prods, substoich, prodstoich;
                  netstoich=nothing, only_use_rate=false, kwargs...)

      (isnothing(prods)&&isnothing(subs)) && error("A reaction requires a non-nothing substrate or product vector.")
      (isnothing(prodstoich)&&isnothing(substoich)) && error("Both substrate and product stochiometry inputs cannot be nothing.")
      if isnothing(subs)
        subs = Vector{Operation}()
        (substoich!=nothing) && error("If substrates are nothing, substrate stiocihometries have to be so too.")
        substoich = typeof(prodstoich)()
    end
    if isnothing(prods)
        prods = Vector{Operation}()
        (prodstoich!=nothing) && error("If products are nothing, product stiocihometries have to be so too.")
        prodstoich = typeof(substoich)()
    end
    ns = isnothing(netstoich) ? get_netstoich(subs, prods, substoich, prodstoich) : netstoich
    Reaction(rate, subs, prods, substoich, prodstoich, ns, only_use_rate)
end


# three argument constructor assumes stoichiometric coefs are one and integers
function Reaction(rate, subs, prods; kwargs...)

    sstoich = isnothing(subs) ? nothing : ones(Int,length(subs))
    pstoich = isnothing(prods) ? nothing : ones(Int,length(prods))
    Reaction(rate, subs, prods, sstoich, pstoich; kwargs...)
end

# calculates the net stoichiometry of a reaction as a vector of pairs (sub,substoich)
function get_netstoich(subs, prods, sstoich, pstoich)
    # stoichiometry as a Dictionary
    nsdict = Dict{Variable,eltype(sstoich)}(sub.op => -sstoich[i] for (i,sub) in enumerate(subs))
    for (i,p) in enumerate(prods)
        coef = pstoich[i]
        prod = p.op
        @inbounds nsdict[prod] = haskey(nsdict, prod) ? nsdict[prod] + coef : coef
    end

    # stoichiometry as a vector
    ns = [el for el in nsdict if el[2] != zero(el[2])]

    ns
end

"""
$(TYPEDEF)

A system of chemical reactions.

# Fields
$(FIELDS)

# Example
Continuing from the example in the [`Reaction`](@ref) definition:
```
rs = ReactionSystem(rxs, t, [A,B,C,D], k)
```
"""
struct ReactionSystem <: AbstractSystem
    """The reactions defining the system."""
    eqs::Vector{Reaction}
    """Independent variable (usually time)."""
    iv::Variable
    """Dependent (state) variables representing amount of each species."""
    states::Vector{Variable}
    """Parameter variables."""
    ps::Vector{Variable}
    """The name of the system"""
    name::Symbol
    """systems: The internal systems"""
    systems::Vector{ReactionSystem}
end

function ReactionSystem(eqs, iv, species, params; systems = ReactionSystem[],
                                                  name = gensym(:ReactionSystem))


    isempty(species) && error("ReactionSystems require at least one species.")
    paramvars = isempty(params) ? Variable[] : convert.(Variable, params)
    ReactionSystem(eqs, iv, convert.(Variable,species), paramvars, name, systems)
end

# Calculate the ODE rate law
function oderatelaw(rx)
    @unpack rate, substrates, substoich, only_use_rate = rx
    rl = rate
    if !only_use_rate
        coef = one(eltype(substoich))
        for (i,stoich) in enumerate(substoich)
            coef *= factorial(stoich)
            rl   *= isone(stoich) ? substrates[i] : substrates[i]^stoich
        end
        (!isone(coef)) && (rl /= coef)
    end
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
                signedrl = (stoich > zero(stoich)) ? rl : -rl
                rhs      = isone(abs(stoich)) ? signedrl : stoich * rl
            else
                Δspec = isone(abs(stoich)) ? rl : abs(stoich) * rl
                rhs   = (stoich > zero(stoich)) ? (eqs[i].rhs + Δspec) : (eqs[i].rhs - Δspec)
            end
            eqs[i] = Equation(eqs[i].lhs, rhs)
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
            i            = species_to_idx[spec]
            signedrlsqrt = (stoich > zero(stoich)) ? rlsqrt : -rlsqrt
            eqs[i,j]     = isone(abs(stoich)) ? signedrlsqrt : stoich * rlsqrt
        end
    end
    eqs
end

function var2op(var)
    Operation(var,Vector{Expression}())
end

# Calculate the Jump rate law (like ODE, but uses X instead of X(t).
# The former generates a "MethodError: objects of type Int64 are not callable" when trying to solve the problem.
function jumpratelaw(rx; rxvars=get_variables(rx.rate))
    @unpack rate, substrates, substoich, only_use_rate = rx
    rl = rate
    for op in rxvars
        rl = substitute(rl, op => var2op(op.op))
    end
    if !only_use_rate
        for (i,stoich) in enumerate(substoich)
            rl *= isone(stoich) ? var2op(substrates[i].op) : Operation(binomial,[var2op(substrates[i].op),stoich])
        end
    end
    rl
end

# if haveivdep=false then time dependent rates will still be classified as mass action
"""
```julia
ismassaction(rx, rs; rxvars = get_variables(rx.rate),
                              haveivdep = any(var -> isequal(rs.iv,convert(Variable,var)), rxvars))
```

True if a given reaction is of mass action form, i.e. `rx.rate` does not depend
on any chemical species that correspond to states of the system, and does not depend
explicitly on the independent variable (usually time).

# Arguments
- `rx`, the [`Reaction`](@ref).
- `rs`, a [`ReactionSystem`](@ref) containing the reaction.
- Optional: `rxvars`, `Variable`s which are not in `rxvars` are ignored as possible dependencies.
- Optional: `haveivdep`, `true` if the [`Reaction`](@ref) `rate` field explicitly depends on the independent variable.
"""
function ismassaction(rx, rs; rxvars = get_variables(rx.rate),
                              haveivdep = any(var -> isequal(rs.iv,convert(Variable,var)), rxvars))
    # if no dependencies must be zero order
    if isempty(rxvars)
        return true
    else
        return !(haveivdep || rx.only_use_rate || any(convert(Variable,rxv) in states(rs) for rxv in rxvars))
    end
end

function assemble_jumps(rs)
    eqs = Vector{Union{ConstantRateJump, MassActionJump, VariableRateJump}}()

    for rx in equations(rs)
        rxvars    = (rx.rate isa Operation) ? get_variables(rx.rate) : Operation[]
        haveivdep = any(var -> isequal(rs.iv,convert(Variable,var)), rxvars)
        if ismassaction(rx, rs; rxvars=rxvars, haveivdep=haveivdep)
            reactant_stoch = isempty(rx.substoich) ? [0 => 1] : [var2op(sub.op) => stoich for (sub,stoich) in zip(rx.substrates,rx.substoich)]
            coef           = isempty(rx.substoich) ? one(eltype(rx.substoich)) : prod(stoich -> factorial(stoich), rx.substoich)
            rate           = isone(coef) ? rx.rate : rx.rate/coef
            net_stoch      = [Pair(var2op(p[1]),p[2]) for p in rx.netstoich]
            push!(eqs, MassActionJump(rate, reactant_stoch, net_stoch, scale_rates=false))
        else
            rl     = jumpratelaw(rx, rxvars=rxvars)
            affect = Vector{Equation}()
            for (spec,stoich) in rx.netstoich
                push!(affect, var2op(spec) ~ var2op(spec) + stoich)
            end
            if haveivdep
                push!(eqs, VariableRateJump(rl,affect))
            else
                push!(eqs, ConstantRateJump(rl,affect))
            end
        end
    end
    eqs
end

"""
```julia
Base.convert(::Type{<:ODESystem},rs::ReactionSystem)
```
Convert a ReactionSystem to an ODESystem.
"""
function Base.convert(::Type{<:ODESystem},rs::ReactionSystem)
    eqs = assemble_drift(rs)
    ODESystem(eqs,rs.iv,rs.states,rs.ps,name=rs.name,
              systems=convert.(ODESystem,rs.systems))
end

"""
```julia
Base.convert(::Type{<:SDESystem},rs::ReactionSystem)
```

Convert a ReactionSystem to a SDESystem.
"""
function Base.convert(::Type{<:SDESystem},rs::ReactionSystem)
    eqs = assemble_drift(rs)
    noiseeqs = assemble_diffusion(rs)
    SDESystem(eqs,noiseeqs,rs.iv,rs.states,rs.ps,
              name=rs.name,systems=convert.(SDESystem,rs.systems))
end

"""
```julia
Base.convert(::Type{<:JumpSystem},rs::ReactionSystem)
```

Convert a ReactionSystem to a JumpSystem.
"""
function Base.convert(::Type{<:JumpSystem},rs::ReactionSystem)
    eqs = assemble_jumps(rs)
    JumpSystem(eqs,rs.iv,rs.states,rs.ps,name=rs.name,
              systems=convert.(JumpSystem,rs.systems))
end


"""
```julia
Base.convert(::Type{<:NonlinearSystem},rs::ReactionSystem)
```

Convert a ReactionSystem to a  NonlinearSystem.
"""
function Base.convert(::Type{<:NonlinearSystem},rs::ReactionSystem)
    states_swaps = map(states -> Operation(states,[var2op(rs.iv)]), rs.states)
    eqs = map(eq -> 0 ~ make_sub!(eq,states_swaps),getproperty.(assemble_drift(rs),:rhs))
    NonlinearSystem(eqs,rs.states,rs.ps,name=rs.name,
              systems=convert.(NonlinearSystem,rs.systems))
end

# Used for Base.convert(::Type{<:NonlinearSystem},rs::ReactionSystem) only, should likely be removed.
function make_sub!(eq,states_swaps)
	for (i,arg) in enumerate(eq.args)
		if any(isequal.(states_swaps,arg))
			eq.args[i] = var2op(arg.op)
		else
			make_sub!(arg,states_swaps)
		end
	end
	return eq
end

### Converts a reaction system to ODE or SDE problems ###


# ODEProblem from AbstractReactionNetwork
function DiffEqBase.ODEProblem(rs::ReactionSystem, u0::Union{AbstractArray, Number}, tspan, p, args...; kwargs...)
    u0 = typeof(u0) <: Array{<:Pair} ? u0 : Pair.(rs.states,u0)
    p = typeof(p) <: Array{<:Pair} ? p : Pair.(rs.ps,p)
    return ODEProblem(convert(ODESystem,rs),u0,tspan,p, args...; kwargs...)
end

# SDEProblem from AbstractReactionNetwork
function DiffEqBase.SDEProblem(rs::ReactionSystem, u0::Union{AbstractArray, Number}, tspan, p, args...; kwargs...)
    u0 = typeof(u0) <: Array{<:Pair} ? u0 : Pair.(rs.states,u0)
    p = typeof(p) <: Array{<:Pair} ? p : Pair.(rs.ps,p)
    p_matrix = zeros(length(rs.states), length(rs.eqs))
    return SDEProblem(convert(SDESystem,rs),u0,tspan,p,args...; noise_rate_prototype=p_matrix,kwargs...)
end

# DiscreteProblem from AbstractReactionNetwork
function DiffEqBase.DiscreteProblem(rs::ReactionSystem, u0::Union{AbstractArray, Number}, tspan::Tuple, p=nothing, args...; kwargs...)
    u0 = typeof(u0) <: Array{<:Pair} ? u0 : Pair.(rs.states,u0)
    p = typeof(p) <: Array{<:Pair} ? p : Pair.(rs.ps,p)
    return DiscreteProblem(convert(JumpSystem,rs), u0,tspan,p, args...; kwargs...)
end

# JumpProblem from AbstractReactionNetwork
function DiffEqJump.JumpProblem(rs::ReactionSystem, prob, aggregator, args...; kwargs...)
    return JumpProblem(convert(JumpSystem,rs), prob, aggregator, args...; kwargs...)
end

# SteadyStateProblem from AbstractReactionNetwork
function DiffEqBase.SteadyStateProblem(rs::ReactionSystem, u0::Union{AbstractArray, Number}, p, args...; kwargs...)
    #u0 = typeof(u0) <: Array{<:Pair} ? u0 : Pair.(rs.states,u0)
    #p = typeof(p) <: Array{<:Pair} ? p : Pair.(rs.ps,p)
    return SteadyStateProblem(ODEFunction(convert(ODESystem,rs)),u0,p, args...; kwargs...)
end
