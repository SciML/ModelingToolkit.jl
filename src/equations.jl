struct DiffEqSystem <: AbstractSystem
    eqs::Vector{Operation}
    ivs::Vector{Variable}
    dvs::Vector{Variable}
    vs::Vector{Variable}
    ps::Vector{Variable}
end
DiffEqSystem(eqs) = DiffEqSystem(eqs,Variable[],Variable[],Variable[],Variable[])

struct NonlinearSystem <: AbstractSystem
    eqs::Vector{Operation}
    vs::Vector{Variable}
    ps::Vector{Variable}
end
NonlinearSystem(eqs) = NonlinearSystem(eqs,Variable[],Variable[])

export DiffEqSystem,NonlinearSystem
