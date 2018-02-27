struct DiffEqSystem <: AbstractSystem
    eqs::Vector{Statement}
    ivs::Vector{IndependentVariable}
    dvs::Vector{DependentVariable}
    vs::Vector{Variable}
    ps::Vector{Parameter}
end

struct NonlinearSystem <: AbstractSystem
    eqs::Vector{Statement}
    vs::Vector{Variable}
    ps::Vector{Parameter}
end

export DiffEqSystem,NonlinearSystem
