abstract type AbstractOptimalControlProblem{uType, tType, isinplace} <:
              SciMLBase.AbstractODEProblem{uType, tType, isinplace} end

struct OptimalControlSolution
    model::Any
    sol::ODESolution
    input_sol::Union{Nothing, ODESolution}
end

function JuMPControlProblem end
function InfiniteOptControlProblem end
function CasADiControlProblem end
function PyomoControlProblem end

function warn_overdetermined(sys, u0map)
    constraintsys = get_constraintsystem(sys)
    if !isnothing(constraintsys)
        (length(constraints(constraintsys)) + length(u0map) > length(unknowns(sys))) &&
            @warn "The control problem is overdetermined. The total number of conditions (# constraints + # fixed initial values given by u0map) exceeds the total number of states. The solvers will default to doing a nonlinear least-squares optimization."
    end
end
