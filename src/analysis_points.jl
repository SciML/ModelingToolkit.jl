using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq
using ModelingToolkit: get_eqs, vars, @set!, get_iv

## Implementation
Base.@kwdef mutable struct AnalysisPoint
    in = nothing
    out = nothing
    name::Symbol
end

function AnalysisPoint(in, out; name)
    AnalysisPoint(in, out, name)
end

function AnalysisPoint(name)
    AnalysisPoint(; name)
end

Base.show(io::IO, ap::AnalysisPoint) = show(io, MIME"text/plain"(), ap)
function Base.show(io::IO, ::MIME"text/plain", ap::AnalysisPoint)
    if get(io, :compact, false)
        print(io, "AnalysisPoint($(ap.in.u), $(ap.out.u); name=$(ap.name))")
    else
        print(io, "AnalysisPoint(")
        printstyled(io, ap.name, color=:cyan)
        if ap.in !== nothing && ap.out !== nothing
            print(io, " from ")
            printstyled(io, ap.in.u, color=:green)
            print(io, " to ")
            printstyled(io, ap.out.u, color=:blue)
        end
        print(io, ")")
    end
end

function ModelingToolkit.connect(in, ap::AnalysisPoint, out)
    ap.in = in
    ap.out = out
    return 0 ~ ap
end

ModelingToolkit.vars(ap::AnalysisPoint; op = Differential) = vars(connect(ap.in, ap.out); op)

"Replace analysis points with the identity connection connect(ap.in, ap.out). This is used when a system containing analysis points is simulated, in which case analysis points have no effect."
function expand_analysis_points(sys)
    new_eqs = map(get_eqs(sys)) do eq
        eq.rhs isa AnalysisPoint || (return eq)
        ap = eq.rhs
        connect(ap.in, ap.out)
    end
    @set! sys.eqs = new_eqs
    sys
end

function get_sensitivity(sys, ap; kwargs...)
    t = get_iv(sys)
    @variables d(t)=0 # Perturbantion serving as input to sensivity transfer function
    new_eqs = map(get_eqs(sys)) do eq
        eq.rhs == ap || (return eq)
        ap.out.u ~ ap.in.u + d # This assumes that the connector as an internal vaiable named u
    end
    # new_eqs = reduce(vcat, new_eqs)
    @set! sys.eqs = new_eqs
    @set! sys.states = [states(sys); d]
    ModelingToolkit.linearize(sys, [d], [ap.out.u]; kwargs...)
end


function get_comp_sensitivity(sys, ap; kwargs...)
    t = get_iv(sys)
    @variables d(t)=0 # Perturbantion serving as input to sensivity transfer function
    new_eqs = map(get_eqs(sys)) do eq
        eq.rhs == ap || (return eq)
        ap.out.u +d ~ ap.in.u # This assumes that the connector as an internal vaiable named u
    end
    # new_eqs = reduce(vcat, new_eqs)
    @set! sys.eqs = new_eqs
    @set! sys.states = [states(sys); d]
    ModelingToolkit.linearize(sys, [d], [ap.in.u]; kwargs...)
end

## Test




@named P = FirstOrder(k=1, T=1)
@named C = Gain(-1)

ap = AnalysisPoint(:plant_input)

eqs = [
    connect(P.output, C.input)
    connect(C.output, ap, P.input)
]

t = ModelingToolkit.get_iv(P)

sys = ODESystem(eqs, t, systems=[P,C], name=:hej)


ssys = structural_simplify(expand_analysis_points(sys))
prob = ODEProblem(ssys, [P.x => 1], (0, 10))
sol = solve(prob, Rodas5())
@test norm(sol[1]) >= 1
@test norm(sol[end]) < 1e-6 # This fails without the feedback through C
# plot(sol)

matrices, _ = get_sensitivity(sys, ap)
@test matrices.A[] == -2
@test matrices.B[]*matrices.C[] == -1 # either one negative
@test matrices.D[] == 1


matrices, _ = get_comp_sensitivity(sys, ap)
@test matrices.A[] == -2
@test matrices.B[]*matrices.C[] == 1 # both positive
@test matrices.D[] == 0

#=
# Equivalent code using ControlSystemsBase
using ControlSystemsBase
P = tf(1.0, [1, 1])
C = 1
S = sensitivity(P, C)      # or feedback(1, P*C)
T = comp_sensitivity(P, C) # or feedback(P*C)
=#