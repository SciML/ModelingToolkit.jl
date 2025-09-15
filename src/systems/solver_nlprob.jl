function generate_ODENLStepData(sys::System, u0, p, mm = calculate_massmatrix(sys),
        nlstep_compile::Bool = true, nlstep_scc::Bool = false)
    nlsys, outer_tmp, inner_tmp = inner_nlsystem(sys, mm, nlstep_compile)
    state = ProblemState(; u = u0, p)
    op = Dict()
    op[ODE_GAMMA[1]] = one(eltype(u0))
    op[ODE_GAMMA[2]] = one(eltype(u0))
    op[ODE_GAMMA[3]] = one(eltype(u0))
    op[ODE_C] = zero(eltype(u0))
    op[outer_tmp] = zeros(eltype(u0), size(outer_tmp))
    op[inner_tmp] = zeros(eltype(u0), size(inner_tmp))
    for v in [unknowns(nlsys); parameters(nlsys)]
        haskey(op, v) && continue
        op[v] = getsym(sys, v)(state)
    end
    nlprob = if nlstep_scc
        SCCNonlinearProblem(nlsys, op; build_initializeprob = false)
    else
        NonlinearProblem(nlsys, op; build_initializeprob = false)
    end

    subsetidxs = [findfirst(isequal(y), unknowns(sys)) for y in unknowns(nlsys)]
    set_gamma_c = setsym(nlsys, (ODE_GAMMA..., ODE_C))
    set_outer_tmp = setsym(nlsys, outer_tmp)
    set_inner_tmp = setsym(nlsys, inner_tmp)
    nlprobmap = generate_nlprobmap(sys, nlsys)

    return SciMLBase.ODENLStepData(
        nlprob, subsetidxs, set_gamma_c, set_outer_tmp, set_inner_tmp, nlprobmap)
end

const ODE_GAMMA = @parameters γ₁ₘₜₖ, γ₂ₘₜₖ, γ₃ₘₜₖ
const ODE_C = only(@parameters cₘₜₖ)

function get_outer_tmp(n::Int)
    only(@parameters outer_tmpₘₜₖ[1:n])
end

function get_inner_tmp(n::Int)
    only(@parameters inner_tmpₘₜₖ[1:n])
end

function inner_nlsystem(sys::System, mm, nlstep_compile::Bool)
    dvs = unknowns(sys)
    eqs = full_equations(sys)
    t = get_iv(sys)
    N = length(dvs)
    @assert length(eqs) == N
    @assert mm == I || size(mm) == (N, N)
    rhss = [eq.rhs for eq in eqs]
    gamma1, gamma2, gamma3 = ODE_GAMMA
    c = ODE_C
    outer_tmp = get_outer_tmp(N)
    inner_tmp = get_inner_tmp(N)

    subrules = Dict([v => unwrap(gamma2*v + inner_tmp[i]) for (i, v) in enumerate(dvs)])
    subrules[t] = unwrap(c)
    new_rhss = map(Base.Fix2(fast_substitute, subrules), rhss)
    new_rhss = collect(outer_tmp) .+ gamma1 .* new_rhss .- gamma3 * mm * dvs
    new_eqs = [0 ~ rhs for rhs in new_rhss]

    new_dvs = unknowns(sys)
    new_ps = [parameters(sys); [gamma1, gamma2, gamma3, c, inner_tmp, outer_tmp]]
    nlsys = System(new_eqs, new_dvs, new_ps; name = :nlsys)
    nlsys = if nlstep_compile
        mtkcompile(nlsys; split = is_split(sys))
    else
        complete(nlsys; split = is_split(sys))
    end
    return nlsys, outer_tmp, inner_tmp
end

struct NLStep_probmap{F}
    f::F
end

function (nlp::NLStep_probmap)(buffer, nlsol)
    nlp.f(buffer, state_values(nlsol), parameter_values(nlsol))
end

function (nlp::NLStep_probmap)(nlsol)
    nlp.f(state_values(nlsol), parameter_values(nlsol))
end

function generate_nlprobmap(sys::System, nlsys::System)
    return NLStep_probmap(build_explicit_observed_function(nlsys, unknowns(sys)))
end
