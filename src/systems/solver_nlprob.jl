function MTKBase.generate_ODENLStepData(
        sys::System, u0, p, mm = calculate_massmatrix(sys),
        nlstep_compile::Bool = true, nlstep_scc::Bool = false;
        jac::Bool = false, limit_bounds::Bool = true
    )
    # Under `nlstep_scc` the whole corrector mechanism is unavailable, and the bounds clamp
    # is automatic rather than requested: deriving it only to reject the problem would break
    # bounded models that decompose perfectly well. A declared `limited` quantity still
    # errors below, since there the user did ask for limiting.
    limit_bounds &= !nlstep_scc
    nlsys, outer_tmp, inner_tmp = inner_nlsystem(sys, mm, nlstep_compile; limit_bounds)
    return assemble_nlstep_data(sys, nlsys, outer_tmp, inner_tmp, u0, p, nlstep_scc, jac)
end

function MTKBase.generate_DAENLStepData(
        sys::System, u0, p, mm = calculate_massmatrix(sys),
        nlstep_compile::Bool = true, nlstep_scc::Bool = false;
        jac::Bool = false
    )
    # The fully implicit stage system reuses the ODE unknowns as increment symbols just as
    # the mass-matrix one does, so `bounds` metadata is no more applicable here. Limiting is
    # not wired through this path yet, and `inner_dae_nlsystem` derives no box, which keeps
    # the DAE behaviour as it was before `limited` existed.
    nlsys, outer_tmp, inner_tmp = inner_dae_nlsystem(sys, mm, nlstep_compile)
    return assemble_nlstep_data(sys, nlsys, outer_tmp, inner_tmp, u0, p, nlstep_scc, jac)
end

function assemble_nlstep_data(
        sys::System, nlsys::System, outer_tmp, inner_tmp, u0, p,
        nlstep_scc::Bool, jac::Bool
    )
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
    # Forward the outer function's `jac` flag so the inner teared nonlinear
    # problem also carries an analytic Jacobian. This lets NonlinearSolve.jl
    # skip the per-iteration AD/FD Jacobian computation on the inner Newton.
    nlprob = if nlstep_scc
        if getmetadata(nlsys, MTKBase.LimitedCtx, nothing) !== nothing
            throw(
                ArgumentError(
                    "`limited` quantities are not supported with `nlstep_scc = true`: the " *
                        "stage limiters compile into the `postcondition` of a single " *
                        "`NonlinearProblem`, which the SCC decomposition splits apart."
                )
            )
        end
        SCCNonlinearProblem(nlsys, op; build_initializeprob = false, jac)
    else
        # The stage unknowns are the Newton increments `z`, not the states: an unknown's
        # `bounds` metadata is a box on the physical `γ₂ z + inner_tmp`, so deriving a
        # static `lb`/`ub` from it here would constrain the wrong variable — and would put
        # the problem into NonlinearSolve's transformed coordinates, which the raw stage
        # values the stepper writes into `nlprob.u0` do not live in. `limit_bounds`
        # delivers the box to the stage solve correctly instead.
        NonlinearProblem(
            nlsys, op; build_initializeprob = false, jac, lb = nothing, ub = nothing
        )
    end

    subsetidxs = [findfirst(isequal(y), unknowns(sys)) for y in unknowns(nlsys)]
    set_gamma_c = setsym(nlsys, (ODE_GAMMA..., ODE_C))
    set_outer_tmp = setsym(nlsys, outer_tmp)
    set_inner_tmp = setsym(nlsys, inner_tmp)
    nlprobmap = generate_nlprobmap(sys, nlsys)

    return SciMLBase.ODENLStepData(
        nlprob, subsetidxs, set_gamma_c, set_outer_tmp, set_inner_tmp, nlprobmap
    )
end

const ODE_GAMMA = @parameters γ₁ₘₜₖ, γ₂ₘₜₖ, γ₃ₘₜₖ
const ODE_C = only(@parameters cₘₜₖ)

function get_outer_tmp(n::Int)
    return only(@parameters outer_tmpₘₜₖ[1:n])
end

function get_inner_tmp(n::Int)
    return only(@parameters inner_tmpₘₜₖ[1:n])
end

function inner_nlsystem(sys::System, mm, nlstep_compile::Bool; limit_bounds::Bool = true)
    dvs = unknowns(sys)
    eqs = full_equations(sys)
    t = get_iv(sys)
    N = length(dvs)
    @assert length(eqs) == N
    @assert mm isa UniformScaling || size(mm) == (N, N)
    rhss = [eq.rhs for eq in eqs]
    gamma1, gamma2, gamma3 = ODE_GAMMA
    c = ODE_C
    outer_tmp = get_outer_tmp(N)
    inner_tmp = get_inner_tmp(N)

    subrules = Dict([v => unwrap(gamma2 * v + inner_tmp[i]) for (i, v) in enumerate(dvs)])
    subrules[t] = unwrap(c)
    new_rhss = map(Base.Fix2(substitute, subrules), rhss)
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
    # `mtkcompile` strips `limited` from a time-dependent system but records what it
    # stripped: the stage system is where those limiters become a `postcondition`.
    nlsys = MTKBase.attach_stage_limiters(sys, nlsys, subrules; limit_bounds)
    return nlsys, outer_tmp, inner_tmp
end

"""
    $(TYPEDSIGNATURES)

Build the stage system of the fully implicit form `0 = F(du, u, p, t)`, i.e.
`g(z, p') = F(γ₁ z + outer_tmp, γ₂ z + inner_tmp, p, c)`.

`full_equations` leaves a compiled system's differential equations as `D(x) ~ rhs` and its
algebraic ones as `0 ~ rhs`, so the residual `rhs - lhs` that `DAEFunction` generates is
`rhs - mm * du`: row `i` of `mm` picks out the `du` component that equation differentiates,
if any. Substituting the derivative argument is therefore subtracting `mm * (γ₁ z +
outer_tmp)`. An algebraic component has an all-zero `mm` column, so its `outer_tmp` entry
drops out — correct, since `F` does not depend on that component's derivative.

`γ₃` scales the `M z` term of the mass-matrix stage system; the fully implicit form has no
such term, so it is absent here and consumers must leave it at `1`. It stays in the
parameter list only so that `set_γ_c` takes the same four values on both paths.
"""
function inner_dae_nlsystem(sys::System, mm, nlstep_compile::Bool)
    dvs = unknowns(sys)
    eqs = full_equations(sys)
    t = get_iv(sys)
    N = length(dvs)
    @assert length(eqs) == N
    @assert mm isa UniformScaling || size(mm) == (N, N)
    rhss = [eq.rhs for eq in eqs]
    gamma1, gamma2, gamma3 = ODE_GAMMA
    c = ODE_C
    outer_tmp = get_outer_tmp(N)
    inner_tmp = get_inner_tmp(N)

    subrules = Dict([v => unwrap(gamma2 * v + inner_tmp[i]) for (i, v) in enumerate(dvs)])
    subrules[t] = unwrap(c)
    new_rhss = map(Base.Fix2(substitute, subrules), rhss)
    new_rhss = new_rhss .- mm * (gamma1 .* dvs .+ collect(outer_tmp))
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
    return nlp.f(buffer, state_values(nlsol), parameter_values(nlsol))
end

function (nlp::NLStep_probmap)(nlsol)
    return nlp.f(state_values(nlsol), parameter_values(nlsol))
end

function generate_nlprobmap(sys::System, nlsys::System)
    return NLStep_probmap(
        build_explicit_observed_function(
            nlsys, unknowns(sys), GeneratedFunctionOptions(; expression = Val{false})
        )
    )
end
