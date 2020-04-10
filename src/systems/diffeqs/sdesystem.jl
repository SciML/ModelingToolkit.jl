struct SDESystem <: AbstractODESystem
    """The expressions defining the drift term."""
    eqs::Vector{Equation}
    """The expressions defining the diffusion term."""
    noiseeqs
    """Independent variable."""
    iv::Variable
    """Dependent (state) variables."""
    dvs::Vector{Variable}
    """Parameter variables."""
    ps::Vector{Variable}
    """
    Time-derivative matrix. Note: this field will not be defined until
    [`calculate_tgrad`](@ref) is called on the system.
    """
    tgrad::RefValue{Vector{Expression}}
    """
    Jacobian matrix. Note: this field will not be defined until
    [`calculate_jacobian`](@ref) is called on the system.
    """
    jac::RefValue{Matrix{Expression}}
    """
    Wfact matrix. Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact::RefValue{Matrix{Expression}}
    """
    Wfact_t matrix. Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact_t::RefValue{Matrix{Expression}}
    """
    Name: the name of the system
    """
    name::Symbol
end

function SDESystem(deqs::AbstractVector{<:Equation}, neqs, iv, dvs, ps; name = gensym(:SDESystem))
    dvs′ = [clean(dv) for dv ∈ dvs]
    ps′ = [clean(p) for p ∈ ps]
    iv′ = clean(iv)
    tgrad = RefValue(Vector{Expression}(undef, 0))
    jac = RefValue(Matrix{Expression}(undef, 0, 0))
    Wfact   = RefValue(Matrix{Expression}(undef, 0, 0))
    Wfact_t = RefValue(Matrix{Expression}(undef, 0, 0))
    SDESystem(deqs, neqs, iv′, dvs′, ps′, tgrad, jac, Wfact, Wfact_t, name)
end

function generate_diffusion_function(sys::SDESystem, dvs = sys.dvs, ps = sys.ps, expression = Val{true}; kwargs...)
    dvs′ = [clean(dv) for dv ∈ dvs]
    ps′ = [clean(p) for p ∈ ps]
    return build_function(sys.noiseeqs, dvs′, ps′, (sys.iv.name,), ODEToExpr(sys), expression; kwargs...)
end

"""
$(SIGNATURES)

Create an `SDEFunction` from the [`SDESystem`](@ref). The arguments `dvs` and `ps`
are used to set the order of the dependent variable and parameter vectors,
respectively.
"""
function DiffEqBase.SDEFunction{iip}(sys::SDESystem, dvs = sys.dvs, ps = sys.ps;
                                     version = nothing, tgrad=false,
                                     jac = false, Wfact = false) where {iip}
    f_oop,f_iip = generate_function(sys, dvs, ps, Val{false})
    g_oop,g_iip = generate_diffusion_function(sys, dvs, ps, Val{false})

    f(u,p,t) = f_oop(u,p,t)
    f(du,u,p,t) = f_iip(du,u,p,t)
    g(u,p,t) = g_oop(u,p,t)
    g(du,u,p,t) = g_iip(du,u,p,t)

    if tgrad
        tgrad_oop,tgrad_iip = generate_tgrad(sys, dvs, ps, Val{false})
        _tgrad(u,p,t) = tgrad_oop(u,p,t)
        _tgrad(J,u,p,t) = tgrad_iip(J,u,p,t)
    else
        _tgrad = nothing
    end

    if jac
        jac_oop,jac_iip = generate_jacobian(sys, dvs, ps, Val{false})
        _jac(u,p,t) = jac_oop(u,p,t)
        _jac(J,u,p,t) = jac_iip(J,u,p,t)
    else
        _jac = nothing
    end

    if Wfact
        tmp_Wfact,tmp_Wfact_t = generate_factorized_W(sys, dvs, ps, true, Val{false})
        Wfact_oop, Wfact_iip = tmp_Wfact
        Wfact_oop_t, Wfact_iip_t = tmp_Wfact_t
        _Wfact(u,p,dtgamma,t) = Wfact_oop(u,p,dtgamma,t)
        _Wfact(W,u,p,dtgamma,t) = Wfact_iip(W,u,p,dtgamma,t)
        _Wfact_t(u,p,dtgamma,t) = Wfact_oop_t(u,p,dtgamma,t)
        _Wfact_t(W,u,p,dtgamma,t) = Wfact_iip_t(W,u,p,dtgamma,t)
    else
        _Wfact,_Wfact_t = nothing,nothing
    end

    SDEFunction{iip}(f,g,jac=_jac,
                      tgrad = _tgrad,
                      Wfact = _Wfact,
                      Wfact_t = _Wfact_t,
                      syms = Symbol.(sys.dvs))
end

function DiffEqBase.SDEFunction(sys::SDESystem, args...; kwargs...)
    SDEFunction{true}(sys, args...; kwargs...)
end
