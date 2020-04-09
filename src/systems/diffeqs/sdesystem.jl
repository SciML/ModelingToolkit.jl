struct SDESystem <: AbstractODESystem
    """The expressions defining the drift term."""
    eqs::Vector{ODEExpr}
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
end

function SDESystem(deqs::AbstractVector{ODEExpr}, neqs, iv, dvs, ps)
    tgrad = RefValue(Vector{Expression}(undef, 0))
    jac = RefValue(Matrix{Expression}(undef, 0, 0))
    Wfact   = RefValue(Matrix{Expression}(undef, 0, 0))
    Wfact_t = RefValue(Matrix{Expression}(undef, 0, 0))
    SDESystem(deqs, neqs, iv, dvs, ps, tgrad, jac, Wfact, Wfact_t)
end

function SDESystem(deqs::AbstractVector{<:Equation}, neqs, iv, dvs, ps)
    _dvs = [deq.op for deq ∈ dvs]
    _iv = iv.op
    _ps = [p.op for p ∈ ps]
    SDESystem(getindex.(convert.(ODEExpr,deqs),2), neqs, _iv, _dvs, _ps)
end

function generate_diffusion_function(sys::SDESystem, dvs = sys.dvs, ps = sys.ps, expression = Val{true}; kwargs...)
    dvs′ = [clean(dv) for dv ∈ dvs]
    ps′ = [clean(p) for p ∈ ps]
    return build_function(sys.noiseeqs, dvs′, ps′, (sys.iv.name,), ODEToExpr(sys), expression; kwargs...)
end
