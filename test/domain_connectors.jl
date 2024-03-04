using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Test

@connector function HydraulicPort(; p_int, name)
    pars = @parameters begin
        ρ
        β
        μ
    end

    vars = @variables begin
        p(t) = p_int
        dm(t), [connect = Flow]
    end

    ODESystem(Equation[], t, vars, pars; name, defaults = [dm => 0])
end

@connector function HydraulicFluid(;
        density = 997,
        bulk_modulus = 2.09e9,
        viscosity = 0.0010016,
        name)
    pars = @parameters begin
        ρ = density
        β = bulk_modulus
        μ = viscosity
    end

    vars = @variables begin
        dm(t), [connect = Flow]
    end

    eqs = [
        dm ~ 0
    ]

    ODESystem(eqs, t, vars, pars; name, defaults = [dm => 0])
end

function FixedPressure(; p, name)
    pars = @parameters begin
        p = p
    end

    vars = []

    systems = @named begin
        port = HydraulicPort(; p_int = p)
    end

    eqs = [
        port.p ~ p
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

function FixedVolume(; vol, p_int, name)
    pars = @parameters begin
        p_int = p_int
        vol = vol
    end

    systems = @named begin
        port = HydraulicPort(; p_int)
    end

    vars = @variables begin
        rho(t) = port.ρ
        drho(t) = 0
    end

    # let
    dm = port.dm
    p = port.p

    eqs = [D(rho) ~ drho
           rho ~ port.ρ * (1 + p / port.β)
           dm ~ drho * vol]

    ODESystem(eqs, t, vars, pars; name, systems)
end

function Valve2Port(; p_s_int, p_r_int, p_int, name)
    pars = @parameters begin
        p_s_int = p_s_int
        p_r_int = p_r_int
        p_int = p_int
        x_int = 0
        scale = 1.0

        k = 0.1
    end

    systems = @named begin
        HS = HydraulicPort(; p_int = p_s_int)
        HR = HydraulicPort(; p_int = p_r_int)
        port = HydraulicPort(; p_int)
    end

    vars = @variables begin
        x(t) = x_int
    end

    # let (flow) ---------
    Δp_s = HS.p - port.p
    Δp_r = port.p - HR.p

    x̃ = abs(x / scale)
    Δp̃_s = abs(Δp_s)
    Δp̃_r = abs(Δp_r)

    flow(Δp̃) = (k) * (Δp̃) * (x̃)

    # 
    eqs = [domain_connect(port, HS, HR)
           port.dm ~ -ifelse(x >= 0, +flow(Δp̃_s), -flow(Δp̃_r))
           HS.dm ~ ifelse(x >= 0, port.dm, 0)
           HR.dm ~ ifelse(x < 0, port.dm, 0)]

    ODESystem(eqs, t, vars, pars; name, systems)
end

function System(; name)
    vars = []
    pars = []
    systems = @named begin
        fluid = HydraulicFluid(; density = 500, bulk_modulus = 1e9, viscosity = 0.01)
        src = FixedPressure(; p = 200)
        rtn = FixedPressure(; p = 0)
        valve = Valve2Port(; p_s_int = 200, p_r_int = 0, p_int = 100)
        vol = FixedVolume(; vol = 0.1, p_int = 100)
    end
    eqs = [domain_connect(fluid, src.port)
           connect(src.port, valve.HS)
           connect(rtn.port, valve.HR)
           connect(vol.port, valve.port)
           valve.x ~ sin(2π * t * 10)]

    return ODESystem(eqs, t, vars, pars; systems, name)
end

@named odesys = System()
esys = ModelingToolkit.expand_connections(odesys)
@test length(equations(esys)) == length(unknowns(esys))

csys = complete(odesys)

sys = structural_simplify(odesys)
@test length(equations(sys)) == length(unknowns(sys))

sys_defs = ModelingToolkit.defaults(sys)
@test Symbol(sys_defs[csys.vol.port.ρ]) == Symbol(csys.fluid.ρ)
