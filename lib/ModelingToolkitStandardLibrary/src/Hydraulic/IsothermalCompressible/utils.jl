# regPow(x, a, delta = 0.01) = x * (x * x + delta * delta)^((a - 1) / 2);
function regPow(x, a, delta = 0.01)
    ifelse(abs(x / delta) >= 1, sign(x) * abs(x / delta)^a * delta^a, (delta^a * x) / delta)
end
regRoot(x, delta = 0.01) = regPow(x, 0.5, delta)

"""
    HydraulicPort(; name)

Connector port for hydraulic components.

# States:
- `p`: [Pa] gauge total pressure
- `dm`: [kg/s] mass flow
"""
@connector function HydraulicPort(; name)
    pars = @parameters begin
        ρ
        β
        μ
        n
        let_gas
        ρ_gas
        p_gas
    end

    vars = @variables begin
        p(t), [guess = 0]
        dm(t), [guess = 0, connect = Flow]
    end

    System(Equation[], t, vars, pars; name)
end

"""
    HydraulicFluid(; density = 997, bulk_modulus = 2.09e9, viscosity = 0.0010016, gas_density = 0.0073955, gas_pressure = -1000, n = 1, let_gas = 1, name)

Fluid parameter setter for isothermal compressible fluid domain.  Defaults given for water at 20°C and 0Pa gage (1atm absolute) reference pressure. Density is modeled using the Tait equation of state.  For pressures below the reference pressure, density is linearly interpolated to the gas state (when `let_gas` is set to 1), this helps prevent pressures from going below the reference pressure.

# Parameters:

- `ρ`: [kg/m^3] fluid density at 0Pa reference gage pressure (set by `density` argument)
- `Β`: [Pa] fluid bulk modulus describing the compressibility (set by `bulk_modulus` argument)
- `μ`: [Pa*s] or [kg/m-s] fluid dynamic viscosity  (set by `viscosity` argument)
- `n`: density exponent
- `let_gas`: set to 1 to allow fluid to transition from liquid to gas (for density calculation only)
- `ρ_gas`: [kg/m^3] density of fluid in gas state at reference gage pressure `p_gas` (set by `gas_density` argument)
- `p_gas`: [Pa] reference pressure (set by `gas_pressure` argument)
"""
@connector function HydraulicFluid(; density = 997, bulk_modulus = 2.09e9,
        viscosity = 0.0010016, gas_density = 0.0073955,
        gas_pressure = -1000, n = 1, let_gas = 1, name)
    pars = @parameters begin
        ρ = density
        β = bulk_modulus
        μ = viscosity
        n = n
        let_gas = let_gas
        ρ_gas = gas_density
        p_gas = gas_pressure
    end

    vars = @variables begin
        dm(t), [guess = 0, connect = Flow]
    end

    eqs = [
        dm ~ 0
    ]

    System(eqs, t, vars, pars; name)
end

function transition(x1, x2, y1, y2, x)
    u = (x - x1) / (x2 - x1)
    blend = u^2 * (3 - 2 * u)
    return (1 - blend) * y1 + blend * y2
end

f_laminar(shape_factor, Re) = shape_factor * regPow(Re, -1, 0.1) #regPow used to avoid dividing by 0, min value is 0.1
f_turbulent(shape_factor, Re) = (shape_factor / 64) / (0.79 * log(Re) - 1.64)^2

"""
    friction_factor(dm, area, d_h, viscosity, shape_factor)

Calculates the friction factor ``f`` for fully developed flow in a tube such that ``Δp = f \\cdot \\rho \\frac{u^2}{2} \\frac{l}{d_h}`` where

- ``Δp``: [Pa] is the pressure difference over the tube length ``l``
- ``\\rho``: [kg/m^3] is the average fluid density
- ``u``: [m/s] is the average fluid velocity
- ``l``: [m] is the tube length

The friction factor is calculated for laminar and turbulent flow with a transition region between Reynolds number 2000 to 3000.  Turbulent flow equation is for smooth tubes, valid for the Reynolds number range up to 5e6.

# Arguments:

- `dm`: [kg/s] mass flow
- `area`: [m^2] tube cross sectional area
- `d_h`: [m] tube hydraulic diameter.  For circular tubes d_h is the tube diameter, otherwise it can be found from `4*area/perimeter`
- `density`: [kg/m^3] fluid density
- `viscosity`: [Pa*s] or [kg/m-s] fluid dynamic viscosity
- `shape_factor`: the constant defining the laminar fully developed constant f*Re related to the shape of the tube cross section

Reference: Introduction to Fluid Mechanics, Fox & McDonald, 5th Edition, equations 8.19 and 8.21
"""
function friction_factor(dm, area, d_h, viscosity, shape_factor)
    # u = abs(dm) / (density * area)
    # Re = density * u * d_h / viscosity

    Re = abs(dm) * d_h / (area * viscosity)

    if Re <= 2000
        return f_laminar(shape_factor, Re)
    elseif 2000 < Re < 3000
        return transition(2000, 3000, f_laminar(shape_factor, Re),
            f_turbulent(shape_factor, Re), Re)
    else
        return f_turbulent(shape_factor, Re)
    end
end
@register_symbolic friction_factor(dm, area, d_h, viscosity, shape_factor)
Symbolics.derivative(::typeof(friction_factor), args, ::Val{1}) = 0
Symbolics.derivative(::typeof(friction_factor), args, ::Val{4}) = 0

density_ref(port) = port.ρ
density_exp(port) = port.n
gas_density_ref(port) = port.ρ_gas
gas_pressure_ref(port) = port.p_gas
bulk_modulus(port) = port.β
viscosity(port) = port.μ

function liquid_density(port, p)
    density_ref(port) *
    regPow(1 + density_exp(port) * p / bulk_modulus(port), 1 / density_exp(port))
end #Tait-Murnaghan equation of state
liquid_density(port) = liquid_density(port, port.p)

# p = beta*(rho/rho_0 - 1)
# (p/beta + 1)*rho_0 = rho

function liquid_pressure(port, rho)
    (rho / density_ref(port) - 1) * bulk_modulus(port)
end

function gas_density(port, p)
    slope = (density_ref(port) - gas_density_ref(port)) / (0 - gas_pressure_ref(port))
    b = density_ref(port)

    return b + p * slope
end

function gas_pressure(port, rho)
    slope = (0 - gas_pressure_ref(port)) / (density_ref(port) - gas_density_ref(port))
    b = 0

    return b + rho * slope
end

function full_density(port, p)
    ifelse(port.let_gas == 1,
        ifelse(p >= 0, liquid_density(port, p), gas_density(port, p)),
        liquid_density(port, p))
end
full_density(port) = full_density(port, port.p)

function full_pressure(port, rho)
    ifelse(port.let_gas == 1,
        ifelse(
            rho >= density_ref(port), liquid_pressure(port, rho), gas_pressure(port, rho)),
        liquid_pressure(port, rho)
    )
end
