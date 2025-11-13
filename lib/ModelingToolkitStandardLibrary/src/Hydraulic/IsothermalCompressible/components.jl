
"""
    Cap(; name)

Caps a hydraulic port to prevent mass flow in or out.

# Connectors:
- `port`: hydraulic port
"""
@mtkmodel Cap begin
    @variables begin
        p(t), [guess = 0]
    end

    @components begin
        port = HydraulicPort()
    end

    @equations begin
        port.p ~ p
        port.dm ~ 0
    end
end

"""
    Open(; name)

Provides an "open" boundary condition for a hydraulic port such that mass flow `dm` is non-zero.  This is opposite from an un-connected hydraulic port or the `Cap` boundary component which sets the mass flow `dm` to zero.

# Connectors:
- `port`: hydraulic port
"""
@mtkmodel Open begin
    @variables begin
        p(t), [guess = 0]
        dm(t), [guess = 0]
    end

    @components begin
        port = HydraulicPort()
    end

    @equations begin
        port.p ~ p
        port.dm ~ dm
    end
end

"""
    TubeBase(add_inertia = true, variable_length = true; area, length_int, head_factor = 1, perimeter = 2 * sqrt(area * pi), shape_factor = 64, name)

Variable length internal flow model of the fully developed incompressible flow friction.  Includes optional inertia term when `add_inertia = true` to model wave propagation.  Hydraulic ports have equal flow but variable pressure.  Density is averaged over the pressures, used to calculated average flow velocity and flow friction.

# States:
- `x`: [m] length of the pipe
- `ddm`: [kg/s^2] Rate of change of mass flow rate in control volume.

# Parameters:
- `area`: [m^2] tube cross sectional area
- `length_int`: [m] initial tube length
- `perimeter`: [m] perimeter of the pipe cross section (needed only for non-circular pipes)
- `shape_factor`: shape factor, see `friction_factor` function
- `head_factor`: effective length multiplier, used to account for addition friction from flow development and additional friction such as pipe bends, entrance/exit lossses, etc.

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
"""
@component function TubeBase(add_inertia = true, variable_length = true;
        area,
        length_int,
        head_factor = 1,
        perimeter = 2 * sqrt(area * pi),
        shape_factor = 64,
        name)
    pars = @parameters begin
        area = area
        length_int = length_int
        perimeter = perimeter
        shape_factor = shape_factor
        head_factor = head_factor
    end

    @variables begin
        x(t), [guess = length_int]
        ddm(t) = 0
    end

    vars = []
    if variable_length
        push!(vars, x)
        c = x
    else
        c = length_int
    end
    add_inertia && push!(vars, ddm)

    systems = @named begin
        port_a = HydraulicPort()
        port_b = HydraulicPort()
    end

    # let ----------------------
    Δp = port_a.p - port_b.p
    dm = port_a.dm

    d_h = 4 * area / perimeter

    # Opting for a more numerically stable constant density (use head factor to compensate if needed)
    ρ = density_ref(port_a)  # (full_density(port_a) + full_density(port_b)) / 2
    μ = viscosity(port_a)

    f = friction_factor(dm, area, d_h, μ, shape_factor)
    u = dm / (ρ * area)

    shear = (1 / 2) * ρ * regPow(u, 2) * f * head_factor * (c / d_h)
    inertia = if add_inertia
        (c / area) * ddm
    else
        0
    end

    eqs = [0 ~ port_a.dm + port_b.dm
           domain_connect(port_a, port_b)]

    if variable_length
        push!(eqs, Δp ~ ifelse(c > 0, shear + inertia, zero(c)))
    else
        push!(eqs, Δp ~ shear + inertia)
    end

    if add_inertia
        push!(eqs, D(dm) ~ ddm)
    end

    System(eqs, t, vars, pars; name, systems)
end

"""
    Tube(N, add_inertia=true; p_int, area, length, head_factor=1, perimeter = 2 * sqrt(area * pi), shape_factor = 64, name)

Constant length internal flow model discretized by `N` (`FixedVolume`: `N`, `TubeBase`:`N-1`) which models the fully developed flow friction, compressibility (when `N>1`), and inertia effects when `add_inertia = true`.  See `TubeBase` and `FixedVolume` for more information.

# Parameters:
- `p_int`: [Pa] initial pressure
- `area`: [m^2] tube cross sectional area
- `length`: [m] real length of the tube
- `perimeter`: [m] perimeter of the pipe cross section (needed only for non-circular pipes)
- `shape_factor`: shape factor, see `friction_factor` function
- `head_factor`: effective length multiplier, used to account for addition friction from flow development and additional friction such as pipe bends, entrance/exit lossses, etc.

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
"""
@component function Tube(N, add_inertia = true; area, length, head_factor = 1,
        perimeter = 2 * sqrt(area * pi),
        shape_factor = 64, p_int, name)
    @assert(N>0,
        "the Tube component must be defined with at least 1 segment (i.e. N>0), found N=$N")

    #TODO: How to set an assert effective_length >= length ??
    pars = @parameters begin
        area = area
        length = length
        head_factor = head_factor
        perimeter = perimeter
        shape_factor = shape_factor
        p_int = p_int
    end

    vars = []

    ports = @named begin
        port_a = HydraulicPort()
        port_b = HydraulicPort()
    end

    pipe_bases = []
    for i in 1:N
        x = TubeBase(add_inertia, false; name = Symbol("p$i"),
            shape_factor = ParentScope(shape_factor),
            area = ParentScope(area),
            length_int = N > 1 ? ParentScope(length) / (N - 1) : ParentScope(length),
            head_factor = ParentScope(head_factor),
            perimeter = ParentScope(perimeter))
        push!(pipe_bases, x)
    end

    eqs = [connect(pipe_bases[1].port_a, port_a)
           connect(pipe_bases[end].port_b, port_b)]

    volumes = []
    for i in 1:(N - 1)
        x = FixedVolume(; name = Symbol("v$i"),
            vol = ParentScope(area) * ParentScope(length) / (N - 1),
            p_int = ParentScope(p_int))
        push!(volumes, x)
        push!(eqs,
            connect(x.port, pipe_bases[i].port_b, pipe_bases[i + 1].port_a))
    end

    return System(eqs, t, vars, pars; name, systems = [ports; pipe_bases; volumes])
end

"""
    FlowDivider(; n, name)

Reduces the flow from `port_a` to `port_b` by `n`.  Useful for modeling parallel tubes efficiently by placing a `FlowDivider` on each end of a tube.

# Parameters:
- `n`: divide flow from `port_a` to `port_b` by `n`

# Connectors:
- `port_a`: full flow hydraulic port
- `port_b`: part flow hydraulic port
"""
@mtkmodel FlowDivider begin

    #TODO: assert n >= 1

    @parameters begin
        n = n
    end

    @variables begin
        dm_a(t), [guess = 0]
        dm_b(t), [guess = 0]
    end

    @components begin
        port_a = HydraulicPort()
        port_b = HydraulicPort()
        open = Open()
    end

    @equations begin
        connect(port_a, port_b, open.port)
        dm_a ~ port_a.dm
        dm_b ~ dm_a / n
        open.dm ~ dm_a - dm_b # extra flow dumps into an open port
        # port_b.dm ~ dm_b # divided flow goes to port_b
    end
end

@component function ValveBase(
        reversible = false; minimum_area = 0, Cd, Cd_reverse = Cd, name)
    pars = @parameters begin
        Cd = Cd
        Cd_reverse = Cd_reverse
        minimum_area = minimum_area
    end

    systems = @named begin
        port_a = HydraulicPort()
        port_b = HydraulicPort()
    end

    vars = @variables begin
        area(t), [guess = 0]
        y(t), [guess = 0]
    end

    # let
    # Opting for a more numerically stable constant density (use head factor to compensate if needed)
    ρ = density_ref(port_a) #(full_density(port_a) + full_density(port_b)) / 2

    x = if reversible
        area
    else
        ifelse(area > minimum_area, area, minimum_area)
    end

    # let ------
    Δp = port_a.p - port_b.p
    dm = port_a.dm
    c = if reversible
        Cd
    else
        ifelse(Δp > 0, Cd, Cd_reverse)
    end

    eqs = [0 ~ port_a.dm + port_b.dm
           domain_connect(port_a, port_b)
           dm ~ regRoot(2 * Δp * ρ / c) * x # I think this should be reformulated as: regRoot(2 DP rho) c x
           y ~ x]

    System(eqs, t, vars, pars; name, systems)
end

"""
    Valve(reversible = false; p_a_int, p_b_int, area_int, Cd, Cd_reverse = Cd, minimum_area = 0, name)

Valve with `area` input and discharge coefficient `Cd` defined by https://en.wikipedia.org/wiki/Discharge_coefficient.  The `Cd_reverse` parameter allows for directional flow restriction, making it possible to define a check valve.

# Parameters:
- `p_a_int`: [Pa] initial pressure for `port_a`
- `p_b_int`: [Pa] initial pressure for `port_b`
- `area_int`: [m^2] initial valve opening
- `Cd`: discharge coefficient flowing from `a → b`
- `Cd_reverse`: discharge coefficient flowing from `b → a`
- `minimum_area`: when `reversible = false` applies a forced minimum area

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
- `area`: real input setting the valve `area`.  When `reversible = true`, negative input reverses flow direction, otherwise a floor of `minimum_area` is enforced.
"""
@component function Valve(reversible = false;
        Cd, Cd_reverse = Cd,
        minimum_area = 0,
        name)
    pars = @parameters begin
        Cd = Cd
        Cd_reverse = Cd_reverse
        minimum_area = minimum_area
    end

    systems = @named begin
        port_a = HydraulicPort()
        port_b = HydraulicPort()
        area = RealInput()
        base = ValveBase(reversible; Cd, Cd_reverse,
            minimum_area)
    end

    vars = []

    eqs = [connect(base.port_a, port_a)
           connect(base.port_b, port_b)
           base.area ~ area.u]

    System(eqs, t, vars, pars; name, systems)
end

@component function VolumeBase(; area, dead_volume = 0, p_int, x_int,
        name)
    pars = @parameters begin
        area = area
        dead_volume = dead_volume
        p_int = p_int
    end

    systems = @named begin
        port = HydraulicPort()
    end

    vars = @variables begin
        x(t) = x_int
        dx(t), [guess = 0]
        rho(t), [guess = liquid_density(port)]
        m(t), [guess = 0]
        vol(t)
    end

    # let
    dm = port.dm
    p = port.p

    eqs = [vol ~ dead_volume + area * x
           D(x) ~ dx
           D(m) ~ dm
           #rho ~ full_density(port, p)
           p ~ full_pressure(port, rho) # see https://github.com/SciML/OrdinaryDiffEq.jl/issues/2561
           m ~ rho * vol]

    initialization_eqs = [p ~ p_int]

    System(eqs, t, vars, pars; name, systems, initialization_eqs)
end

"""
    FixedVolume(; p_int, vol, name)

Fixed fluid volume.

# Parameters:
- `p_int`: [Pa] initial pressure
- `vol`: [m^3] fixed volume

# Connectors:
- `port`: hydraulic port
"""
@component function FixedVolume(; vol, name, p_int)
    pars = @parameters begin
        vol = vol
        p_int = p_int
    end

    systems = @named begin
        port = HydraulicPort(;)
    end

    vars = @variables begin
        rho(t), [guess = liquid_density(port)]
        m(t), [guess = vol * liquid_density(port)]
        p(t) = p_int
    end

    # let
    dm = port.dm

    eqs = [D(m) ~ dm
           #    rho ~ full_density(port, p)
           p ~ full_pressure(port, rho) # see https://github.com/SciML/OrdinaryDiffEq.jl/issues/2561
           p ~ port.p
           m ~ rho * vol]

    System(eqs, t, vars, pars; name, systems)
end

"""
    Volume(; x, dx=0, p, drho=0, dm=0, area, direction = 1, name)

Volume with moving wall with `flange` connector for converting hydraulic energy to 1D mechanical.  The `direction` argument aligns the mechanical port with the hydraulic port, useful when connecting two dynamic volumes together in oppsing directions to create an actuator.

```
     ┌─────────────────┐ ───
     │                 │  ▲
                       │  │
dm ────►               │  │ area
                       │  │
     │                 │  ▼
     └─────────────────┤ ───
                       │
                       └─► x (= ∫ flange.v * direction)
```

# Features:
- volume discretization with flow resistance and inertia: use `N` to control number of volume and resistance elements.  Set `N=0` to turn off volume discretization. See `TubeBase` for more information about flow resistance.
- minimum volume flow shutoff with damping and directional resistance.  Use `reversible=false` when problem defines volume position `x` and solves for `dm` to prevent numerical instability.

# Parameters:
## volume
- `p`: [Pa] initial pressure
- `area`: [m^2] moving wall area
- `x`: [m] initial wall position
- `dx=0`: [m/s] initial wall velocity
- `drho=0`: [kg/m^3/s] initial density derivative
- `dm=0`: [kg/s] initial flow

- `direction`: [+/-1] applies the direction conversion from the `flange` to `x`

# Connectors:
- `port`: hydraulic port
- `flange`: mechanical translational port

See also [`FixedVolume`](@ref), [`DynamicVolume`](@ref)
"""
@component function Volume(;
        #parameters
        area,
        direction = +1,
        x_int,
        name)
    pars = @parameters begin
        area = area
        x_int = x_int
    end

    vars = @variables begin
        x(t) = x_int
        dx(t), [guess = 0]
        p(t), [guess = 0]
        f(t), [guess = 0]
        rho(t), [guess = 0]
        m(t), [guess = 0]
        dm(t), [guess = 0]
    end

    systems = @named begin
        port = HydraulicPort()
        flange = MechanicalPort()
        damper = ValveBase(reversible;
            Cd,
            Cd_reverse,
            minimum_area)
    end

    systems = @named begin
        port = HydraulicPort()
        flange = MechanicalPort()
        damper = ValveBase(reversible;
            Cd,
            Cd_reverse,
            minimum_area)
    end

    eqs = [
           # connectors
           port.p ~ p
           port.dm ~ dm
           flange.v * direction ~ dx
           flange.f * direction ~ -f
           # differentials
           D(x) ~ dx
           D(m) ~ dm

           # physics
           # rho ~ full_density(port, p)
           p ~ full_pressure(port, rho) # see https://github.com/SciML/OrdinaryDiffEq.jl/issues/2561
           f ~ p * area
           m ~ rho * x * area]

    System(eqs, t, vars, pars; name, systems)
end

"""
    DynamicVolume(reversible = false; p_int,  area, x_int = 0, x_max, x_min = 0, x_damp = x_min, direction = +1, perimeter = 2 * sqrt(area * pi), shape_factor = 64, head_factor = 1, Cd = 1e2, Cd_reverse = Cd, name)

Volume with moving wall with `flange` connector for converting hydraulic energy to 1D mechanical.  The `direction` argument aligns the mechanical port with the hydraulic port, useful when connecting two dynamic volumes together in oppsing directions to create an actuator.

```
     ┌─────────────────┐ ───
     │                 │  ▲
                       │  │
dm ────►               │  │ area
                       │  │
     │                 │  ▼
     └─────────────────┤ ───
                       │
                       └─► x (= ∫ flange.v * direction)
```

# Features:
- minimum volume flow shutoff with damping and directional resistance.  Use `reversible=false` when problem defines volume position `x` and solves for `dm` to prevent numerical instability.

# Parameters:
## volume
- `p_int`: [Pa] initial pressure
- `area`: [m^2] moving wall area
- `x_max`: [m] max wall position, needed for volume discretization to apply the correct volume sizing as a function of `x`
- `x_min`: [m] wall position that shuts off flow and prevents negative volume.
- `x_damp`: [m] wall position that initiates a linear damping region before reaching full flow shut off.  Helps provide a smooth end stop.

- `direction`: [+/-1] applies the direction conversion from the `flange` to `x`

## flow resistance
- `perimeter`: [m] perimeter of the cross section (needed only for non-circular volumes)
- `shape_factor`: shape factor, see `friction_factor` function
- `head_factor`: effective length multiplier, used to account for addition friction from flow development and additional friction such as pipe bends, entrance/exit lossses, etc.

## flow shut off and damping
- `Cd`: discharge coefficient for flow out of the volume.  *Note: area is 1m² when valve is fully open.  Ensure this does not induce unwanted flow resistance.*
- `Cd_reverse`: discharge coefficient for flow into the volume. Use a lower value to allow easy wall release, in some cases the wall can "stick".


# Connectors:
- `port`: hydraulic port
- `flange`: mechanical translational port
"""
@component function DynamicVolume(reversible = false;
        area,
        x_int = 0,
        x_max,
        x_min = 0,
        x_damp = x_min,
        direction = +1,

        # Tube
        perimeter = 2 * sqrt(area * pi),
        shape_factor = 64,
        head_factor = 1,
        p_int,

        # Valve
        Cd = 1e2,
        Cd_reverse = Cd,
        minimum_area = 0,

        # Damping
        d = 0, name)
    @assert (direction == +1)||(direction == -1) "direction argument must be +/-1, found $direction"

    #TODO: How to set an assert effective_length >= length ??
    pars = @parameters begin
        area = area

        x_int = x_int
        x_max = x_max
        x_min = x_min
        x_damp = x_damp

        perimeter = perimeter
        shape_factor = shape_factor
        head_factor = head_factor
        p_int = p_int

        Cd = Cd
        Cd_reverse = Cd_reverse
        minimum_area = minimum_area

        d = d
    end

    vars = @variables begin
        x(t) = x_int
        vol(t), [guess = x_int * area]
    end

    systems = @named begin
        port = HydraulicPort(;)
        flange = MechanicalPort(;)
        damper = ValveBase(reversible;
            Cd,
            Cd_reverse,
            minimum_area)
        moving_volume = VolumeBase(;
            area,
            dead_volume = area * x_int,
            p_int,
            x_int = 0)
    end

    ratio = (x - x_min) / (x_damp - x_min)

    damper_area = if reversible
        one(x)
    else
        ifelse(x >= x_damp, one(x), ifelse((x < x_damp) & (x > x_min), ratio, zero(x)))
    end

    dx = moving_volume.dx
    p = moving_volume.port.p

    eqs = [vol ~ x * area
           D(x) ~ flange.v * direction
           damper.area ~ damper_area
           connect(port, damper.port_b)
           connect(moving_volume.port, damper.port_a)
           dx ~ flange.v * direction
           p * area - dx * d ~ -flange.f * direction]

    return System(eqs, t, vars, pars; name, systems)
end

"""
    SpoolValve(reversible = false; x_int, Cd, d, name)

Spool valve with `x` valve opening input as mechanical flange port and `d` diameter of orifice. See `Valve` for more information.

# Parameters:
- `x_int`: [m] initial valve opening
- `d`: [m] orifice diameter
- `Cd`: discharge coefficient flowing from `a → b`

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
- `flange`: mechanical translational port

See [`Valve`](@ref) for more information.
"""
@component function SpoolValve(reversible = false; Cd, d, x_int, name)
    pars = @parameters begin
        d = d
        Cd = Cd
        x_int = x_int
    end

    systems = @named begin
        port_a = HydraulicPort(;)
        port_b = HydraulicPort(;)
        flange = MechanicalPort()
        valve = ValveBase(reversible; Cd)
    end

    vars = @variables begin
        x(t) = x_int
        dx(t), [guess = 0]
    end

    eqs = [D(x) ~ dx
           flange.v ~ dx
           flange.f ~ 0 #TODO: model flow force
           connect(valve.port_a, port_a)
           connect(valve.port_b, port_b)
           valve.area ~ x * 2π * d]

    System(eqs, t, vars, pars; name, systems)
end

"""
    SpoolValve2Way(reversible = false; m, g, x_int, Cd, d, name)

2-ways spool valve with 4 ports and spool mass. Fluid flow direction S → A and B → R when `x` is positive and S → B and A → R when `x` is negative.

# Parameters:
- `m`: [kg] mass of the  spool
- `g`: [m/s²] gravity field acting on the spool, positive value acts in the positive direction
- `x_int`: [m] initial valve opening
- `d`: [m] orifice diameter
- `Cd`: discharge coefficient flowing from `s → a` and `b → r`

# Connectors:
- `port_s`: hydraulic port
- `port_a`: hydraulic port
- `port_b`: hydraulic port
- `port_r`: hydraulic port
- `flange`: mechanical translational port

See [`SpoolValve`](@ref) for more information.
"""
@component function SpoolValve2Way(reversible = false; m, g, Cd, d, x_int, name)
    pars = @parameters begin
        m = m
        g = g

        d = d

        Cd = Cd

        x_int = x_int
        # dx_int = dx_int
    end

    vars = []

    systems = @named begin
        vSA = SpoolValve(reversible; Cd, d, x_int)
        vBR = SpoolValve(reversible; Cd, d, x_int)

        port_s = HydraulicPort(;)
        port_a = HydraulicPort(;)
        port_b = HydraulicPort(;)
        port_r = HydraulicPort(;)

        mass = Mass(; m = m, g = g)

        flange = MechanicalPort()
    end

    eqs = [connect(vSA.port_a, port_s)
           connect(vSA.port_b, port_a)
           connect(vBR.port_a, port_b)
           connect(vBR.port_b, port_r)
           connect(vSA.flange, vBR.flange, mass.flange, flange)]

    initialization_eqs = [
        mass.s ~ x_int
    # mass.v ~ dx_int
    ]

    System(eqs, t, vars, pars; name, systems, initialization_eqs)
end

"""
    Actuator(N, add_inertia = true, reversible = false;
        p_a_int,
        p_b_int,
        area_a,
        area_b,
        perimeter_a = 2 * sqrt(area_a * pi),
        perimeter_b = 2 * sqrt(area_b * pi),
        length_a_int,
        length_b_int,
        shape_factor_a = 64,
        shape_factor_b = 64,
        head_factor_a = 1,
        head_factor_b = 1,
        m,
        g,
        x_int = 0,
        minimum_volume_a = 0,
        minimum_volume_b = 0,
        damping_volume_a = minimum_volume_a,
        damping_volume_b = minimum_volume_b,
        Cd = 1e4,
        Cd_reverse = Cd,
        p_a_int,
        p_b_int,
        name)

Actuator made of two DynamicVolumes connected in opposite direction with body mass attached.

# Features:
- volume discretization with flow resistance and inertia: use `N` to control number of volume and resistance elements.  Set `N=0` to turn off volume discretization. See `TubeBase` for more information about flow resistance.
- minimum volume flow shutoff with damping and directional resistance.  Use `reversible=false` when problem defines volume position `x` and solves for `dm` to prevent numerical instability.

# Parameters:
## volume
- `p_a_int`: [Pa] initial pressure for `port_a`
- `p_b_int`: [Pa] initial pressure for `port_b`
- `area_a`: [m^2] moving wall area of volume `A`
- `area_b`: [m^2] moving wall area of volume `B`
- `length_a_int`: [m] initial wall position for `A`
- `length_b_int`: [m] initial wall position for `b`

## mass
- `m`: [kg] mass of the body
- `g`: [m/s²] gravity field acting on the mass, positive value acts in the positive direction
- `x_int`: [m] initial flange position

## flow resistance
- `perimeter_a`: [m] perimeter of the cross section `A` (needed only for non-circular volumes)
- `perimeter_b`: [m] perimeter of the cross section `B` (needed only for non-circular volumes)
- `shape_factor_a`: shape factor of `A`, see `friction_factor` function
- `shape_factor_b`: shape factor of `B`, see `friction_factor` function
- `head_factor_a`: effective length multiplier for `A`, used to account for addition friction from flow development and additional friction such as pipe bends, entrance/exit lossses, etc.
- `head_factor_b`: effective length multiplier for `B`, used to account for addition friction from flow development and additional friction such as pipe bends, entrance/exit lossses, etc.

## flow shut off and damping
- `minimum_volume_a`: [m^3] minimum volume `A` that shuts off flow and prevents negative volume.
- `minimum_volume_b`: [m^3] minimum volume `B` that shuts off flow and prevents negative volume.
- `damping_volume_a`: [m^3] volume of `A` that initiates a linear damping region before reaching full flow shut off.  Helps provide a smooth end stop.
- `damping_volume_b`: [m^3] volume of `B` that initiates a linear damping region before reaching full flow shut off.  Helps provide a smooth end stop.
- `Cd`: discharge coefficient for flow out of the volume.  *Note: area is 1m² when valve is fully open.  Ensure this does not induce unwanted flow resistance.*
- `Cd_reverse`: discharge coefficient for flow into the volume. Use a lower value to allow easy wall release, in some cases the wall can "stick".


# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
- `flange`: mechanical translational port
"""
@component function Actuator(reversible = false;
        area_a,
        area_b,
        perimeter_a = 2 * sqrt(area_a * pi),
        perimeter_b = 2 * sqrt(area_b * pi),
        length_a_int,
        length_b_int,
        shape_factor_a = 64,
        shape_factor_b = 64,
        head_factor_a = 1,
        head_factor_b = 1,
        m,
        g,
        x_int = 0,
        dx_int = 0,
        minimum_volume_a = 0,
        minimum_volume_b = 0,
        damping_volume_a = minimum_volume_a,
        damping_volume_b = minimum_volume_b,
        Cd = 1e4,
        Cd_reverse = Cd,
        d = 0,
        p_a_int,
        p_b_int,
        name)
    pars = @parameters begin
        area_a = area_a
        area_b = area_b
        perimeter_a = perimeter_a
        perimeter_b = perimeter_b
        shape_factor_a = shape_factor_a
        shape_factor_b = shape_factor_b
        head_factor_a = head_factor_a
        head_factor_b = head_factor_b
        x_int = x_int
        dx_int = dx_int
        length_a_int = length_a_int
        length_b_int = length_b_int
        minimum_volume_a = minimum_volume_a
        minimum_volume_b = minimum_volume_b
        damping_volume_a = damping_volume_a
        damping_volume_b = damping_volume_b
        Cd = Cd
        Cd_reverse = Cd_reverse
        m = m
        g = g
        d = d
        p_a_int = p_a_int
        p_b_int = p_b_int
    end

    vars = @variables begin
        x(t) = x_int
        dx(t) = dx_int
    end

    total_length = length_a_int + length_b_int

    #TODO: include effective_length
    systems = @named begin
        vol_a = DynamicVolume(reversible; direction = +1,
            area = area_a,
            x_int = length_a_int,
            x_max = total_length,
            x_min = minimum_volume_a / area_a,
            x_damp = damping_volume_a / area_a,
            perimeter = perimeter_a,
            shape_factor = shape_factor_a,
            head_factor = head_factor_a,
            Cd,
            Cd_reverse,
            d,
            p_int = p_a_int)

        vol_b = DynamicVolume(reversible; direction = -1,
            area = area_b,
            x_int = length_b_int,
            x_max = total_length,
            x_min = minimum_volume_b / area_b,
            x_damp = damping_volume_b / area_b,
            perimeter = perimeter_b,
            shape_factor = shape_factor_b,
            head_factor = head_factor_b,
            Cd,
            Cd_reverse,
            d,
            p_int = p_b_int)
        mass = Mass(; m, g)
        port_a = HydraulicPort()
        port_b = HydraulicPort()
        flange = MechanicalPort()
    end

    eqs = [connect(vol_a.port, port_a)
           connect(vol_b.port, port_b)
           connect(vol_a.flange, vol_b.flange, mass.flange, flange)
           D(x) ~ dx
           dx ~ vol_a.flange.v]

    initialization_eqs = [
        mass.s ~ x_int
    ]

    System(eqs, t, vars, pars; name, systems, initialization_eqs)
end

"""
    Orifice()

A valve in fixed position, with parameters for area and the discharge coefficient (fitting the form Effective Area = area x Cd)  

```
     ┌ 
     │                   
         ▲
dm ────►  effective area
         ▼
     │           
     └
```

# Features:
- 

# Parameters:
## volume
- `area`: [m^2] physical area
- `cd`: [unitless] discharge coefficient

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
"""

@mtkmodel Orifice begin
    @parameters begin
        orifice_area = 0.00094
        Cd = 0.6 # TODO Cd here is defined differently from Valve(). 
        # Here it follows the form Effective Orifice Area = Cd x Physical Orifice Area
        # The Valve component should be updated too.
    end
    @components begin
        area = Constant(k = orifice_area)
        valve = Valve(Cd = 1 / (Cd * Cd))
        port_a = HydraulicPort()
        port_b = HydraulicPort()
    end
    @equations begin
        connect(valve.area, area.output)
        connect(valve.port_a, port_a)
        connect(valve.port_b, port_b)
    end
end
