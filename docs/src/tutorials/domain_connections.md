# [Domains](@id domains)
## Basics
A domain in ModelingToolkit.jl is a network of connected components that share properties of the medium in the network.  For example, a collection of hydraulic components connected together will have a fluid medium.  Using the domain feature, one only needs to define and set the fluid medium properties once, in one component, rather than at each component.  The way this works in ModelingToolkit.jl is by defining a connector (with Through/Flow and Across variables) with parameters defining the medium of the domain.  Then a second connector is defined, with the same parameters, and the same Through/Flow variable, which acts as the setter.  For example, a hydraulic domain may have a hydraulic connector, `HydraulicPort`, that defines a fluid medium with density (`ρ`), viscosity (`μ`), and a bulk modulus (`β`), a through/flow variable mass flow (`dm`) and an across variable pressure (`p`).

```@example domain
using ModelingToolkit

@parameters t
D = Differential(t)

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
```

The fluid medium setter for `HydralicPort` may be defined as `HydraulicFluid` with the same parameters and through/flow variable.  But now, the parameters can be set through the function keywords.

```@example domain
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
        dm ~ 0,
    ]

    ODESystem(eqs, t, vars, pars; name, defaults = [dm => 0])
end
```

Now, we can connect a `HydraulicFluid` component to any `HydraulicPort` connector, and the parameters of all `HydraulicPort`'s in the network will be automatically set.  Let's consider a simple example, connecting a pressure source component to a volume component.  Note that we don't need to define density for the volume component, it's supplied by the `HydraulicPort` (`port.ρ`).

```@example domain
@component function FixedPressure(; p, name)
    pars = @parameters p = p
    systems = @named begin port = HydraulicPort(; p_int = p) end

    eqs = [port.p ~ p]

    ODESystem(eqs, t, [], pars; name, systems)
end

@component function FixedVolume(; vol, p_int, name)
    pars = @parameters begin
        p_int = p_int
        vol = vol
    end

    systems = @named begin port = HydraulicPort(; p_int) end

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
```
When the system is defined we can generate a fluid component and connect it to the system.  Here `fluid` is connected to the `src.port`, but it could also be connected to `vol.port`, any connection in the network is fine.  Note: we can visualize the system using `ModelingToolkitDesigner.jl`, where a dashed line is used to show the `fluid` connection to represent a domain connection that is only transporting parameters and not states.

```@example domain
@component function System(; name)
    systems = @named begin
        src = FixedPressure(; p=200e5)
        vol = FixedVolume(; vol=0.1, p_int=200e5)

        fluid = HydraulicFluid(; density=876)
    end

    eqs = [
        connect(fluid, src.port)
        connect(src.port, vol.port)
    ]

    ODESystem(eqs, t, [], []; systems, name)
end

@named odesys = System()

using ModelingToolkitDesigner
path = joinpath(@__DIR__, "domain_connections") 
design = ODESystemDesign(odesys, path);
ModelingToolkitDesigner.view(design, false)
```

To see how the domain works, we can examine the set parameter values for each of the ports `src.port` and `vol.port`.  First we assemble the system using `structural_simplify()` and then check the default value of `vol.port.ρ`, whichs points to the setter value `fluid₊ρ`.  Likewise, `src.port.ρ`, will also point to the setter value `fluid₊ρ`.  Therefore, there is now only 1 defined density value `fluid₊ρ` which sets the density for the connected network.  

```@example domain
sys = structural_simplify(odesys)
ModelingToolkit.defaults(sys)[complete(odesys).vol.port.ρ]
```

If we have a more complicated system, for example a hydraulic actuator, with a separated fluid on both sides of the piston, it's possible we might have 2 separate domain networks.  In this case we can connect 2 separate fluids, or the same fluid, to both networks.  

```@example domain
@component function Actuator(; p_int, mass, area, name)
    pars = @parameters begin
        p_int = p_int
        mass = mass
        area = area
    end

    systems = @named begin 
        port_a = HydraulicPort(; p_int) 
        port_b = HydraulicPort(; p_int) 
    end

    vars = @variables begin
        x(t) = 0
        dx(t) = 0
        ddx(t) = 0
    end

    eqs = [
        D(x) ~ dx
        D(dx) ~ ddx
        mass*ddx ~ (port_a.p - port_b.p)*area
        port_a.dm ~ +(port_a.ρ)*dx*area
        port_b.dm ~ -(port_b.ρ)*dx*area
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

@component function ActuatorSystem2(; name)
    systems = @named begin
        src_a = FixedPressure(; p=200e5)
        src_b = FixedPressure(; p=200e5)
        act = Actuator(; p_int=200e5, mass=1000, area=0.1)

        fluid_a = HydraulicFluid(; density=876)
        fluid_b = HydraulicFluid(; density=999)
    end

    eqs = [
        connect(fluid_a, src_a.port)
        connect(fluid_b, src_b.port)
        connect(src_a.port, act.port_a)
        connect(src_b.port, act.port_b)
    ]

    ODESystem(eqs, t, [], []; systems, name)
end

@named actsys2 = ActuatorSystem2()

design2 = ODESystemDesign(actsys2, path);
ModelingToolkitDesigner.view(design2, false)


@component function ActuatorSystem1(; name)
    systems = @named begin
        src_a = FixedPressure(; p=200e5)
        src_b = FixedPressure(; p=200e5)
        act = Actuator(; p_int=200e5, mass=1000, area=0.1)

        fluid = HydraulicFluid(; density=876)
    end

    eqs = [
        connect(fluid, src_a.port)
        connect(fluid, src_b.port)

        connect(src_a.port, act.port_a)
        connect(src_b.port, act.port_b)
    ]

    ODESystem(eqs, t, [], []; systems, name)
end

@named actsys1 = ActuatorSystem1()

design1 = ODESystemDesign(actsys1, path);
ModelingToolkitDesigner.view(design1, false)
```

