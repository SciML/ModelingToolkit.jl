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
nothing #hide
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
nothing #hide
```

Now, we can connect a `HydraulicFluid` component to any `HydraulicPort` connector, and the parameters of all `HydraulicPort`'s in the network will be automatically set.  Let's consider a simple example, connecting a pressure source component to a volume component.  Note that we don't need to define density for the volume component, it's supplied by the `HydraulicPort` (`port.ρ`).

```@example domain
@component function FixedPressure(; p, name)
    pars = @parameters p = p
    systems = @named begin
        port = HydraulicPort(; p_int = p)
    end

    eqs = [port.p ~ p]

    ODESystem(eqs, t, [], pars; name, systems)
end

@component function FixedVolume(; vol, p_int, name)
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
nothing #hide
```

When the system is defined we can generate a fluid component and connect it to the system.  Here `fluid` is connected to the `src.port`, but it could also be connected to `vol.port`, any connection in the network is fine.  Note: we can visualize the system using `ModelingToolkitDesigner.jl`, where a dashed line is used to show the `fluid` connection to represent a domain connection that is only transporting parameters and not states.

```@example domain
@component function System(; name)
    systems = @named begin
        src = FixedPressure(; p = 200e5)
        vol = FixedVolume(; vol = 0.1, p_int = 200e5)

        fluid = HydraulicFluid(; density = 876)
    end

    eqs = [connect(fluid, src.port)
        connect(src.port, vol.port)]

    ODESystem(eqs, t, [], []; systems, name)
end

@named odesys = System()
nothing #hide
```

```@setup domain
# code to generate diagrams...
# using ModelingToolkitDesigner 
# path = raw"C:\Work\Assets\ModelingToolkit.jl\domain_connections" 
# design = ODESystemDesign(odesys, path); 

# using CairoMakie
# CairoMakie.set_theme!(Theme(;fontsize=12))
# fig = ModelingToolkitDesigner.view(design, false)
# save(joinpath(path, "odesys.svg"), fig; resolution=(300,300))
```

![odesys](https://github.com/SciML/ModelingToolkit.jl/assets/40798837/d19fbcf4-781c-4743-87b7-30bed348ff98)

To see how the domain works, we can examine the set parameter values for each of the ports `src.port` and `vol.port`.  First we assemble the system using `structural_simplify()` and then check the default value of `vol.port.ρ`, whichs points to the setter value `fluid₊ρ`.  Likewise, `src.port.ρ`, will also point to the setter value `fluid₊ρ`.  Therefore, there is now only 1 defined density value `fluid₊ρ` which sets the density for the connected network.

```@repl domain
sys = structural_simplify(odesys)
ModelingToolkit.defaults(sys)[complete(odesys).vol.port.ρ]
```

## Multiple Domain Networks

If we have a more complicated system, for example a hydraulic actuator, with a separated fluid on both sides of the piston, it's possible we might have 2 separate domain networks.  In this case we can connect 2 separate fluids, or the same fluid, to both networks.  First a simple actuator is defined with 2 ports.

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

    eqs = [D(x) ~ dx
        D(dx) ~ ddx
        mass * ddx ~ (port_a.p - port_b.p) * area
        port_a.dm ~ +(port_a.ρ) * dx * area
        port_b.dm ~ -(port_b.ρ) * dx * area]

    ODESystem(eqs, t, vars, pars; name, systems)
end
nothing #hide
```

A system with 2 different fluids is defined and connected to each separate domain network.

```@example domain
@component function ActuatorSystem2(; name)
    systems = @named begin
        src_a = FixedPressure(; p = 200e5)
        src_b = FixedPressure(; p = 200e5)
        act = Actuator(; p_int = 200e5, mass = 1000, area = 0.1)

        fluid_a = HydraulicFluid(; density = 876)
        fluid_b = HydraulicFluid(; density = 999)
    end

    eqs = [connect(fluid_a, src_a.port)
        connect(fluid_b, src_b.port)
        connect(src_a.port, act.port_a)
        connect(src_b.port, act.port_b)]

    ODESystem(eqs, t, [], []; systems, name)
end

@named actsys2 = ActuatorSystem2()
nothing #hide
```

```@setup domain
# design = ODESystemDesign(actsys2, path);
# fig = ModelingToolkitDesigner.view(design, false)
# save(joinpath(path, "actsys2.svg"), fig; resolution=(500,300))
```

![actsys2](https://github.com/SciML/ModelingToolkit.jl/assets/40798837/8ed50035-f6ac-48cb-a585-1ef415154a02)

After running `structural_simplify()` on `actsys2`, the defaults will show that `act.port_a.ρ` points to `fluid_a₊ρ` and `act.port_b.ρ` points to `fluid_b₊ρ`.  This is a special case, in most cases a hydraulic system will have only 1 fluid, however this simple system has 2 separate domain networks.  Therefore, we can connect a single fluid to both networks.  This does not interfere with the mathematical equations of the system, since no states are connected.

```@example domain
@component function ActuatorSystem1(; name)
    systems = @named begin
        src_a = FixedPressure(; p = 200e5)
        src_b = FixedPressure(; p = 200e5)
        act = Actuator(; p_int = 200e5, mass = 1000, area = 0.1)

        fluid = HydraulicFluid(; density = 876)
    end

    eqs = [connect(fluid, src_a.port)
        connect(fluid, src_b.port)
        connect(src_a.port, act.port_a)
        connect(src_b.port, act.port_b)]

    ODESystem(eqs, t, [], []; systems, name)
end

@named actsys1 = ActuatorSystem1()
nothing #hide
```

```@setup domain
# design = ODESystemDesign(actsys1, path);
# fig = ModelingToolkitDesigner.view(design, false)
# save(joinpath(path, "actsys1.svg"), fig; resolution=(500,300))
```

![actsys1](https://github.com/SciML/ModelingToolkit.jl/assets/40798837/054404eb-dbb7-4b85-95c0-c9503d0c4d00)

## Special Connection Cases (`domain_connect()`)

In some cases a component will be defined with 2 connectors of the same domain, but they are not connected.  For example the `Restrictor` defined here gives equations to define the behavior of how the 2 connectors `port_a` and `port_b` are physically connected.

```@example domain
@component function Restrictor(; name, p_int)
    pars = @parameters begin
        K = 0.1
        p_int = p_int
    end

    systems = @named begin
        port_a = HydraulicPort(; p_int)
        port_b = HydraulicPort(; p_int)
    end

    eqs = [port_a.dm ~ (port_a.p - port_b.p) * K
        0 ~ port_a.dm + port_b.dm]

    ODESystem(eqs, t, [], pars; systems, name)
end
nothing #hide
```

Adding the `Restrictor` to the original system example will cause a break in the domain network, since a `connect(port_a, port_b)` is not defined.

```@example domain
@component function RestrictorSystem(; name)
    systems = @named begin
        src = FixedPressure(; p = 200e5)
        res = Restrictor(; p_int = 200e5)
        vol = FixedVolume(; vol = 0.1, p_int = 200e5)

        fluid = HydraulicFluid(; density = 876)
    end

    eqs = [connect(fluid, src.port)
        connect(src.port, res.port_a)
        connect(res.port_b, vol.port)]

    ODESystem(eqs, t, [], []; systems, name)
end

@named ressys = RestrictorSystem()
sys = structural_simplify(ressys)
nothing #hide
```

```@setup domain
# design = ODESystemDesign(ressys, path);
# fig = ModelingToolkitDesigner.view(design, false)
# save(joinpath(path, "ressys.svg"), fig; resolution=(500,300))
```

![ressys](https://github.com/SciML/ModelingToolkit.jl/assets/40798837/3740f0e2-7324-4c1f-af8b-eba02cfece81)

When `structural_simplify()` is applied to this system it can be seen that the defaults are missing for `res.port_b` and `vol.port`.

```@repl domain
ModelingToolkit.defaults(sys)[complete(ressys).res.port_a.ρ]
ModelingToolkit.defaults(sys)[complete(ressys).res.port_b.ρ]
ModelingToolkit.defaults(sys)[complete(ressys).vol.port.ρ]
```

To ensure that the `Restrictor` component does not disrupt the domain network, the [`domain_connect()`](@ref) function can be used, which explicitly only connects the domain network and not the states.

```@example domain
@component function Restrictor(; name, p_int)
    pars = @parameters begin
        K = 0.1
        p_int = p_int
    end

    systems = @named begin
        port_a = HydraulicPort(; p_int)
        port_b = HydraulicPort(; p_int)
    end

    eqs = [domain_connect(port_a, port_b) # <-- connect the domain network
        port_a.dm ~ (port_a.p - port_b.p) * K
        0 ~ port_a.dm + port_b.dm]

    ODESystem(eqs, t, [], pars; systems, name)
end

@named ressys = RestrictorSystem()
sys = structural_simplify(ressys)
nothing #hide
```

Now that the `Restrictor` component is properly defined using `domain_connect()`, the defaults for `res.port_b` and `vol.port` are properly defined.

```@repl domain
ModelingToolkit.defaults(sys)[complete(ressys).res.port_a.ρ]
ModelingToolkit.defaults(sys)[complete(ressys).res.port_b.ρ]
ModelingToolkit.defaults(sys)[complete(ressys).vol.port.ρ]
```
