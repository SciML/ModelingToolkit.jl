@connector HeatPort begin
    @parameters begin
        T_guess = 273.15 + 20
        Q_flow_guess = 0.0
    end

    @variables begin
        T(t), [guess = T_guess]
        Q_flow(t), [guess = Q_flow_guess, connect = Flow]
    end
end
Base.@doc """
    HeatPort(; T = nothing, T_guess = 273.15 + 20, Q_flow = nothing, Q_flow_guess = 0.0, name)

Port for a thermal system.
# Parameters: 
- `T_guess`: [K] Initial guess for the temperature of the port (set to 273.15 + 20).
- `Q_flow_guess`: [W] Initial guess for the heat flow rate at the port (set to 0.0).

# States:
- `T`: [K] Temperature of the port. Guess set to `T_guess`. Passing a value for `T` will set its default.
- `Q_flow`: [W] Heat flow rate at the port. Guess set to `Q_flow_guess`. Passing a value for `Q_flow` will set its default.
""" HeatPort

"""
    Element1D(; name, dT = 0.0, Q_flow = 0.0)

This partial model contains the basic connectors and variables to allow heat transfer models to be created that do not
store energy. This model defines and includes equations for the temperature drop across the element, `dT`, and the heat
flow rate through the element from `port_a` to `port_b`, `Q_flow`.

# States:

  - `dT`:  [`K`] Temperature difference across the component a.T - b.T. It accepts an initial value, which defaults to 0.0.
  - `Q_flow`: [`W`] Heat flow rate from port a -> port b. It accepts an initial value, which defaults to 0.0.

# Connectors:

`port_a`
`port_b`
"""
@mtkmodel Element1D begin
    @components begin
        port_a = HeatPort()
        port_b = HeatPort()
    end
    @variables begin
        dT(t), [guess = 0.0]
        Q_flow(t), [guess = 0.0]
    end
    @equations begin
        dT ~ port_a.T - port_b.T
        port_a.Q_flow ~ Q_flow
        port_a.Q_flow + port_b.Q_flow ~ 0
    end
end

"""
    ConvectiveElement1D(; name, dT = 0.0, Q_flow = 0.0)

This partial model contains the basic connectors and variables to allow heat
transfer models to be created that do not store energy. This model defines and
includes equations for the temperature drop across the element, `dT`, and the heat
flow rate through the element from `solid` to `fluid`, `Q_flow`.

# States:

  - `dT`:  [`K`] Temperature difference across the component `solid.T` - `fluid.T`. It accepts an initial value, which defaults to 0.0.
  - `Q_flow`: [`W`] Heat flow rate from `solid` -> `fluid`. It accepts an initial value, which defaults to 0.0.

# Connectors:

`solid`
`fluid`
"""
@mtkmodel ConvectiveElement1D begin
    @components begin
        solid = HeatPort()
        fluid = HeatPort()
    end
    @variables begin
        dT(t), [guess = 0.0]
        Q_flow(t), [guess = 0.0]
    end
    @equations begin
        dT ~ solid.T - fluid.T
        solid.Q_flow ~ Q_flow
        solid.Q_flow + fluid.Q_flow ~ 0
    end
end
