"""
    Force(; name)

Linear 1D force input source

# Connectors:

  - `flange`: 1-dim. translational flange
  - `f`: real input
"""
@mtkmodel Force begin
    @components begin
        flange = MechanicalPort(; v = 0.0)
        f = RealInput()
    end

    @equations begin
        flange.f ~ -f.u
    end
end

"""
    Position(solves_force = true; name)

Linear 1D position input source.  Set `solves_force=false` to force input force to 0 (i.e. only the position is given, the respective force needed is already provided elsewhere in the model).  

# Connectors:

  - `flange`: 1-dim. translational flange
  - `s`: real input
"""
@component function Position(solves_force = true; name)
    vars = []

    systems = @named begin
        flange = MechanicalPort(; v = 0)
        s = RealInput()
    end

    eqs = [
        D(s.u) ~ flange.v
    ]

    !solves_force && push!(eqs, 0 ~ flange.f)

    System(eqs, t, vars, [];
        name, systems)
end

"""
    Velocity(solves_force = true; name)

Linear 1D position input source.  Set `solves_force=false` to force input force to 0 (i.e. only the velocity is given, the respective force needed is already provided elsewhere in the model).  

# Connectors:

  - `flange`: 1-dim. translational flange
  - `v`: real input
"""
@component function Velocity(solves_force = true; name)
    systems = @named begin
        flange = MechanicalPort(; v = 0)
        v = RealInput()
    end

    eqs = [
        v.u ~ flange.v
    ]

    !solves_force && push!(eqs, 0 ~ flange.f)

    System(eqs, t, [], []; name, systems)
end

"""
Acceleration(solves_force = true; name)

Linear 1D position input source.  Set `solves_force=false` to force input force to 0 (i.e. only the acceleration is given, the respective force needed is already provided elsewhere in the model).  

# Connectors:

  - `flange`: 1-dim. translational flange
  - `a`: real input
"""
@component function Acceleration(solves_force = true; name)
    systems = @named begin
        flange = MechanicalPort(; v = 0)
        a = RealInput()
    end

    vars = @variables v(t)

    eqs = [v ~ flange.v
           D(v) ~ a.u]

    !solves_force && push!(eqs, 0 ~ flange.f)

    System(eqs, t, vars, []; name, systems)
end
