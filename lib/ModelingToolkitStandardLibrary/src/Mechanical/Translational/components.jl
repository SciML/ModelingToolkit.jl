"""
    Free(; name)

Use to close a system that has un-connected `MechanicalPort`'s where the force should not be zero (i.e. you want to solve for the force to produce the given movement of the port)

# Connectors:

  - `flange`: 1-dim. translational flange
"""
@mtkmodel Free begin
    @components begin
        flange = MechanicalPort()
    end
    @variables begin
        f(t)
    end
    @equations begin
        flange.f ~ f
    end
end

"""
    Fixed(; name)

Fixes a flange position (velocity = 0)

# Connectors:

  - `flange`: 1-dim. translational flange
"""
@mtkmodel Fixed begin
    @components begin
        flange = MechanicalPort()
    end
    @equations begin
        flange.v ~ 0
    end
end

"""
    Mass(; name, m, g = 0)

Sliding mass with inertia

# Parameters:

  - `m`: [kg] mass of sliding body
  - `g = 0`: [m/s^2] [m/s²] gravity field acting on the mass, positive value acts in the positive direction
  

# States:

  - `v`: [m/s] absolute linear velocity of sliding mass
  - `s`: [m] absolute position of sliding mass (optional with parameter s)

# Connectors:

  - `flange`: 1-dim. translational flange
"""
@component function Mass(; name, m, g = 0)
    pars = @parameters begin
        m = m
        g = g
    end
    @named flange = MechanicalPort()

    vars = @variables begin
        s(t), [guess = 0]
        v(t), [guess = 0]
        f(t), [guess = 0]
    end

    eqs = [flange.v ~ v
           flange.f ~ f
           D(s) ~ v
           D(v) ~ f / m + g]

    return compose(System(eqs, t, vars, pars; name = name),
        flange)
end

const REL = Val(:relative)

"""
    Spring(; name, k, delta_s = 0.0,  va=0.0, v_b_0=0.0)

Linear 1D translational spring

# Parameters:

  - `k`: [N/m] Spring constant
  - `delta_s`: initial spring stretch
  - `va`: [m/s] Initial value of absolute linear velocity at flange_a (default 0 m/s)
  - `v_b_0`: [m/s] Initial value of absolute linear velocity at flange_b (default 0 m/s)

# Connectors:

  - `flange_a`: 1-dim. translational flange on one side of spring
  - `flange_b`: 1-dim. translational flange on opposite side of spring
"""
@component function Spring(; name, k)
    Spring(REL; name, k)
end # default

@component function Spring(::Val{:relative}; name, k)
    pars = @parameters begin
        k = k
    end
    vars = @variables begin
        delta_s(t), [guess = 0]
        f(t), [guess = 0]
    end

    @named flange_a = MechanicalPort()
    @named flange_b = MechanicalPort()

    eqs = [D(delta_s) ~ flange_a.v - flange_b.v
           f ~ k * delta_s
           flange_a.f ~ +f
           flange_b.f ~ -f]
    return compose(System(eqs, t, vars, pars; name = name),
        flange_a,
        flange_b) #flange_a.f => +k*delta_s, flange_b.f => -k*delta_s
end

const ABS = Val(:absolute)
@component function Spring(::Val{:absolute}; name, k, l = 0)
    pars = @parameters begin
        k = k
        l = l
    end
    vars = @variables begin
        sa(t), [guess = 0]
        sb(t), [guess = 0]
        f(t), [guess = 0]
    end

    @named flange_a = MechanicalPort()
    @named flange_b = MechanicalPort()

    eqs = [D(sa) ~ flange_a.v
           D(sb) ~ flange_b.v
           f ~ k * (sa - sb - l) #delta_s
           flange_a.f ~ +f
           flange_b.f ~ -f]
    return compose(System(eqs, t, vars, pars; name = name),
        flange_a,
        flange_b) #, flange_a.f => k * (flange_a__s - flange_b__s - l)
end

"""
    Damper(; name, d, flange_a.v = 0.0, flange_b.v = 0.0)

Linear 1D translational damper

# Parameters:

  - `d`: [N.s/m] Damping constant

# Connectors:

  - `flange_a`: 1-dim. translational flange on one side of damper. Initial value of state `v` is set to 0.0 m/s.
  - `flange_b`: 1-dim. translational flange on opposite side of damper. Initial value of state `v` is set to 0.0 m/s.
"""
@mtkmodel Damper begin
    @parameters begin
        d
    end
    @variables begin
        v(t), [guess = 0]
        f(t), [guess = 0]
    end

    @components begin
        flange_a = MechanicalPort()
        flange_b = MechanicalPort()
    end

    @equations begin
        v ~ flange_a.v - flange_b.v
        f ~ v * d
        flange_a.f ~ +f
        flange_b.f ~ -f
    end
end
