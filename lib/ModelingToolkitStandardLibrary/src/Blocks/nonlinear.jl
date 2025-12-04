_clamp(u, u_min, u_max) = max(min(u, u_max), u_min)
_dead_zone(u, u_min, u_max) = ifelse(u > u_max, u - u_max, ifelse(u < u_min, u - u_min, 0))

"""
    Limiter(;name, y_max, y_min = y_max > 0 ? -y_max : -Inf)

Limit the range of a signal.

# Parameters:

  - `y_max`: Maximum of output signal
  - `y_min`: Minimum of output signal

# Connectors:

  - `input`
  - `output`
"""
@component function Limiter(; name, y_max, y_min = y_max > 0 ? -y_max : -Inf)
    @symcheck y_max ≥ y_min || throw(ArgumentError("`y_min` must be smaller than `y_max`"))
    m = (y_max + y_min) / 2
    siso = SISO(u_start = m, y_start = m, name = :siso) # Default signals to center of saturation to minimize risk of saturation while linearizing etc.
    @unpack u, y = siso
    pars = @parameters y_max=y_max [description="Maximum allowed output of Limiter $name"] y_min=y_min [
        description="Minimum allowed output of Limiter $name"
    ]
    eqs = [
        y ~ _clamp(u, y_min, y_max)
    ]
    extend(System(eqs, t, [], pars; name = name), siso)
end

"""
    DeadZone(; name, u_max, u_min = -u_max)

The DeadZone block defines a region of zero output.
If the input is within `u_min` ... `u_max`, the output is zero. Outside of this zone, the output is a linear function of the input with a slope of 1.

```
       y▲
        │     /
        │    /
  u_min │   /
─────|──┼──|───────► u
    /   │   u_max
   /    │
  /     │
```

# Parameters:

  - `u_max`: Upper limit of dead zone
  - `u_min`: Lower limit of dead zone

# Connectors:

  - `input`
  - `output`
"""
@mtkmodel DeadZone begin
    @parameters begin
        u_max, [description = "Upper limit of dead zone of DeadZone"]
        u_min = -u_max, [description = "Lower limit of dead zone of DeadZone"]
    end
    begin
        if !ModelingToolkit.isvariable(u_max)
            u_max ≥ u_min || throw(ArgumentError("`u_min` must be smaller than `u_max`"))
        end
    end

    @extend u, y = siso = SISO()

    @equations begin
        y ~ _dead_zone(u, u_min, u_max)
    end
end

"""
    SlewRateLimiter(; name, y_start, rising = 1.0, falling = -rising, Td = 0.001)

Limits the slew rate of a signal.
Initial value of state `Y` can be set with `int.y`

# Parameters:

  - `rising`: Maximum rising slew rate
  - `falling`: Maximum falling slew rate
  - `Td`: [s] Derivative time constant
  - `y_start`: Initial value of `y` state of SISO

# Connectors:

  - `input`
  - `output`
"""
@component function SlewRateLimiter(;
        name, y_start = 0.0, rising = 1.0, falling = -rising, Td = 0.001)
    pars = @parameters begin
        rising = rising, [description = "Maximum rising slew rate of SlewRateLimiter"]
        falling = falling, [description = "Derivative time constant of SlewRateLimiter"]
        Td = Td, [description = "Derivative time constant"]
        y_start = y_start
    end

    getdefault(rising) ≥ getdefault(falling) ||
        throw(ArgumentError("`rising` must be smaller than `falling`"))
    getdefault(Td) > 0 ||
        throw(ArgumentError("Time constant `Td` must be strictly positive"))

    @named siso = SISO(; y_start)
    @unpack y, u = siso

    eqs = [
        D(y) ~ max(min((u - y) / Td, rising), falling)
    ]

    initialization_eqs = [
        y ~ y_start
    ]

    return extend(System(eqs, t, [], pars; name, initialization_eqs), siso)
end
