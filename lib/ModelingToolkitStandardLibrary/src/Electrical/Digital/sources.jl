"""
```julia
PulseDiff(; name, Val = 1, dt = 0.1)
```

# States

  - `val(t)`
    Output value of the source

# Connectors

  - `d`
    Output [`DigitalPin`](@ref)
"""
function PulseDiff(; name, Val = 1, dt = 0.1)
    @named d = DigitalPin()
    @variables val(t)
    D = ModelingToolkit.Difference(t; dt = dt)

    eqs = [D(val) ~ Val
           val ~ d.val]

    System(eqs, t, [val], [], systems = [d], defaults = Dict(Val => 0), name = name)
end

"""
```julia
Set(; name)
```

Source that outputs a constant signal of `1`.

# Connectors

  - `d`
    Output [`DigitalPin`](@ref)
"""
function Set(; name)
    @named d = DigitalPin()

    eqs = [
        d.val ~ 1
    ]
    System(eqs, t, [], [], systems = [d], name = name)
end

"""
```julia
Reset(; name)
```

Source that outputs a constant signal of `1`

# Connectors

  - `d`
    Output [`DigitalPin`](@ref)
"""
function Reset(; name)
    @named d = DigitalPin()

    eqs = [
        d.val ~ 0
    ]
    System(eqs, t, [], [], systems = [d], name = name)
end

"""
```julia
Pulse(; name, duty_cycle = 0.5, T = 1.0)
```

Pulse output with specified `duty_cycle` and time period (`T`)

# Connectors

  - `d`
    Output [`DigitalPin`](@ref)
"""
function Pulse(; name, duty_cycle = 0.5, T = 1.0)
    @named d = DigitalPin()

    eqs = [
        d.val ~ IfElse.ifelse(t % T > duty_cycle * T, 1, 0)
    ]
    System(eqs, t, [], [], systems = [d], name = name)
end
