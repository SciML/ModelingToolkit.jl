"""
```julia
Not(; name)
```

NOT gate in 9-level logic.

# Connectors

  - `x`
    Input [`DigitalPin`](@ref)
  - `y`
    Output [`DigitalPin`](@ref)
"""
function Not(; name)
    @named x = DigitalPin()
    @named y = DigitalPin()

    eqs = [x.i ~ y.i
           y.val ~ _not(x.val)]
    System(eqs, t, [], [], systems = [x, y], name = name)
end

"""
```julia
And(; name, N = 2)
```

AND gate in 9-level logic, with `N` inputs

# Connectors

  - `x1`, `x2`, ...
    `N` input [`DigitalPin`](@ref)s
  - `y`
    Output [`DigitalPin`](@ref)
"""
function And(; name, N = 2)
    x = map(1:N) do i
        DigitalPin(name = Symbol(:x, i))
    end
    @named y = DigitalPin()

    vals = [k.val for k in x]
    eqs = [y.val ~ _and(vals...)
           y.i ~ sum(k -> k.i, x)]
    System(eqs, t, [], [], systems = [x..., y], name = name)
end

"""
```julia
Nand(; name, N = 2)
```

NAND gate in 9-level logic, with `N` inputs

# Connectors

  - `x1`, `x2`, ...
    `N` input [`DigitalPin`](@ref)s
  - `y`
    Output [`DigitalPin`](@ref)
"""
function Nand(; name, N = 2)
    x = map(1:N) do i
        DigitalPin(name = Symbol(:x, i))
    end
    @named y = DigitalPin()

    vlist = [k.val for k in x]
    eqs = [y.val ~ _not(_and(vlist...))
           y.i ~ sum(k -> k.i, x)]
    System(eqs, t, [], [], systems = [x..., y], name = name)
end

"""
```julia
Or(; name, N = 2)
```

OR gate in 9-level logic, with `N` inputs

# Connectors

  - `x1`, `x2`, ...
    `N` input [`DigitalPin`](@ref)s
  - `y`
    Output [`DigitalPin`](@ref)
"""
function Or(; name, N = 2)
    x = map(1:N) do i
        DigitalPin(name = Symbol(:x, i))
    end
    @named y = DigitalPin()

    vals = [k.val for k in x]
    eqs = [y.val ~ _or(vals...)
           y.i ~ sum(k -> k.i, x)]
    System(eqs, t, [], [], systems = [x..., y], name = name)
end

"""
```julia
Nor(; name, N = 2)
```

NOR gate in 9-level logic, with `N` inputs

# Connectors

  - `x1`, `x2`, ...
    `N` input [`DigitalPin`](@ref)s
  - `y`
    Output [`DigitalPin`](@ref)
"""
function Nor(; name, N = 2)
    x = map(1:N) do i
        DigitalPin(name = Symbol(:x, i))
    end
    @named y = DigitalPin()

    vlist = [k.val for k in x]
    eqs = [y.val ~ _not(_or(vlist...))
           y.i ~ sum(k -> k.i, x)]
    System(eqs, t, [], [], systems = [x..., y], name = name)
end

"""
```julia
Xor(; name, N = 2)
```

XOR gate in 9-level logic, with `N` inputs

# Connectors

  - `x1`, `x2`, ...
    `N` input [`DigitalPin`](@ref)s
  - `y`
    Output [`DigitalPin`](@ref)
"""
function Xor(; name, N = 2)
    x = map(1:N) do i
        DigitalPin(name = Symbol(:x, i))
    end
    @named y = DigitalPin()

    vals = [k.val for k in x]
    eqs = [y.val ~ _xor(vals...)
           y.i ~ sum(k -> k.i, x)]
    System(eqs, t, [], [], systems = [x..., y], name = name)
end

"""
```julia
Xnor(; name, N = 2)
```

XNOR gate in 9-level logic, with `N` inputs

# Connectors

  - `x1`, `x2`, ...
    `N` input [`DigitalPin`](@ref)s
  - `y`
    Output [`DigitalPin`](@ref)
"""
function Xnor(; name, N = 2)
    x = map(1:N) do i
        DigitalPin(name = Symbol(:x, i))
    end
    @named y = DigitalPin()

    vlist = [k.val for k in x]
    eqs = [y.val ~ _not(_xor(vlist...))
           y.i ~ sum(k -> k.i, x)]
    System(eqs, t, [], [], systems = [x..., y], name = name)
end
