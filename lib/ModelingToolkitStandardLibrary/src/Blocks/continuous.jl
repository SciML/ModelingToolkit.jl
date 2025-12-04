"""
    Integrator(;name, k = 1, x = 0.0)

Outputs `y = ∫k*u dt`, corresponding to the transfer function ``1/s``.
Initial value of integrator state ``x`` can be set with `x`

# Connectors:

  - `input`
  - `output`

# Parameters:

  - `k`: Gain of integrator

# Unknowns:

  - `x`: State of Integrator. Defaults to 0.0.
"""
@mtkmodel Integrator begin
    @extend u, y = siso = SISO()
    @variables begin
        x(t) = 0.0, [description = "State of Integrator"]
    end
    @parameters begin
        k = 1, [description = "Gain"]
    end
    @equations begin
        D(x) ~ k * u
        y ~ x
    end
end

"""
    Derivative(; name, k = 1, T, x = 0.0)

Outputs an approximate derivative of the input. The transfer function of this block is

```
k      k        ks  
─ - ─────── = ────── 
T   sT² + T   sT + 1
```

and a state-space realization is given by `ss(-1/T, 1/T, -k/T, k/T)`
where `T` is the time constant of the filter.
A smaller `T` leads to a more ideal approximation of the derivative.

Initial value of the state ``x`` can be set with `x`.

# Parameters:

  - `k`: Gain
  - `T`: [s] Time constant (T>0 required; T=0 is ideal derivative block)

# Unknowns:

  - `x`: Unknown of Derivative. Defaults to 0.0.

# Connectors:

  - `input`
  - `output`
"""
@mtkmodel Derivative begin
    @extend u, y = siso = SISO()
    @variables begin
        x(t) = 0.0, [description = "Derivative-filter state"]
    end
    @parameters begin
        T = T, [description = "Time constant"]
        k = 1, [description = "Gain"]
    end
    begin
        @symcheck T > 0 ||
                  throw(ArgumentError("Time constant `T` has to be strictly positive"))
    end
    @equations begin
        D(x) ~ (u - x) / T
        y ~ (k / T) * (u - x)
    end
end

"""
    FirstOrder(; name, k = 1.0, T, x = 0.0, lowpass = true)

A first-order filter with a single real pole at `s = -1/T` and gain `k`. If `lowpass=true` (default), the transfer function
is given by ``Y(s)/U(s) = ``


```
   k
───────
sT + 1
```

and if `lowpass=false`, by

```
sT + 1 - k
──────────
  sT + 1
```

Initial value of the state `x` can be set with `x`

# Parameters:

  - `k`: Gain
  - `T`: [s] Time constant (T>0 required)

# Connectors:

  - `input`
  - `output`

See also [`SecondOrder`](@ref)
"""
@mtkmodel FirstOrder begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        lowpass = true
    end
    @variables begin
        x(t) = 0.0, [description = "State of FirstOrder filter"]
    end
    @parameters begin
        T = T, [description = "Time constant"]
        k = 1.0, [description = "Gain"]
    end
    begin
        @symcheck T > 0 ||
                  throw(ArgumentError("Time constant `T` has to be strictly positive"))
    end
    @equations begin
        D(x) ~ (k * u - x) / T
        lowpass ? y ~ x : y ~ k * u - x
    end
end

"""
    SecondOrder(; name, k = 1.0, w = 1.0, d = 1.0, x = 0.0, xd = 0.0)

A second-order filter with gain `k`, a bandwidth of `w` rad/s and relative damping `d`. The transfer function
is given by `Y(s)/U(s) = `

```
      k*w^2
─────────────────
s² + 2d*w*s + w^2
```

Critical damping corresponds to `d=1`, which yields the fastest step response without overshoot, `d < 1` results in an underdamped filter while `d > 1` results in an overdamped filter.
`d = 1/√2` corresponds to a Butterworth filter of order 2 (maximally flat frequency response).
Initial value of the state `x` can be set with `x`, and of derivative state `xd` with `xd`.

# Parameters:

  - `k`: Gain
  - `w`: [`rad/s`] Angular frequency
  - `d`: Damping

# Connectors:

  - `input`
  - `output`
"""
@mtkmodel SecondOrder begin
    @extend u, y = siso = SISO()
    @variables begin
        x(t), [description = "State of SecondOrder filter", guess = 0.0]
        xd(t), [description = "Derivative state of SecondOrder filter", guess = 0.0]
    end
    @parameters begin
        k = 1.0, [description = "Gain"]
        w = 1.0, [description = "Bandwidth (angular frequency)"]
        d = 1.0, [description = "Relative damping"]
    end
    @equations begin
        D(x) ~ xd
        D(xd) ~ w * (w * (k * u - x) - 2 * d * xd)
        y ~ x
    end
end

"""
    PI(;name, k = 1.0, T = 1.0, int.x = 0.0)

Textbook version of a PI-controller without actuator saturation and anti-windup measure.
The proportional gain can be set with `k`
Initial value of integrator state `x` can be set with `int.x`

The PI controller is implemented on standard form:
```math
U(s) = k (1 + \\dfrac{1}{sT}) E(S)
```

# Parameters:
  - `k`: Proportional gain
  - `T`: [s] Integrator time constant (T>0 required)

# Connectors:

  - `err_input`
  - `ctr_output`

See also [`LimPI`](@ref)
"""
@mtkmodel PI begin
    @parameters begin
        k = 1.0, [description = "Proportional gain"]
        T = 1.0, [description = "Integrator time constant"]
    end
    begin
        @symcheck T > 0 ||
                  throw(ArgumentError("Time constant `T` has to be strictly positive"))
    end
    @components begin
        err_input = RealInput() # control error
        ctr_output = RealOutput() # control signal
        gainPI = Gain(; k)
        addPI = Add()
        int = Integrator(k = 1 / T, x = 0.0)
    end
    @equations begin
        connect(err_input, addPI.input1)
        connect(addPI.output, gainPI.input)
        connect(gainPI.output, ctr_output)
        connect(err_input, int.input)
        connect(int.output, addPI.input2)
    end
end

"""
    PID(;name, k=1, Ti=false, Td=false, Nd=10, int__x=0, der__x=0)

Text-book version of a PID-controller without actuator saturation and anti-windup measure.

# Parameters:

  - `k`: Gain
  - `Ti`: [s] Integrator time constant (Ti>0 required). If set to false, no integral action is used.
  - `Td`: [s] Derivative time constant (Td>0 required). If set to false, no derivative action is used.
  - `Nd`: [s] Time constant for the derivative approximation (Nd>0 required; Nd=0 is ideal derivative).
  - `int__x`: Initial value for the integrator.
  - `der__x`: Initial value for the derivative state.

# Connectors:

  - `err_input`
  - `ctr_output`

See also [`LimPID`](@ref)
"""
@component function PID(; name, k = 1, Ti = false, Td = false, Nd = 10, int__x = 0,
        der__x = 0)
    with_I = !isequal(Ti, false)
    with_D = !isequal(Td, false)
    @named err_input = RealInput() # control error
    @named ctr_output = RealOutput() # control signal
    @symcheck Ti ≥ 0 ||
              throw(ArgumentError("Ti out of bounds, got $(Ti) but expected Ti ≥ 0"))
    @symcheck Td ≥ 0 ||
              throw(ArgumentError("Td out of bounds, got $(Td) but expected Td ≥ 0"))
    @symcheck Nd > 0 ||
              throw(ArgumentError("Nd out of bounds, got $(Nd) but expected Nd > 0"))

    pars = @parameters begin
        k = k, [description = "Proportional gain"]
        Ti = Ti, [description = "Integrator time constant"]
        Td = Td, [description = "Derivative time constant"]
        Nd = Nd, [description = "Derivative limit"]
    end

    @named gainPID = Gain(; k)
    @named addPID = Add3()
    if with_I
        @named int = Integrator(k = 1 / Ti, x = int__x)
    else
        @named Izero = Constant(k = 0)
    end
    if with_D
        @named der = Derivative(k = Td, T = 1 / Nd, x = der__x)
    else
        @named Dzero = Constant(k = 0)
    end
    sys = [err_input, ctr_output, gainPID, addPID]
    if with_I
        push!(sys, int)
    else
        push!(sys, Izero)
    end
    if with_D
        push!(sys, der)
    else
        push!(sys, Dzero)
    end
    eqs = [
        connect(err_input, addPID.input1),
        connect(addPID.output, gainPID.input),
        connect(gainPID.output, ctr_output)
    ]
    if with_I
        push!(eqs, connect(err_input, int.input))
        push!(eqs, connect(int.output, addPID.input2))
    else
        push!(eqs, connect(Izero.output, addPID.input2))
    end
    if with_D
        push!(eqs, connect(err_input, der.input))
        push!(eqs, connect(der.output, addPID.input3))
    else
        push!(eqs, connect(Dzero.output, addPID.input3))
    end
    System(eqs, t, [], pars; name = name, systems = sys)
end

"""
    LimPI(; name, k = 1.0, T, Ta, int__x = 0.0, u_max = 1.0, u_min = -u_max)

Text-book version of a PI-controller with actuator saturation and anti-windup measure.

The PI controller is implemented on standard form
```math
u(t) = sat(k (e(t) + ∫\\dfrac{1}{T}e(t) dt) )
```
The simplified expression above is given without the anti-windup protection.

# Parameters:

  - `k`: Proportional gain
  - `T`: [s] Integrator time constant (T>0 required)
  - `Ta`: [s] Tracking time constant (Ta>0 required)

# Connectors:

  - `err_input`
  - `ctr_output`
"""
@component function LimPI(; name, k = 1, T, u_max, u_min = -u_max, Ta, int__x = 0.0)
    @symcheck Ta > 0 ||
              throw(ArgumentError("Time constant `Ta` has to be strictly positive"))
    @symcheck T > 0 || throw(ArgumentError("Time constant `T` has to be strictly positive"))
    @symcheck u_max ≥ u_min || throw(ArgumentError("u_min must be smaller than u_max"))
    pars = @parameters begin
        k = k, [description = "Proportional gain"]
        T = T, [description = "Integrator time constant"]
        Ta = Ta, [description = "Tracking time constant"]
        u_max = u_max, [description = "Upper saturation limit"]
        u_min = u_min, [description = "Lower saturation limit"]
    end
    @named err_input = RealInput() # control error
    @named ctr_output = RealOutput() # control signal
    @named gainPI = Gain(; k)
    @named addPI = Add()
    @named addTrack = Add()
    @named int = Integrator(k = 1 / T, x = int__x)
    @named limiter = Limiter(y_max = u_max, y_min = u_min)
    @named addSat = Add(k1 = 1, k2 = -1)
    @named gainTrack = Gain(k = 1 / Ta)
    sys = [err_input, ctr_output, gainPI, addPI, int, addTrack, limiter, addSat, gainTrack]
    eqs = [
        connect(err_input, addPI.input1),
        connect(addPI.output, gainPI.input),
        connect(gainPI.output, limiter.input),
        connect(limiter.output, ctr_output),
        connect(limiter.input, addSat.input2),
        connect(limiter.output, addSat.input1),
        connect(addSat.output, gainTrack.input),
        connect(err_input, addTrack.input1),
        connect(gainTrack.output, addTrack.input2),
        connect(addTrack.output, int.input),
        connect(int.output, addPI.input2)
    ]
    System(eqs, t, [], pars; name = name, systems = sys)
end

"""
    LimPID(; k, Ti=false, Td=false, wp=1, wd=1, Ni, Nd=12, u_max=Inf, u_min=-u_max, gains = false, name)

Proportional-Integral-Derivative (PID) controller with output saturation, set-point weighting and integrator anti-windup.

The equation for the control signal is roughly

```
k(ep + 1/Ti * ∫e + Td * d/dt(ed))
e = u_r - u_y
ep = wp*u_r - u_y
ed = wd*u_r - u_y
```

where the transfer function for the derivative includes additional filtering, see `? Derivative` for more details.

# Parameters:

  - `k`: Proportional gain
  - `Ti`: [s] Integrator time constant. Set to `false` to turn off integral action.
  - `Td`: [s] Derivative time constant. Set to `false` to turn off derivative action.
  - `wp`: [0,1] Set-point weighting in the proportional part.
  - `wd`: [0,1] Set-point weighting in the derivative part.
  - `Nd`: [1/s] Derivative limit, limits the derivative gain to Nd/Td. Reasonable values are ∈ [8, 20]. A higher value gives a better approximation of an ideal derivative at the expense of higher noise amplification.
  - `Ni`: `Ni*Ti` controls the time constant `Ta` of anti-windup tracking. A common (default) choice is `Ta = √(Ti*Td)` which is realized by `Ni = √(Td / Ti)`. Anti-windup can be effectively turned off by setting `Ni = Inf`.
  - `gains`: If `gains = true`, `Ti` and `Td` will be interpreted as gains with a fundamental PID transfer function on parallel form `ki=Ti, kd=Td, k + ki/s + kd*s`.

# Connectors:

  - `reference`
  - `measurement`
  - `ctr_output`
"""
@component function LimPID(; name, k = 1, Ti = false, Td = false, wp = 1, wd = 1,
        Ni = Ti == 0 ? Inf : √(max(Td / Ti, 1e-6)),
        Nd = 10,
        u_max = Inf,
        u_min = u_max > 0 ? -u_max : -Inf,
        gains = false,
        int__x = 0.0,
        der__x = 0.0)
    with_I = !isequal(Ti, false)
    with_D = !isequal(Td, false)
    with_AWM = Ni != Inf
    if gains
        Ti = k / Ti
        Td = Td / k
    end
    @symcheck Ti ≥ 0 ||
              throw(ArgumentError("Ti out of bounds, got $(Ti) but expected Ti ≥ 0"))
    @symcheck Td ≥ 0 ||
              throw(ArgumentError("Td out of bounds, got $(Td) but expected Td ≥ 0"))
    @symcheck u_max ≥ u_min || throw(ArgumentError("u_min must be smaller than u_max"))
    @symcheck Nd > 0 ||
              throw(ArgumentError("Nd out of bounds, got $(Nd) but expected Nd > 0"))

    pars = @parameters begin
        k = k, [description = "Proportional gain"]
        Ti = Ti, [description = "Integrator time constant"]
        Td = Td, [description = "Derivative time constant"]
        wp = wp, [description = "Set-point weighting in the proportional part"]
        wd = wd, [description = "Set-point weighting in the derivative part"]
        Ni = Ni, [description = "Anti-windup tracking gain"]
        Nd = Nd, [description = "Derivative limit"]
        u_max = u_max, [description = "Upper saturation limit"]
        u_min = u_min, [description = "Lower saturation limit"]
    end
    @named reference = RealInput()
    @named measurement = RealInput()
    @named ctr_output = RealOutput() # control signal
    @named addP = Add(k1 = wp, k2 = -1)
    @named gainPID = Gain(; k)
    @named addPID = Add3()
    @named limiter = Limiter(y_max = u_max, y_min = u_min)
    if with_I
        if with_AWM
            @named addI = Add3(k1 = 1, k2 = -1, k3 = 1)
            @named addSat = Add(k1 = 1, k2 = -1)
            @named gainTrack = Gain(k = 1 / (k * Ni))
        else
            @named addI = Add(k1 = 1, k2 = -1)
        end
        @named int = Integrator(k = 1 / Ti, x = int__x)
    else
        @named Izero = Constant(k = 0)
    end
    if with_D
        @named der = Derivative(k = Td, T = 1 / Nd, x = der__x)
        @named addD = Add(k1 = wd, k2 = -1)
    else
        @named Dzero = Constant(k = 0)
    end

    sys = [reference, measurement, ctr_output, addP, gainPID, addPID, limiter]
    if with_I
        if with_AWM
            push!(sys, [addSat, gainTrack]...)
        end
        push!(sys, [addI, int]...)
    else
        push!(sys, Izero)
    end
    if with_D
        push!(sys, [addD, der]...)
    else
        push!(sys, Dzero)
    end

    eqs = [
        connect(reference, addP.input1),
        connect(measurement, addP.input2),
        connect(addP.output, addPID.input1),
        connect(addPID.output, gainPID.input),
        connect(gainPID.output, limiter.input),
        connect(limiter.output, ctr_output)
    ]
    if with_I
        push!(eqs, connect(reference, addI.input1))
        push!(eqs, connect(measurement, addI.input2))
        if with_AWM
            push!(eqs, connect(limiter.input, addSat.input2))
            push!(eqs, connect(limiter.output, addSat.input1))
            push!(eqs, connect(addSat.output, gainTrack.input))
            push!(eqs, connect(gainTrack.output, addI.input3))
        end
        push!(eqs, connect(addI.output, int.input))
        push!(eqs, connect(int.output, addPID.input3))
    else
        push!(eqs, connect(Izero.output, addPID.input3))
    end
    if with_D
        push!(eqs, connect(reference, addD.input1))
        push!(eqs, connect(measurement, addD.input2))
        push!(eqs, connect(addD.output, der.input))
        push!(eqs, connect(der.output, addPID.input2))
    else
        push!(eqs, connect(Dzero.output, addPID.input2))
    end

    System(eqs, t, [], pars; name = name, systems = sys)
end

"""
    StateSpace(A, B, C, D = 0; x = zeros(size(A,1)), u0 = zeros(size(B,2)), y0 = zeros(size(C,1)), name)

A linear, time-invariant state-space system on the form.

```math
\\begin{aligned}
ẋ &= Ax + Bu \\\\
y &= Cx + Du
\\end{aligned}
```

Transfer functions can also be simulated by converting them to a StateSpace form.

`y0` and `u0` can be used to set an operating point, providing them changes the dynamics from an LTI system to the affine system

```math
\\begin{aligned}
ẋ &= Ax + B(u - u0) \\\\
y &= Cx + D(u - u0) + y0
\\end{aligned}
```

For a nonlinear system

```math
\\begin{aligned}
ẋ &= f(x, u) \\\\
y &= h(x, u)
\\end{aligned}
```

linearized around the operating point `x₀, u₀`, we have `y0, u0 = h(x₀, u₀), u₀`.
"""
@component function StateSpace(; A, B, C, D = nothing, x = zeros(size(A, 1)), name,
        u0 = zeros(size(B, 2)), y0 = zeros(size(C, 1)))
    nx, nu, ny = size(A, 1), size(B, 2), size(C, 1)
    size(A, 2) == nx || error("`A` has to be a square matrix.")
    size(B, 1) == nx || error("`B` has to be of dimension ($nx x $nu).")
    size(C, 2) == nx || error("`C` has to be of dimension ($ny x $nx).")
    if B isa AbstractVector
        B = reshape(B, length(B), 1)
    end
    if isnothing(D) || iszero(D)
        D = zeros(ny, nu)
    else
        size(D) == (ny, nu) || error("`D` has to be of dimension ($ny x $nu).")
    end
    @named input = RealInput(nin = nu)
    @named output = RealOutput(nout = ny)
    @variables x(t)[1:nx]=x [
        description = "State variables of StateSpace system $name"
    ]
    # pars = @parameters A=A B=B C=C D=D # This is buggy
    eqs = [ # FIXME: if array equations work
        [Differential(t)(x[i]) ~
         sum(A[i, k] * x[k] for k in 1:nx) +
         sum(B[i, j] * (input.u[j] - u0[j]) for j in 1:nu)
         for i in 1:nx]..., # cannot use D here
        [output.u[j] ~
         sum(C[j, i] * x[i] for i in 1:nx) +
         sum(D[j, k] * (input.u[k] - u0[k]) for k in 1:nu) + y0[j]
         for j in 1:ny]...
    ]
    compose(System(eqs, t, vcat(x...), [], name = name), [input, output])
end

StateSpace(A, B, C, D = nothing; kwargs...) = StateSpace(; A, B, C, D, kwargs...)

symbolic_eps(t) = eps(t)
@register_symbolic symbolic_eps(t)

"""
    TransferFunction(; b, a, name)

A single input, single output, linear time-invariant system provided as a transfer-function.
```
Y(s) = b(s) / a(s)  U(s)
```
where `b` and `a` are vectors of coefficients of the numerator and denominator polynomials, respectively, ordered such that the coefficient of the highest power of `s` is first.

The internal state realization is on controller canonical form, with state variable `x`, output variable `y` and input variable `u`. For numerical robustness, the realization used by the integrator is scaled by the last entry of the `a` parameter. The internally scaled state variable is available as `x_scaled`.

To set the initial state, it's recommended to set the initial condition for `x`, and let that of `x_scaled` be computed automatically.

# Parameters:
- `b`: Numerator polynomial coefficients, e.g., `2s + 3` is specified as `[2, 3]`
- `a`: Denominator polynomial coefficients, e.g., `s² + 2ωs + ω^2` is specified as `[1, 2ω, ω^2]`

# Connectors:
  - `input`
  - `output`

See also [`StateSpace`](@ref) which handles MIMO systems, as well as [ControlSystemsMTK.jl](https://juliacontrol.github.io/ControlSystemsMTK.jl/stable/) for an interface between [ControlSystems.jl](https://juliacontrol.github.io/ControlSystems.jl/stable/) and ModelingToolkit.jl for advanced manipulation of transfer functions and linear statespace systems. For linearization, see [`linearize`](@ref) and [Linear Analysis](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/linear_analysis/).
"""
@component function TransferFunction(; b = [1], a = [1, 1], name)
    nb = length(b)
    na = length(a)
    nb <= na ||
        error("Transfer function is not proper, the numerator must not be longer than the denominator")
    nx = na - 1
    nbb = max(0, na - nb)

    @named begin
        input = RealInput()
        output = RealOutput()
    end

    @parameters begin
        b[1:nb] = b,
        [
            description = "Numerator coefficients of transfer function (e.g., 2s + 3 is specified as [2,3])"
        ]
        a[1:na] = a,
        [
            description = "Denominator coefficients of transfer function (e.g., `s² + 2ωs + ω^2` is specified as [1, 2ω, ω^2])"
        ]
        bb[1:(nbb + nb)] = [zeros(nbb); b]
    end
    d = bb[1] / a[1]# , [description = "Direct feedthrough gain"]

    a = collect(a)
    a_end = ifelse(a[end] > 100 * symbolic_eps(sqrt(a' * a)), a[end], 1.0)

    pars = [collect(b); a; collect(bb)]
    @variables begin
        x(t)[1:nx] = zeros(nx),
        [description = "State of transfer function on controller canonical form"]
        x_scaled(t)[1:nx] = collect(x) * a_end, [description = "Scaled vector x"]
        u(t), [description = "Input of transfer function"]
        y(t), [description = "Output of transfer function"]
    end

    x = collect(x)
    x_scaled = collect(x_scaled)
    bb = collect(bb)

    sts = [x; x_scaled; y; u]

    if nx == 0
        eqs = [y ~ d * u]
    else
        eqs = Equation[D(x_scaled[1]) ~ (-a[2:na]'x_scaled + a_end * u) / a[1]
                       D.(x_scaled[2:nx]) .~ x_scaled[1:(nx - 1)]
                       y ~ ((bb[2:na] - d * a[2:na])'x_scaled) / a_end + d * u
                       x .~ x_scaled ./ a_end]
    end
    push!(eqs, input.u ~ u)
    push!(eqs, output.u ~ y)
    compose(System(eqs, t, sts, pars; name = name), input, output)
end
