using Test
using Latexify
using ModelingToolkit

### Tips for generating latex tests:
### Latexify has an unexported macro:
###
### Latexify.@generate_test
###
### which generates a test using a given latexify function.
### For example:
###
### Latexify.@generate_test latexify([1, 2, 3], [4, 5, 6]; env=:mdtable)
###
### This puts a ready-made test in your clipboard which you can paste into the
### test file.
###
### Just be sure to remove all such macros before you commit a change since it
### will cause issues with Travis.

@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x)*D(x-y)/D(z),
       0 ~ σ*x*(ρ-z)/10-y,
       D(z) ~ x*y^(2//3) - β*z]



# Latexify.@generate_test latexify(eqs)
@test latexify(eqs) == replace(
raw"\begin{align}
\frac{dx(t)}{dt} =& \frac{\sigma \left( y(t) - x(t) \right) \frac{d\left( x(t) - y(t) \right)}{dt}}{\frac{dz(t)}{dt}} \\
0 =& \frac{\sigma x(t) \left( \rho - z(t) \right)}{10} - y(t) \\
\frac{dz(t)}{dt} =& x(t) y(t)^{\frac{2}{3}} - \beta z(t)
\end{align}
", "\r\n"=>"\n")


@variables u[1:3](t)
@parameters p[1:3]
eqs = [D(u[1]) ~ p[3]*(u[2]-u[1]),
       0 ~ p[2]*p[3]*u[1]*(p[1]-u[1])/10-u[2],
       D(u[3]) ~ u[1]*u[2]^(2//3) - p[3]*u[3]]

# Latexify.@generate_test latexify(eqs)
@test latexify(eqs) == replace(
raw"\begin{align}
\frac{du{_1}(t)}{dt} =& p{_3} \left( u{_2}(t) - u{_1}(t) \right) \\
0 =& \frac{p{_2} p{_3} u{_1}(t) \left( p{_1} - u{_1}(t) \right)}{10} - u{_2}(t) \\
\frac{du{_3}(t)}{dt} =& u{_1}(t) u{_2}(t)^{\frac{2}{3}} - p{_3} u{_3}(t)
\end{align}
", "\r\n"=>"\n")


eqs = [D(u[1]) ~ p[3]*(u[2]-u[1]),
       D(u[2]) ~ p[2]*p[3]*u[1]*(p[1]-u[1])/10-u[2],
       D(u[3]) ~ u[1]*u[2]^(2//3) - p[3]*u[3]]

sys = ODESystem(eqs)

# Latexify.@generate_test latexify(eqs; show_iv=false)
@test latexify(eqs; show_iv = false) == replace(
raw"\begin{align}
\frac{du{_1}}{dt} =& p{_3} \left( u{_2} - u{_1} \right) \\
\frac{du{_2}}{dt} =& \frac{p{_2} p{_3} u{_1} \left( p{_1} - u{_1} \right)}{10} - u{_2} \\
\frac{du{_3}}{dt} =& u{_1} u{_2}^{\frac{2}{3}} - p{_3} u{_3}
\end{align}
", "\r\n"=>"\n")

@test latexify(eqs; show_iv=false) == latexify(sys; show_iv=false)


@test latexify(eqs) == replace(
raw"\begin{align}
\frac{du{_1}(t)}{dt} =& p{_3} \left( u{_2}(t) - u{_1}(t) \right) \\
\frac{du{_2}(t)}{dt} =& \frac{p{_2} p{_3} u{_1}(t) \left( p{_1} - u{_1}(t) \right)}{10} - u{_2}(t) \\
\frac{du{_3}(t)}{dt} =& u{_1}(t) u{_2}(t)^{\frac{2}{3}} - p{_3} u{_3}(t)
\end{align}
", "\r\n"=>"\n")



# Latexify.@generate_test latexify(sys; show_iv=false, cdot = true)
@test latexify(sys; show_iv = false, cdot = true) == replace(
raw"\begin{align}
\frac{du{_1}}{dt} =& p{_3} \cdot \left( u{_2} - u{_1} \right) \\
\frac{du{_2}}{dt} =& \frac{p{_2} \cdot p{_3} \cdot u{_1} \cdot \left( p{_1} - u{_1} \right)}{10} - u{_2} \\
\frac{du{_3}}{dt} =& u{_1} \cdot u{_2}^{\frac{2}{3}} - p{_3} \cdot u{_3}
\end{align}
", "\r\n"=>"\n")

# Latexify.@generate_test latexify(sys; show_iv=false, cdot = false)
@test latexify(sys; show_iv = false, cdot = false) == replace(
raw"\begin{align}
\frac{du{_1}}{dt} =& p{_3} \left( u{_2} - u{_1} \right) \\
\frac{du{_2}}{dt} =& \frac{p{_2} p{_3} u{_1} \left( p{_1} - u{_1} \right)}{10} - u{_2} \\
\frac{du{_3}}{dt} =& u{_1} u{_2}^{\frac{2}{3}} - p{_3} u{_3}
\end{align}
", "\r\n"=>"\n")


@test latexify(sys; cdot=true) == latexify(equations(sys); cdot=true)

@parameters t
@variables x(t)
@derivatives D'~t
eqs = [D(x) ~ (1+cos(t))/(1+2*x)]
eqs[1].rhs.args[1].args[2].f isa Function

# Latexify.@generate_test latexify(eqs)
@test latexify(eqs) == replace(
raw"\begin{align}
\frac{dx(t)}{dt} =& \frac{1 + \cos\left( t \right)}{1 + 2 x(t)}
\end{align}
", "\r\n"=>"\n")


@parameters t
@variables x(t)
@derivatives D'~t
eqs = [D(x) ~ D(D(D(x)))]

# Latexify.@generate_test latexify(eqs)
@test latexify(eqs) == replace(
raw"\begin{align}
\frac{dx(t)}{dt} =& \frac{d^{3}x(t)}{dt^{3}}
\end{align}
", "\r\n"=>"\n")

# Latexify.@generate_test latexify(eqs; show_iv=false)
@test latexify(eqs; show_iv = false) == replace(
raw"\begin{align}
\frac{dx}{dt} =& \frac{d^{3}x}{dt^{3}}
\end{align}
", "\r\n"=>"\n")


@parameters t, x
@variables u(..)
@derivatives Dt'~t
@derivatives Dx'~x

eqs = [Dt(u(t, x)) ~ Dx(Dx(u(t, x)))] 
# Latexify.@generate_test latexify(eqs)
@test latexify(eqs) == replace(
raw"\begin{align}
\frac{du(t,x)}{dt} =& \frac{d^{2}u(t,x)}{dx^{2}}
\end{align}
", "\r\n"=>"\n")

eqs = [Dt(u(t, x)) ~ Dx(Dx(Dt(Dx(Dt(u(t, x))))))] 
@test latexify(eqs) == replace(
raw"\begin{align}
\frac{du(t,x)}{dt} =& \frac{d^{5}u(t,x)}{dx^{3}dt^{2}}
\end{align}
", "\r\n"=>"\n")
