using ModelingToolkit, InfiniteOpt, JuMP, Ipopt
using ModelingToolkit: D_nounits as D, t_nounits as t, varmap_to_vars

@mtkmodel Pendulum begin
    @parameters begin
        g = 9.8
        L = 0.4
        K = 1.2
        m = 0.3
    end
    @variables begin
        θ(t) # state
        ω(t) # state
        τ(t) = 0 # input
        y(t) # output
    end
    @equations begin
        D(θ) ~ ω
        D(ω) ~ -g / L * sin(θ) - K / m * ω + τ / m / L^2
        y ~ θ * 180 / π
    end
end
@named model = Pendulum()
model = complete(model)

inputs = [model.τ]
(f_oop, f_ip), dvs, psym, io_sys = ModelingToolkit.generate_control_function(
    model, inputs, split = false)

outputs = [model.y]
f_obs = ModelingToolkit.build_explicit_observed_function(io_sys, outputs; inputs = inputs)

expected_state_order = [model.θ, model.ω]
permutation = [findfirst(isequal(x), expected_state_order) for x in dvs] # This maps our expected state order to the actual state order

##

ub = varmap_to_vars([model.θ => 2pi, model.ω => 10], dvs)
lb = varmap_to_vars([model.θ => -2pi, model.ω => -10], dvs)
xf = varmap_to_vars([model.θ => pi, model.ω => 0], dvs)
nx = length(dvs)
nu = length(inputs)
ny = length(outputs)

##
m = InfiniteModel(optimizer_with_attributes(Ipopt.Optimizer,
    "print_level" => 0, "acceptable_tol" => 1e-3, "constr_viol_tol" => 1e-5, "max_iter" => 1000,
    "tol" => 1e-5, "mu_strategy" => "monotone", "nlp_scaling_method" => "gradient-based",
    "alpha_for_y" => "safer-min-dual-infeas", "bound_mult_init_method" => "mu-based", "print_user_options" => "yes"));

@infinite_parameter(m, τ in [0, 1], num_supports=51,
    derivative_method=OrthogonalCollocation(4)) # Time variable
guess_xs = [t -> pi, t -> 0.1][permutation]
guess_us = [t -> 0.1]
InfiniteOpt.@variables(m,
    begin
        # state variables
        (lb[i] <= x[i = 1:nx] <= ub[i], Infinite(τ), start = guess_xs[i]) # state variables
        -10 <= u[i = 1:nu] <= 10, Infinite(τ), (start = guess_us[i]) # control variables
        0 <= tf <= 10, (start = 5) # Final time
        0.2 <= L <= 0.6, (start = 0.4) # Length parameter
    end)

# Trace the dynamics
x0, p = ModelingToolkit.get_u0_p(io_sys, [model.θ => 0, model.ω => 0], [model.L => L])

xp = f_oop(x, u, p, τ)
cp = f_obs(x, u, p, τ) # Test that it's possible to trace through an observed function

@objective(m, Min, tf)
@constraint(m, [i = 1:nx], x[i](0)==x0[i]) # Initial condition
@constraint(m, [i = 1:nx], x[i](1)==xf[i]) # Terminal state

x_scale = varmap_to_vars([model.θ => 1
                          model.ω => 1], dvs)

# Add dynamics constraints
@constraint(m, [i = 1:nx], (∂(x[i], τ) - tf * xp[i]) / x_scale[i]==0)

optimize!(m)

# Extract the optimal solution
opt_tf = value(tf)
opt_time = opt_tf * value(τ)
opt_x = [value(x[i]) for i in permutation]
opt_u = [value(u[i]) for i in 1:nu]
opt_L = value(L)

# Plot the results
# using Plots
# plot(opt_time, opt_x[1], label = "θ", xlabel = "Time [s]", layout=3)
# plot!(opt_time, opt_x[2], label = "ω", sp=2)
# plot!(opt_time, opt_u[1], label = "τ", sp=3)

using Test
@test opt_x[1][end]≈pi atol=1e-3
@test opt_x[2][end]≈0 atol=1e-3

@test opt_x[1][1]≈0 atol=1e-3
@test opt_x[2][1]≈0 atol=1e-3

@test opt_L≈0.2 atol=1e-3 # Smallest permissible length is optimal
