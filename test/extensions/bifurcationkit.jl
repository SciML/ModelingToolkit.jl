using BifurcationKit, ModelingToolkit, Test

# Checks pitchfork diagram and that there are the correct number of branches (a main one and two children)
let 
    @variables t x(t) y(t)
    @parameters μ α
    eqs = [0 ~ μ*x - x^3 + α*y,
            0 ~ -y]
    @named nsys = NonlinearSystem(eqs, [x, y], [μ, α])

    bif_par = μ
    p_start = [μ => -1.0, α => 1.0]
    u0_guess = [x => 1.0, y => 1.0]
    plot_var = x;

    using BifurcationKit
    bprob = BifurcationProblem(nsys, u0_guess, p_start, bif_par; plot_var=plot_var, jac=false)

    p_span = (-4.0, 6.0)
    opt_newton = NewtonPar(tol = 1e-9, max_iterations = 20)
    opts_br = ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds = 0.01,
        max_steps = 100, nev = 2, newton_options = opt_newton,
        p_min = p_span[1], p_max = p_span[2],
        detect_bifurcation = 3, n_inversion = 4, tol_bisection_eigenvalue = 1e-8, dsmin_bisection = 1e-9);

    bf = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside=true)

    @test length(bf.child) == 2
end