using BifurcationKit, ModelingToolkit, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
# Simple pitchfork diagram, compares solution to native BifurcationKit, checks they are identical.
# Checks using `jac=false` option.
let
    # Creates model.
    @variables x(t) y(t)
    @parameters μ α
    eqs = [0 ~ μ * x - x^3 + α * y,
        0 ~ -y]
    @named nsys = NonlinearSystem(eqs, [x, y], [μ, α])
    nsys = complete(nsys)
    # Creates BifurcationProblem 
    bif_par = μ
    p_start = [μ => -1.0, α => 1.0]
    u0_guess = [x => 1.0, y => 1.0]
    plot_var = x
    bprob = BifurcationProblem(nsys,
        u0_guess,
        p_start,
        bif_par;
        plot_var = plot_var,
        jac = false)

    # Conputes bifurcation diagram.
    p_span = (-4.0, 6.0)
    opts_br = ContinuationPar(max_steps = 500, p_min = p_span[1], p_max = p_span[2])
    bif_dia = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside = true)

    # Computes bifurcation diagram using BifurcationKit directly (without going through MTK).
    function f_BK(u, p)
        x, y = u
        μ, α = p
        return [μ * x - x^3 + α * y, -y]
    end
    bprob_BK = BifurcationProblem(f_BK,
        [1.0, 1.0],
        [-1.0, 1.0],
        (BifurcationKit.@optic _[1]);
        record_from_solution = (x, p; k...) -> x[1])
    bif_dia_BK = bifurcationdiagram(bprob_BK,
        PALC(),
        2,
        (args...) -> opts_br;
        bothside = true)

    # Compares results.
    @test getfield.(bif_dia.γ.branch, :x) ≈ getfield.(bif_dia_BK.γ.branch, :x)
    @test getfield.(bif_dia.γ.branch, :param) ≈ getfield.(bif_dia_BK.γ.branch, :param)
    @test bif_dia.γ.specialpoint[1].x == bif_dia_BK.γ.specialpoint[1].x
    @test bif_dia.γ.specialpoint[1].param == bif_dia_BK.γ.specialpoint[1].param
    @test bif_dia.γ.specialpoint[1].type == bif_dia_BK.γ.specialpoint[1].type
end

# Lotka–Volterra model, checks exact position of bifurcation variable and bifurcation points.
# Checks using ODESystem input.
let
    # Creates a Lotka–Volterra model.
    @parameters α a b
    @variables x(t) y(t) z(t)
    eqs = [D(x) ~ -x + a * y + x^2 * y,
        D(y) ~ b - a * y - x^2 * y]
    @named sys = ODESystem(eqs, t)
    sys = complete(sys)
    # Creates BifurcationProblem
    bprob = BifurcationProblem(sys,
        [x => 1.5, y => 1.0],
        [a => 0.1, b => 0.5],
        b;
        plot_var = x)

    # Computes bifurcation diagram.
    p_span = (0.0, 2.0)
    opt_newton = NewtonPar(tol = 1e-9, max_iterations = 2000)
    opts_br = ContinuationPar(dsmax = 0.05,
        max_steps = 500,
        newton_options = opt_newton,
        p_min = p_span[1],
        p_max = p_span[2],
        n_inversion = 4)
    bif_dia = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside = true)

    # Tests that the diagram has the correct values (x = b)
    all([b.x ≈ b.param for b in bif_dia.γ.branch])

    # Tests that we get two Hopf bifurcations at the correct positions.
    hopf_points = sort(
        getfield.(filter(sp -> sp.type == :hopf, bif_dia.γ.specialpoint),
            :x);
        by = x -> x[1])
    @test length(hopf_points) == 2
    @test hopf_points[1] ≈ [0.41998733080424205, 1.5195495712453098]
    @test hopf_points[2] ≈ [0.7899715592573977, 1.0910379583813192]
end

# Simple fold bifurcation model, checks exact position of bifurcation variable and bifurcation points.
# Checks that default parameter values are accounted for.
# Checks that observables (that depend on other observables, as in this case) are accounted for.
let
    # Creates model, and uses `structural_simplify` to generate observables.
    @parameters μ p=2
    @variables x(t) y(t) z(t)
    eqs = [0 ~ μ - x^3 + 2x^2,
        0 ~ p * μ - y,
        0 ~ y - z]
    @named nsys = NonlinearSystem(eqs, [x, y, z], [μ, p])
    nsys = structural_simplify(nsys)

    # Creates BifurcationProblem.
    bif_par = μ
    p_start = [μ => 1.0]
    u0_guess = [x => 1.0, y => 0.1, z => 0.1]
    plot_var = x
    bprob = BifurcationProblem(nsys, u0_guess, p_start, bif_par; plot_var = plot_var)

    # Computes bifurcation diagram.
    p_span = (-4.3, 12.0)
    opt_newton = NewtonPar(tol = 1e-9, max_iterations = 20)
    opts_br = ContinuationPar(dsmax = 0.05,
        max_steps = 500,
        newton_options = opt_newton,
        p_min = p_span[1],
        p_max = p_span[2],
        n_inversion = 4)
    bif_dia = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside = true)

    # Tests that the diagram has the correct values (x = b)
    all([b.x ≈ 2 * b.param for b in bif_dia.γ.branch])

    # Tests that we get two fold bifurcations at the correct positions.
    fold_points = sort(getfield.(filter(sp -> sp.type == :bp, bif_dia.γ.specialpoint),
        :param))
    @test length(fold_points) == 2
    @test fold_points ≈ [-1.1851851706940317, -5.6734983580551894e-6] # test that they occur at the correct parameter values).
end

let
    @mtkmodel FOL begin
        @parameters begin
            τ # parameters
        end
        @variables begin
            x(t) # dependent variables
            RHS(t)
        end
        @equations begin
            RHS ~ τ + x^2 - 0.1
            D(x) ~ RHS
        end
    end

    @mtkbuild fol = FOL()

    par = [fol.τ => 0.0]
    u0 = [fol.x => -1.0]
    #prob = ODEProblem(fol, u0, (0.0, 1.), par)

    bif_par = fol.τ
    bp = BifurcationProblem(fol, u0, par, bif_par)
    opts_br = ContinuationPar(p_min = -1.0,
        p_max = 1.0)
    bf = bifurcationdiagram(bp, PALC(), 2, opts_br)

    @test bf.γ.specialpoint[1].param≈0.1 atol=1e-4 rtol=1e-4

    # Test with plot variable as observable
    pvar = ModelingToolkit.get_var_to_name(fol)[:RHS]
    bp = BifurcationProblem(fol, u0, par, bif_par; plot_var = pvar)
    opts_br = ContinuationPar(p_min = -1.0,
        p_max = 1.0)
    bf = bifurcationdiagram(bp, PALC(), 2, opts_br)
    @test bf.γ.specialpoint[1].param≈0.1 atol=1e-4 rtol=1e-4
end
