using ModelingToolkit, Test
using Symbolics
using StaticArrays

rhs(eq) = simplify(eq.rhs)

@testset "input_affine_form" begin
    # Test with simple linear system
    @testset "Simple linear system" begin
        @variables x1 x2 u1 u2
        state = [x1, x2]
        inputs = [u1, u2]

        eqs = D.(state) .~ [
            -x1 + 2*x2 + u1,
            x1*x2 - x2 + u1 + 2*u2
        ]

        f, g = input_affine_form(eqs, inputs)

        # Verify reconstruction
        eqs_reconstructed = f + g * inputs
        @test isequal(Symbolics.simplify.(eqs_reconstructed), rhs.(eqs))

        # Check dimensions
        @test length(f) == length(eqs)
        @test size(g) == (length(eqs), length(inputs))
    end

    # Test with Segway dynamics example
    @testset "Segway dynamics" begin
        # Segway parameters
        grav = 9.81
        R = 0.195
        M = 2 * 2.485
        Jc = 2 * 0.0559
        L = 0.169
        m = 44.798
        Jg = 3.836
        m0 = 52.710
        J0 = 5.108
        Km = 2 * 1.262
        bt = 2 * 1.225

        # Dynamics of Segway in Euler-Lagrange form
        Dq(q) = [m0 m*L*cos(q[2]); m*L*cos(q[2]) J0]
        function H(q, q̇)
            return SA[
                -m * L * sin(q[2]) * q̇[2] + bt * (q̇[1] - R * q̇[2]) / R,
                -m * grav * L * sin(q[2]) - bt * (q̇[1] - R * q̇[2])
            ]
        end
        B(q) = SA[Km / R, -Km]

        # Convert to control affine form
        function f_seg(x)
            q, q̇ = x[SA[1, 2]], x[SA[3, 4]]
            return [q̇; -Dq(q) \ H(q, q̇)]
        end
        function g_seg(x)
            q, q̇ = x[SA[1, 2]], x[SA[3, 4]]
            return [SA[0, 0]; Dq(q) \ B(q)]
        end

        # Trace dynamics symbolically
        @variables q1 q2 qd1 qd2 u
        x = [q1; q2; qd1; qd2]
        inputs = [u]
        eqs = D.(x) .~ f_seg(x) + g_seg(x) * u

        # Extract control-affine form
        fe, ge = input_affine_form(eqs, inputs)

        # Test reconstruction
        eqs2 = fe + ge * inputs
        diff = Symbolics.simplify.(eqs2 - rhs.(eqs), expand = true)

        # The difference should be zero or very close to zero symbolically
        # We test numerically since symbolic simplification might not be perfect
        f2, _ = build_function(fe, x, expression = false)
        g2, _ = build_function(ge, x, expression = false)

        for i in 1:10
            x_val = rand(length(x))
            @test f2(x_val) ≈ f_seg(x_val) rtol=1e-10
            @test g2(x_val) ≈ g_seg(x_val) rtol=1e-10
        end
    end

    # Test with multiple inputs
    @testset "Multiple inputs" begin
        @variables x1 x2 x3 u1 u2
        state = [x1, x2, x3]
        inputs = [u1, u2]

        eqs = D.(state) .~ [
            x2,
            x3,
            -x1 - 2*x2 - x3 + u1 + 3*u2
        ]

        f, g = input_affine_form(eqs, inputs)

        # Expected results
        f_expected = [x2, x3, -x1 - 2*x2 - x3]
        g_expected = [0 0; 0 0; 1 3]

        @test isequal(Symbolics.simplify.(f), Symbolics.simplify.(f_expected))

        # Test g matrix elements
        for i in 1:size(g, 1), j in 1:size(g, 2)

            @test isequal(Symbolics.simplify(g[i, j]), g_expected[i, j])
        end
    end

    # Test with nonlinear state dynamics
    @testset "Nonlinear state dynamics" begin
        @variables x1 x2 u
        state = [x1, x2]
        inputs = [u]

        eqs = D.(state) .~ [
            x2,
            -sin(x1) - x2 + u
        ]

        f, g = input_affine_form(eqs, inputs)

        # Expected results
        f_expected = [x2, -sin(x1) - x2]
        g_expected = reshape([0, 1], 2, 1)

        @test isequal(Symbolics.simplify.(f), Symbolics.simplify.(f_expected))
        @test isequal(g, g_expected)
    end
end
