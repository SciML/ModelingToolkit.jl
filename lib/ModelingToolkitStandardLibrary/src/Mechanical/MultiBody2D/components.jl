@mtkmodel Link begin
    @parameters begin
        m
        l
        I
        g
        x1_0 = 0.0
        y1_0 = 0.0
    end

    @variables begin
        (A(t)), [state_priority = 10]
        (dA(t)), [state_priority = 10]
        (ddA(t)), [state_priority = 10]

        fx1(t)
        fy1(t)

        fx2(t)
        fy2(t)

        x1(t) = x1_0
        dx1(t)

        y1(t) = y1_0
        dy1(t)

        x2(t) = l + x1_0
        dx2(t)

        y2(t)
        dy2(t)

        x_cm(t) = l / 2 + x1_0
        dx_cm(t)
        ddx_cm(t)

        y_cm(t)
        dy_cm(t)
        ddy_cm(t)
    end

    @components begin
        TX1 = Flange()
        TY1 = Flange()

        TX2 = Flange()
        TY2 = Flange()
    end

    @equations begin
        D(A) ~ dA
        D(dA) ~ ddA
        D(x1) ~ dx1
        D(y1) ~ dy1
        D(x2) ~ dx2
        D(y2) ~ dy2
        D(x_cm) ~ dx_cm
        D(dx_cm) ~ ddx_cm
        D(y_cm) ~ dy_cm
        D(dy_cm) ~ ddy_cm

        # x forces
        m * ddx_cm ~ fx1 + fx2

        # y forces
        m * ddy_cm ~ m * g + fy1 + fy2

        # torques
        I * ddA ~
        -fy1 * (x2 - x1) / 2 + fy2 * (x2 - x1) / 2 + fx1 * (y2 - y1) / 2 -
        fx2 * (y2 - y1) / 2

        # geometry
        x2 ~ l * cos(A) + x1
        y2 ~ l * sin(A) + y1
        x_cm ~ l * cos(A) / 2 + x1
        y_cm ~ l * sin(A) / 2 + y1
        TX1.f ~ fx1
        TX1.s ~ x1
        TY1.f ~ fy1
        TY1.s ~ y1
        TX2.f ~ fx2
        TX2.s ~ x2
        TY2.f ~ fy2
        TY2.s ~ y2
    end
end
