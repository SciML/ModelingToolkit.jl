using ModelingToolkit, OrdinaryDiffEq, NonlinearSolve, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

@connector Port begin
    p(t)
    dm(t)=0, [connect = Flow]
end

@connector Flange begin
    dx(t)=0
    f(t), [connect = Flow]
end

# Components ----
@mtkmodel Orifice begin
    @parameters begin
        Cₒ=2.7
        Aₒ=0.00094
        ρ₀=1000
        p′=0
    end
    @variables begin
        dm(t)=0
        p₁(t)=p′
        p₂(t)=p′
    end
    @components begin
        port₁ = Port(p=p′)
        port₂ = Port(p=p′)
    end
    begin
        u = dm/(ρ₀*Aₒ)
    end
    @equations begin
        dm ~ +port₁.dm
        dm ~ -port₂.dm
        p₁ ~ port₁.p
        p₂ ~ port₂.p

        p₁ - p₂ ~ (1/2)*ρ₀*u^2*Cₒ
    end
end

@mtkmodel Volume begin
    @parameters begin
        A=0.1
        ρ₀=1000
        β=2e9
        direction=+1
        p′
        x′
    end
    @variables begin
        p(t)=p′
        x(t)=x′
        dm(t)=0
        f(t)=p′ * A
        dx(t)=0
        r(t), [guess = 1000]
        dr(t), [guess = 1000]
    end
    @components begin
        port = Port(p=p′)
        flange = Flange(f=-p′ * A * direction)
    end
    @equations begin
        D(x) ~ dx
        D(r) ~ dr

        p ~ +port.p
        dm ~ +port.dm # mass is entering
        f ~ -flange.f * direction # force is leaving
        dx ~ flange.dx * direction

        r ~ ρ₀*(1 + p/β)
        dm ~ (r*dx*A) + (dr*x*A)
        f ~ p * A
    end
end

@mtkmodel Mass begin
    @parameters begin
        m = 100
        f′
    end
    @variables begin
        f(t)=f′
        x(t)=0
        dx(t)=0
        ẍ(t)=f′/m
    end
    @components begin
        flange = Flange(f=f′)
    end
    @equations begin
        D(x) ~ dx
        D(dx) ~ ẍ

        f ~ flange.f
        dx ~ flange.dx

        m*ẍ ~ f
    end
end

@mtkmodel Actuator begin
    @parameters begin
        p₁′
        p₂′
    end
    begin #constants
        x′=0.5
        A=0.1
    end
    @components begin
        port₁ = Port(p=p₁′)
        port₂ = Port(p=p₂′)
        vol₁ = Volume(p′=p₁′, x′=x′,  direction=-1)
        vol₂ = Volume(p′=p₂′, x′=x′,  direction=+1)
        mass = Mass(f′=(p₂′ - p₁′)*A)
        flange = Flange(f=0)
    end
    @equations begin
        connect(port₁, vol₁.port)
        connect(port₂, vol₂.port)
        connect(vol₁.flange, vol₂.flange, mass.flange, flange)
    end
end

@mtkmodel Source begin
    @parameters begin
        p′
    end
    @components begin
        port = Port(p=p′)
    end
    @equations begin
        port.p ~ p′
    end
end

@mtkmodel Damper begin
    @parameters begin
        c = 1000
    end
    @components begin
        flange = Flange(f=0)
    end
    @equations begin
        flange.f ~ c*flange.dx
    end
end

@mtkmodel System begin
    @components begin
        res₁ = Orifice(p′=300e5)
        res₂ = Orifice(p′=0)
        act = Actuator(p₁′=300e5, p₂′=0)
        src = Source(p′=300e5)
        snk = Source(p′=0)
        dmp = Damper()
    end
    @equations begin
        connect(src.port, res₁.port₁)
        connect(res₁.port₂, act.port₁)
        connect(act.port₂, res₂.port₁)
        connect(res₂.port₂, snk.port)
        connect(dmp.flange, act.flange)
    end
end

@mtkbuild sys = System()
initprob = ModelingToolkit.InitializationProblem(sys)
@test initprob isa NonlinearLeastSquaresProblem
@test length(initprob.u0) == 2
initsol = solve(initprob, reltol = 1e-12, abstol = 1e-12)
@test SciMLBase.successful_retcode(initsol)

@connector Flange begin
    dx(t), [guess = 0]
    f(t), [guess = 0, connect=Flow]
end

@mtkmodel Mass begin
    @parameters begin
        m = 100
    end
    @variables begin
        dx(t)
        f(t), [guess = 0]
    end
    @components begin
        flange = Flange()
    end
    @equations begin
        # connectors
        flange.dx ~ dx
        flange.f ~ -f
        
        # physics
        f ~ m*D(dx)
    end
end

@mtkmodel Damper begin
    @parameters begin
        d = 1
    end
    @variables begin
        dx(t), [guess = 0]
        f(t), [guess = 0]
    end
    @components begin
        flange = Flange()
    end
    @equations begin
        # connectors
        flange.dx ~ dx
        flange.f ~ -f
        
        # physics
        f ~ d*dx
    end
end

@mtkmodel MassDamperSystem begin
    @components begin
        mass = Mass(;dx=100,m=10)
        damper = Damper(;d=1)
    end
    @equations begin
        connect(mass.flange, damper.flange)
    end
end

@mtkbuild sys = MassDamperSystem()
initprob = ModelingToolkit.InitializationProblem(sys)
@test initprob isa NonlinearProblem
initsol = solve(initprob, reltol = 1e-12, abstol = 1e-12)
@test SciMLBase.successful_retcode(initsol)