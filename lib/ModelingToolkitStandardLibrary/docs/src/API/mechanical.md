# ModelingToolkit Standard Library: Mechanical Components

```@contents
Pages = ["mechanical.md"]
Depth = 3
```

## Index

```@index
Pages = ["mechanical.md"]
```

## Rotational Components

```@meta
CurrentModule = ModelingToolkitStandardLibrary.Mechanical.Rotational
```

### Rotational Utils

```@docs
Flange
Support
PartialCompliantWithRelativeStates
PartialElementaryOneFlangeAndSupport2
PartialElementaryTwoFlangesAndSupport2
PartialCompliant
```

### Rotational Core Components

```@docs
Fixed
Inertia
Spring
Damper
SpringDamper
IdealGear
RotationalFriction
```

### Rotational Sources

```@docs
Torque
Speed
Position
```

### Rotational Sensors

```@docs
AngleSensor
SpeedSensor
TorqueSensor
RelSpeedSensor
```

## Translational Components

```@meta
CurrentModule = ModelingToolkitStandardLibrary.Mechanical.Translational
```

### Translational Utils

```@docs
MechanicalPort
```

### Translational Core Components

```@docs
Mass
Spring
Damper
Fixed
```

### Translational Sources

```@docs
Force
Position
Velocity
Acceleration
```

### Translational Sensors

```@docs
ForceSensor
PositionSensor
AccelerationSensor
```
