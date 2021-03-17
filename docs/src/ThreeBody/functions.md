# `ThreeBody` Calculations

Common `ThreeBody` problem calculations are provided through functions. 

## Frame Representations

You can convert between the Inertial and Rotating (Synodic) reference frames through
the `inertial` and `synodic` functions.

```@docs
inertial
synodic
```

## Dimensionalization

Functions to nondimensionalize spacecraft states, and re-dimensionalize spacecraft
states are provided.

```@docs
time_scale_factor,
nondimensionalize_length,
nondimensionalize_velocity,
nondimensionalize_time,
nondimensionalize_mass_parameter,
nondimensionalize,
redimensionalize_length,
redimensionalize_velocity,
redimensionalize_time,
redimensionalize
```

## Halo Orbit Solvers
```@docs
halo_analytic,
halo
```

## Other Common Calculations

```@docs
potential_energy, 
jacobi_constant,
lagrange,
potential_energy_hessian, 
accel,
halo_numerical_tic!,
state_transition_dynamics,
nondimensional_radius,  
inertial, 
synodic
```