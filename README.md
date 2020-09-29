# Astrodynamics.jl
A simple Astrodynamics package, written with Julia!

## Motivation 

This package was created to learn more about Astrodynamics, and will be developed alongside a Graduate Astrodynamics course at the University of Maryland. The package [JuliaSpace/Astrodynamics.jl](https://github.com/JuliaSpace/Astrodynamics.jl) is much more thought out, and is more fully featured. Basically, if you're not me, use that one!

## Usage

#### Two-body Problem

The `TwoBody` module handles Astrodynamics scenarios within the two-body problem. You can make a `TwoBodyOrbit` by specifying a `CelestialBody` (currently only `Earth` and `Sun` are supported), and a Cartesian or Keplerian state.

```Julia
    # Cartesian state to TwoBodyOrbit
    r̅      = [0.0, 11681.0, 0.0]u"km"
    v̅      = [5.134, 4.226, 2.787]u"km/s"
    orbit1 = TwoBodyOrbit(r̅, v̅, Earth)

    # Keplerian state to TwoBodyOrbit
    e      =  0.3      * u"rad"
    a      =  15000    * u"km" + Earth.R
    i      =  10       * u"°"
    Ω      =  0        * u"°"
    ω      =  10       * u"°"
    ν      =  0        * u"°"
    orbit2 =  TwoBodyOrbit(e, a, i, Ω, ω, ν, Earth)

    # This is a true fact!
    orbit1 ≈ orbit2

    # For the rest of this section...
    orbit = orbit1
```

Now you can solve __Kepler's Prediction Problem__,  __propagate__ the satellite's trajectory over a specified intervol in time, and __plot__ the resultant trajectory with `Plots.jl`.

```Julia
# Kepler's Prediction problem
orbit_later = kepler(orbit, orbital_period(orbit))

# Orbit propagation
sols = propagate_twobody(orbit)

# Plotting (with Plots.jl kwargs)
twobody_plot3d(orbit; kwargs...)

# Another true fact!
sols.step[end] ≈ orbit_later
```

You may have noticed the `orbital_period` function. All common two-body problem equations have been included as functions with common arguments,`orbital_period(a, μ)`, and with `Astrodynamics.jl` structure arguments, `orbital_period(orbit)`. The current list of supported functions is shown below.

```Julia
# Look at ~all~ those functions
semimajor_axis
eccentricity
eccentricity_vector
inclination
true_anomoly, 
periapsis_radius
apoapsis_radius
periapsis_velocity
apoapsis_velocity,      
radius
radius_vector 
velocity
velocity_vector
orbital_period
time_since_periapsis
mean_motion
mean_motion_vector
semi_parameter
conic_anomoly
specific_angular_momentum
specific_angular_momentum_vector
specific_energy
```

Not sure how to use one of those? Check the docstrings in Julia's REPL!

```Julia
help?> eccentricity
search: eccentricity eccentricity_vector

  eccentricity(r̅, v̅, μ)
  eccentricity(orbit::TwoBodyOrbit)

  Returns orbital eccentricity, e.
```

#### NBody

The `NBody` module helps to solve the classical gravitational `NBody` problem. This is the baby version - point mass bodies, and no relativity. But it's still useful!

You can make your own `Body` by specifying an initial Cartesian state, and a mass.

```Julia
# It's MY Earth, and I want it now
r̅₁ = [0.0, 0.0, 0.0]u"km"
v̅₁ = [0.0, 0.0, 0.0]u"km/s"
m₁ = Earth.m
myEarth = Body(r̅₁, v̅₁, m₁)

# And we'll need a satellite...
r̅₂ = [0.0, 11681.0, 0.0]u"km"
v̅₂ = [5.134, 4.226, 2.787]u"km/s"
m₂ = 1000.0u"kg"
mySatellite = Body(r̅₂, v̅₂, m₂)
```

A `MultibodySystem` contains many `Bodies`.

```Julia
# Construct a MultibodySystem
sys = MultibodySystem([myEarth, mySatellite])
```

And you can __propagate__ a `MultibodySystem` through time to numerically find the final states for each `Body`! The package `DifferentialEquations.jl` is used for the numerical integration. For all __propagation__ functions in `Astrodynamics.jl`, you can specify `kwargs` as you would for a `DifferentialEquations.jl` `solve` call.

```Julia
sols = propagate_multibody(sys, 10000u"s"; abstol=1e-14, reltol=1e-14)
```

## More to come!

~ Joe



