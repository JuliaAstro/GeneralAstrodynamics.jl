# Astrodynamics.jl
A simple Astrodynamics package, written with Julia!

## Motivation 

This package was created to learn more about Astrodynamics, and will be developed alongside a Graduate Astrodynamics course at the University of Maryland. The package [JuliaSpace/Astrodynamics.jl](https://github.com/JuliaSpace/Astrodynamics.jl) is much more thought out, and is more fully featured. Basically, if you're not me, use that one!

## Credits

\[1\] Vallado, David A. Fundamentals of astrodynamics and applications. Vol. 12. Springer Science & Business Media, 2001.
* All equations and algorithms within `Astrodynamics` are pulled from Vallado's _Fundamentals of Astrodynamics and Applications_, as well as course notes from ENAE 601 (Astrodynamics) at the University of Maryland.

\[2\] [JuliaAstro/AstroBase.jl](https://github.com/JuliaAstro/AstroBase.jl)
* `AstroBase` is referenced as a well thought-out Julia package structure example (I'm new to Julia!), as well as feature ideas.
* For example: my first shot at implementing `TwoBody` used separate structures for Cartesian two-body orbits, and Keplerian two-body orbits. `AstroBase` seems to keep these in one structure - that's way better! 
* Now my `Orbit` structure tracks both Cartesian and Keplerian representations for orbit conditions _in the background_. You provide a Cartesian or Keplerian representation to the `Orbit` constructor, and `TwoBody` handles the transformations behind the scenes.

## Usage

#### Units are Required!

Currently, `Astrodynamics.jl` `reexport`'s `Unitful`, `UnitfulAstro`, and `UnitfulAngles`. Units are required for all `TwoBody` and `NBody` computations. In the future, I'd like to make units optional.

#### Two-body Problem

The `TwoBody` module handles Astrodynamics scenarios within the two-body problem. You can make a `Orbit` by specifying a `CelestialBody` (Sun, Earth, Moon, Mars, etc.), and a Cartesian or Keplerian state.

```Julia
# Cartesian state to Orbit
r̅      = [0.0, 11681.0, 0.0]u"km"
v̅      = [5.134, 4.226, 2.787]u"km/s"
orbit1 = Orbit(r̅, v̅, Earth)

# Keplerian state to Orbit
e      =  0.3      * u"rad"
a      =  15000    * u"km" + Earth.R
i      =  10       * u"°"
Ω      =  0        * u"°"
ω      =  10       * u"°"
ν      =  0        * u"°"
orbit2 =  Orbit(e, a, i, Ω, ω, ν, Earth)

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
sols = propagate(orbit, orbital_period(orbit))

# Plotting (with Plots.jl kwargs)
plot3d(sols; title="Plots.jl keywords work!", xlabel="Woo")

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
apoapsis_velocity
radius 
velocity
orbital_period
time_since_periapsis
mean_motion
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
  eccentricity(orbit::Orbit)

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

A `MultibodySystem` contains an array of `Bodies`.

```Julia
# Construct a MultibodySystem
sys = MultibodySystem([myEarth, mySatellite])
```

And you can __propagate__ a `MultibodySystem` through time to numerically find the final states for each `Body`. The package `DifferentialEquations.jl` is used for the numerical integration. For all __propagation__ functions in `Astrodynamics.jl`, you can specify `kwargs` as you would for a `DifferentialEquations.jl` `solve` call.

```Julia
# Propagate n-body system
sols = propagate(sys, 10000u"s"; abstol=1e-14, reltol=1e-14)
```

As with a two-body `Orbit`, you can also plot each timestep in the n-body propagation.

```Julia
# Plot n-body propagation results
plot3d(sols; title="Plots.jl keywords work!", xlabel="Woo")
```

## More to come!

~ Joe



