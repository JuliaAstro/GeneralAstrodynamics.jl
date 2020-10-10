![RunTests](https://github.com/cadojo/Astrodynamics.jl/workflows/RunTests/badge.svg)
![Documentation](https://github.com/cadojo/Astrodynamics.jl/workflows/Documentation/badge.svg)
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

Check out the [Getting Started](https://cadojo.github.io/Astrodynamics.jl/stable/Overview/usage/#Getting-Started) for usage examples, and more detail about using this package. Some quick examples are shown below!

#### Two-body Problem

The `TwoBody` module handles Astrodynamics scenarios within the two-body problem. 

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

# Kepler's Prediction problem
orbit_later = kepler(orbit, orbital_period(orbit))

# Orbit propagation
sols = propagate(orbit, orbital_period(orbit))

# Plotting (with Plots.jl kwargs)
plot3d(sols; title="Plots.jl keywords work!", xlabel="Woo")

# Another true fact!
sols.step[end] ≈ orbit_later
```

#### NBody

The `NBody` module helps to solve the classical gravitational `NBody` problem. 

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

# Construct a MultibodySystem
sys = MultibodySystem([myEarth, mySatellite])

# Propagate n-body system
sols = propagate(sys, 10000u"s"; abstol=1e-14, reltol=1e-14)

# Plot n-body propagation results
plot3d(sols; title="Plots.jl keywords work!", xlabel="Woo")
```
