[![Tests](https://github.com/cadojo/UnitfulAstrodynamics.jl/workflows/Tests/badge.svg)](https://github.com/cadojo/UnitfulAstrodynamics.jl/actions?query=workflow%3ATests)
[![Docs](https://github.com/cadojo/UnitfulAstrodynamics.jl/workflows/Documentation/badge.svg)](https://cadojo.github.io/UnitfulAstrodynamics.jl/stable)

[![version](https://juliahub.com/docs/UnitfulAstrodynamics/version.svg)](https://juliahub.com/ui/Packages/UnitfulAstrodynamics/uJGLZ)

# UnitfulAstrodynamics.jl
Common astrodynamics calculations with units!

## Motivation 

This package aims to provide a simple interface for common astrodynamics problems. It was created to learn more about Astrodynamics, and will be developed alongside a Graduate Astrodynamics course at the University of Maryland. The packages [JuliaSpace/Astrodynamics.jl](https://github.com/JuliaSpace/Astrodynamics.jl) and [JuliaAstro/AstroBase.jl](https://github.com/JuliaAstro/AstroBase.jl) are more fully featured. I will continue adding features to this package, but for a more complete feature set, use the packages provided by [JuliaSpace](https://github.com/JuliaSpace) and [JuliaAstro](https://github.com/JuliaAstro).

## Credits

\[1\] Vallado, David A. Fundamentals of astrodynamics and applications. Vol. 12. Springer Science & Business Media, 2001.
* All equations and algorithms within `UnitfulAstrodynamics` are pulled from Vallado's _Fundamentals of Astrodynamics and Applications_, as well as course notes from ENAE 601 (Astrodynamics) at the University of Maryland.

\[2\] [JuliaAstro/AstroBase.jl](https://github.com/JuliaAstro/AstroBase.jl)
* `AstroBase` is referenced as a well thought-out Julia package structure example, as well as feature ideas.

\[3\] [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro.jl](https://github.com/JuliaAstro/UnitfulAstro.jl) are used for unit handling.

## Usage

Check out the [Getting Started](https://cadojo.github.io/UnitfulAstrodynamics.jl/stable/#Getting-Started) documentation for code examples, and more detail about using this package. Some quick examples are shown below!

#### Two-body Problem

The `TwoBody` module handles Astrodynamics scenarios within the two-body problem. 

```Julia
# Cartesian state to Orbit
rᵢ = [0.0, 11681.0, 0.0] * u"km"
vᵢ = [5.134, 4.226, 2.787] * u"km/s"
orbit1 = Orbit(rᵢ, vᵢ, Earth)

# Keplerian state to Orbit
e      =  0.3      
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
plot(sols; title="Plots.jl keywords work!", xlabel="Woo")

# Another true fact!
sols.step[end] ≈ orbit_later
```

#### NBody

The `NBody` module helps to solve the classical gravitational `NBody` problem. 

```Julia
# It's MY Earth, and I want it now
r₁ = [0.0, 0.0, 0.0]u"km"
v₁ = [0.0, 0.0, 0.0]u"km/s"
m₁ = Earth.m
myEarth = Body(r₁, v₁, m₁)

# And we'll need a satellite...
r₂ = [0.0, 11681.0, 0.0]u"km"
v₂ = [5.134, 4.226, 2.787]u"km/s"
m₂ = 1000.0u"kg"
mySatellite = Body(r₂, v₂, m₂)

# Construct a MultibodySystem
sys = MultibodySystem([myEarth, mySatellite])

# Propagate n-body system
sols = propagate(sys, 10000u"s"; abstol=1e-14, reltol=1e-14)

# Plot n-body propagation results
plot(sols; title="Plots.jl keywords work!", xlabel="Woo")
```
