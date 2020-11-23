# Getting Started

## Installation

`UnitfulAstrodynamics` is included in Julia's General package registry.

```Julia
# In Julia's REPL
]add UnitfulAstrodynamics

# Or, with `Pkg`
import Pkg
Pkg.add("UnitfulAstrodynamics")
```

## Units are Required!

`UnitfulAstrodynamics.jl` uses `Reexport.jl` to expose `Unitful`, `UnitfulAstro`, and `UnitfulAngles`. Units are required for all `TwoBody`, `ThreeBody`, and `NBody` data structures. Functions often have non-unit equivalents -- 
check the docstrings!

## TwoBody

The `TwoBody` module handles Astrodynamics scenarios within the two-body problem. You can make a `Orbit` by specifying a `CelestialBody` (Sun, Earth, Moon, Mars, etc.), and a Cartesian or Keplerian state.

```Julia
# Cartesian state to Orbit
r      = [0.0, 11681.0, 0.0]u"km"
v      = [5.134, 4.226, 2.787]u"km/s"
orbit1 = Orbit(r, v, Earth)

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
```

Now you can solve __Kepler's Prediction Problem__,  __propagate__ the satellite's trajectory over a specified intervol in time, and __plot__ the resultant trajectory with `Plots.jl`.

```Julia
# Kepler's Prediction problem
orbit_later = kepler(orbit, orbital_period(orbit))

# Lambert's Proplem
v₁, v₂ = lambert(orbit.rᵢ, orbit_later.rᵢ, Earth.μ, orbital_period(orbit), :short)

# Orbit propagation
sols = propagate(orbit, orbital_period(orbit))

# Plotting (with Plots.jl kwargs)
plot(sols; title="Plots.jl keywords work!", xlabel="Woo")

# Another true fact!
sols.step[end] ≈ orbit_later
```

You may have noticed the `orbital_period` function. All common two-body problem equations have been included as functions with common arguments,`orbital_period(a, μ)`, and with `Astrodynamics.jl` structure arguments, `orbital_period(orbit)`. The current list of supported functions is described in [`TwoBody` Calculations](@ref).

Not sure how to use one of those helper functions? Check the docstrings in Julia's REPL!

```Julia
help?> eccentricity
search: eccentricity eccentricity_vector

  eccentricity(r̅, v̅, μ)
  eccentricity(orbit::Orbit)

  Returns orbital eccentricity, e.
```

## ThreeBody

The `ThreeBody` module helps to solve the Circular Restricted `ThreeBody` problem.

```julia
# Hardcode Gravity parameters for the Sun, 
# and the Earth-Moon System
μₛ = 1.32712440018e20u"m^3/s^2"
μₑ = 4.035032351966808e14u"m^3/s^2"

# Dimensional initial conditions for spacecraft
r = [2e9, 7000, 2000]u"km"
v = [0.001, 0.08, 0.02]u"km/s"
t = 500u"d"

# Construct nondimensional state
sys = ThreeBodySystem(AU, μₛ, μₑ, r, v, t);

# Propagate!
sols = propagate(sys)
```


## NBody

The `NBody` module helps to solve the classical gravitational `NBody` problem. This is the baby version - point mass bodies, and no relativity. But it's still useful!

You can make your own `Body` by specifying an initial Cartesian state, and a mass.

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
plot(sols; title="Plots.jl keywords work!", xlabel="Woo")
```
