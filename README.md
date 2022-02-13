[![Tests](https://github.com/cadojo/GeneralAstrodynamics.jl/workflows/Tests/badge.svg)](https://github.com/cadojo/GeneralAstrodynamics.jl/actions?query=workflow%3ATests)
[![Docs](https://github.com/cadojo/GeneralAstrodynamics.jl/workflows/Documentation/badge.svg)](https://cadojo.github.io/GeneralAstrodynamics.jl/dev)

# GeneralAstrodynamics.jl
_Common astrodynamics calculations, with units!_

## JuliaCon Talk

Check out `GeneralAstrodynamics` in action at JuliaCon 2021! The talk [_Going to Jupiter with Julia_](https://www.youtube.com/watch?v=WnvKaUsGv8w) walks through a simple Jupiter mission design while gently introducing astrodynamics, Julia, and `GeneralAstrodynamics`.

## Features

### Restricted Two-body Problem (R2BP)
* Structures for Cartesian and Keplerian states, and R2BP systems
* Functions which implement common R2BP equations
* Kepler and Lambert solvers
* Orbit propagation and plotting

### Circular Restricted Three-body Problem (CR3BP)
* Structures for dimensioned and normalized Cartesian states, and dimensioned and normalized CR3BP systems
* Functions which implement common CR3BP equations
* Analytical and iterative (numerical) Halo orbit solvers
* Unstable and stable Halo orbit manifold computation
* Orbit propagation and plotting
* Zero-velocity curve computation and plotting

### N-body Problem (NBP)
* This was implemented in a previous package version, and is currently being refactored

## Usage

Some quick examples are below!

```julia
# Installation
import Pkg
Pkg.add("GeneralAstrodynamics") # or julia> ]install GeneralAstrodynamics

# Loading
using GeneralAstrodynamics, Unitful

# Construct a R2BP orbit (massless spacecraft 
# moving due to the gravity of one planet)
orbit = let e = 0.4, a = 10_000, i = Ω = ω = ν = 0, planet = Earth
    orbitalstate = KeplerianState(e, a, i, Ω, ω, ν)
    Orbit(orbitalstate, Earth)
end

# Alternatively, use a `CartesianState`
orbit = Orbit(
    CartesianState(randn(6)), # random state vector, [r..., v...]
    Earth
)

# Construct a CR3BP orbit (massless spacecraft moving
# due to the gravity of two planets, both of which
# move in a circle about their common center of mass)
orbit = Orbit(
    CartesianState(randn(6)), # random state vector (again!)
    SunEarth
)

# Propagate any orbit in time (after `using DifferentialEquations`)
using DifferentialEquations
trajectory = propagate(orbit, 10u"d") # unitful times are convenient here!

# Constract a periodic orbit within CR3BP dynamics (Halo orbit),
# and the orbital period `T` (also requires `DifferentialEquations`)
orbit, T = halo(SunEarth; L=1, Az=75_000u"km")

# Construct a manifold which converges to (stable), or 
# diverges from (unstable) the Halo orbit
superslide = manifold(orbit, T; duration=2T, eps=-1e8, direction=Val{:stable})

# Plot any `Trajectory` or `Manifold` (after `using Plots`)
using Plots
plot(trajectory; title="R2BP Trajectory")
plot(propagate(orbit, T); vars=:XY, label="Halo Orbit", aspect_ratio=1)
plot(superslide; vars=:XY, title="Stable Manifold near Earth")
```

In the coming months, the [Getting Started](https://cadojo.github.io/GeneralAstrodynamics.jl/dev/) page will have code examples, and other documentation for fundamental astrodynamics concepts, and `GeneralAstrodynamics` usage. Stay tuned and/or submit pull requests!
