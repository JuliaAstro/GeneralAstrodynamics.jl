[![Tests](https://github.com/cadojo/GeneralAstrodynamics.jl/workflows/Tests/badge.svg)](https://github.com/cadojo/GeneralAstrodynamics.jl/actions?query=workflow%3ATests)
[![Docs](https://github.com/cadojo/GeneralAstrodynamics.jl/workflows/Documentation/badge.svg)](https://cadojo.github.io/GeneralAstrodynamics.jl/)

# `GeneralAstrodynamics.jl`

_Common astrodynamics calculations, with hooks into the SciML ecosystem._

## JuliaCon Talk

Check out `GeneralAstrodynamics` in action at JuliaCon 2021! The talk
[_Going to Jupiter with Julia_](https://www.youtube.com/watch?v=WnvKaUsGv8w)
walks through a simple Jupiter mission design while gently introducing
astrodynamics, Julia, and `GeneralAstrodynamics`.

## Features

### Restricted Two-body Problem (R2BP)

- Structures for Cartesian and Keplerian states, and R2BP systems
- Functions which implement common R2BP equations
- Kepler and Lambert solvers
- Orbit propagation and plotting

### Circular Restricted Three-body Problem (CR3BP)

- Structures for dimensioned and normalized Cartesian states, and dimensioned
  and normalized CR3BP systems
- Functions which implement common CR3BP equations
- Analytical and iterative (numerical) Halo orbit solvers
- Unstable and stable Halo orbit manifold computation
- Orbit propagation and plotting
- Zero-velocity curve computation and plotting

### N-body Problem (NBP)

- This was implemented in a previous package version, and is currently being
  refactored

## Usage

```julia
using GeneralAstrodynamics

orbit = rand(R2BPOrbit)
trajectory = propagate(orbit, orbital_period(orbit))

furnsh(
    de440s(),                   # position and velocity data for nearby planets
    latest_leapseconds_tls(),   # timekeeping, parsing epochs
    gm_de440(),                 # mass parameters for major solar system bodies
    pck00011(),                 # physical properties of major solar system bodies
)

μ = reduced_mass(
  gm("earth"),
  gm("moon"),
)

orbit, T = let
  u, T = halo(μ, 2; amplitude=1e-2)

  CR3BPOrbit(CartesianState(u), CR3BParameters(μ)), T
end

trajectory = propagate(orbit, T)
```
