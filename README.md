[![Tests](https://github.com/cadojo/GeneralAstrodynamics.jl/workflows/Tests/badge.svg)](https://github.com/cadojo/GeneralAstrodynamics.jl/actions?query=workflow%3ATests)
[![Docs](https://github.com/cadojo/GeneralAstrodynamics.jl/workflows/Documentation/badge.svg)](https://cadojo.github.io/GeneralAstrodynamics.jl/dev)

# GeneralAstrodynamics.jl
Common astrodynamics calculations, with units!

## ⚠️ Under Construction ⚠️
This package is currently being __completely refactored__! This is good news, 
since the new implementation will do some cool things with Julia's type system.
The documentation and tests may frequently fail on the `main` branch on
this repo as I work out the kinks. No unstable versions will be pushed 
to the general registry. All `stable` documentation and release versions
are still stable and work.

Here's a sneak peak at some features to come! The plot below shows a family of Halo orbits about Earth-Moon L1. Upcoming
package features will take advantage of manifolds about Halo orbits like these to find low-cost transfer designs
for interplanetary missions! 🚀

#### Analytical and Numerical Halo Orbit Solvers
![analytical_v_numerical](https://user-images.githubusercontent.com/12131808/117874829-7e006080-b26f-11eb-97df-cefaba8e2087.png)

#### Earth-Moon Halo Orbit Family
![earth_moon_halo](https://user-images.githubusercontent.com/12131808/117874868-8a84b900-b26f-11eb-9ef2-54e6658a261e.png)

#### CR3BP Manifolds
![unstable_manifold](https://user-images.githubusercontent.com/12131808/117874975-a4be9700-b26f-11eb-8db0-80c90240156b.png)

## Features
* Restricted two-body problem equations, states, propagation, and plotting
* Restricted three-body problem equations, states, propagation, and iterative Halo orbit solvers
* N-body problem equations, states, propagation, and plotting
* A collection of fairly accurate planetary constants from our solar system (pulled from [SPICE](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/) kernals)

More to come! In the near term, additional features will include...
* Manifold-based transfer equations and states within the circular restricted three-body problem
* Hohmann-based transfer equations and states within the restricted two-body problem
* Zero-velocity curve plots for circular restricted three-body problem trajectories
* Stability analysis for circular restricted three-body problem states


## Usage

Check out the [Getting Started](https://cadojo.github.io/GeneralAstrodynamics.jl/stable/Overview/getting-started/#Getting-Started) documentation for code examples, and more detail about using this package. 
