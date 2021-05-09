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

<img width="713" alt="A plot of Halo orbits about Earth-Moon L1" src="https://user-images.githubusercontent.com/12131808/114203607-517ebf00-9926-11eb-8d77-16c7fc2b303a.png">

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

## Motivation 

This package aims to provide a simple interface for common astrodynamics problems. It was created to learn more about Astrodynamics, and will be developed alongside a Graduate Astrodynamics course at the University of Maryland. The packages [JuliaSpace/Astrodynamics.jl](https://github.com/JuliaSpace/Astrodynamics.jl) and [JuliaAstro/AstroBase.jl](https://github.com/JuliaAstro/AstroBase.jl) are more fully featured. I will continue adding features to this package, but for a more complete feature set, use the packages provided by [JuliaSpace](https://github.com/JuliaSpace) and [JuliaAstro](https://github.com/JuliaAstro).

## Credits

\[1\] Vallado, David A. Fundamentals of astrodynamics and applications. Vol. 12. Springer Science & Business Media, 2001.
* Many equations and algorithms within `GeneralAstrodynamics` are pulled from Vallado's _Fundamentals of Astrodynamics and Applications_, as well as course notes from ENAE 601 (Astrodynamics) at the University of Maryland.

\[2\] [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro.jl](https://github.com/JuliaAstro/UnitfulAstro.jl) are used for unit handling.

\[3\]

## Usage

Check out the [Getting Started](https://cadojo.github.io/GeneralAstrodynamics.jl/stable/Overview/getting-started/#Getting-Started) documentation for code examples, and more detail about using this package. 
