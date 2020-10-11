## Package Overview

`UnitfulAstrodynamics.jl` is a `Unitful` Astrodynamics library, which includes `TwoBody` and `NBody` problem calculations, as well as orbit propagation and plotting!

## Motivation 

This package was created to learn more about Astrodynamics, and will be developed alongside a Graduate Astrodynamics course at the University of Maryland. The package [JuliaSpace/Astrodynamics.jl](https://github.com/JuliaSpace/Astrodynamics.jl) is much more thought out, and is more fully featured. Basically, if you're not me, use that one!

## Credits

\[1\] Vallado, David A. Fundamentals of astrodynamics and applications. Vol. 12. Springer Science & Business Media, 2001.
* All equations and algorithms within `Astrodynamics` are pulled from Vallado's _Fundamentals of Astrodynamics and Applications_, as well as course notes from ENAE 601 (Astrodynamics) at the University of Maryland.

\[2\] [JuliaAstro/AstroBase.jl](https://github.com/JuliaAstro/AstroBase.jl)
* `AstroBase` is referenced as a well thought-out Julia package structure example (I'm new to Julia!), as well as feature ideas.
* For example: my first shot at implementing `TwoBody` used separate structures for Cartesian two-body orbits, and Keplerian two-body orbits. `AstroBase` seems to keep these in one structure - that's way better! 
* Now my `Orbit` structure tracks both Cartesian and Keplerian representations for orbit conditions _in the background_. You provide a Cartesian or Keplerian representation to the `Orbit` constructor, and `TwoBody` handles the transformations behind the scenes.
