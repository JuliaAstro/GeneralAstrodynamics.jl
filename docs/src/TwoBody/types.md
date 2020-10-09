## `TwoBody` Data Structures

The `TwoBody` module contains two key structures: `Orbit` and `CelestialBody`.

An `Orbit` is the core structure for `TwoBody` calculations. It contains an orbital state (both Cartesian and the equivalent Keplerian representation), and a central body.

The central body within `Orbit` is of type `CelestialBody`. Common bodies in our solar system have been added for convenience, as described in [Default `CelestialBodies`](@ref), but you can also make your own.


```@docs
Orbit
CelestialBody
```

## Abstract Types and Pre-defined Parameters

The first section in `TwoBody` documentation described the core [`TwoBody` Data Structures](@ref). Each data structure has an abstract parent type. All `Orbit` structures extend `TwoBodySystem`. In addition, all `Orbit` structures are paremeterized by their conic section, which is of type `AbstractConic`. All conic sections are pre-defined structures: `Circular`, `Elliptical`, `Parabolic`, `Hyperbolic`, and the `Invalid` conic is used to describe invalid orbital states (such as providing a `NaN` value to an `Orbit` constructor).

```@docs
TwoBodySystem
AbstractConic
Circular
Elliptical
Parabolic
Hyperbolic
Invalid
InvalidOrbit
```

## Default `CelestialBodies`

For convenence, the following common bodies in _our_ solar system have already been defined!

* `Sun`
* `Mercury` 
* `Venus `
* `Earth` 
* `Moon`
* `Luna`
* `Mars`
* `Jupiter` 
* `Saturn`
* `Uranus`
* `Neptune`
* `Pluto`
