"""
A `ModelingToolkit.ODESystem` for atmospheric entry. Currently, only exponential atmosphere
models are provided! The output model is cached with `Memoize.jl`. Planet-specific
parameters default to Earth values.

The order of the states follows: `[γ, v, r, θ]`.

The order of the parameters follows: `[R, P, H, m, A, C, μ]`

# Extended Help

This model describes how an object moves through an exponential atmosphere, above a
spherical planet.

### Usage

```julia
model = PlanarEntry()
```
"""
function PlanarEntry(; name=:PlanarEntry)

    @variables t

    @variables γ(t) [description = "flight path angle in degrees"]
    @variables v(t) [description = "airspeed in meters per second"]
    @variables r(t) [description = "polar distance relative to planet center in meters"]
    @variables θ(t) [description = "polar angle relative to planet horizontal in degrees"]
    @parameters R [description = "spherical planet radius"]
    @parameters P [description = "atmospheric density at sea level in kilograms per meter cubed"]
    @parameters H [description = "scale factor for exponential atmosphere in meters"]
    @parameters m [description = "entry vehicle mass in kilograms"]
    @parameters A [description = "entry vehicle surface area in square meters"]
    @parameters C [description = "lift to drag ratio in meters squared per quartic second"]
    @parameters μ [description = "planet mass parameter in meters per second cubed"]
    δ = Differential(t)

    β = m / (C * A)
    vc = √(μ / r)
    g₀ = μ / R^2
    g = g₀ * (R / r)^2
    h = r - R
    ρ = P * exp(-h / H)
    Dₘ = (ρ / 2) * v^2 / β
    Lₘ = C / Dₘ

    eqs = [
        δ(γ) ~ (one(v) / v) * (Lₘ - (1 - (v / vc)^2) * g * cos(γ)),
        δ(v) ~ -Dₘ - g * sin(γ),
        δ(r) ~ v * sin(γ),
        δ(θ) ~ (v / r) * cos(γ)
    ]

    model = ODESystem(
        eqs, t, [γ, v, r, θ], [R, P, H, m, A, C, μ]; name=name
    )

    return model

end


"""
Returns an `ODEFunction` for Planar Entry dynamics. Results are cached with `Memoize.jl`.

The order of the states follows: `[γ, v, r, θ]`.

The order of the parameters follows: `[R, P, H, m, A, C, μ]`

# Extended Help

### Usage

The `name` keyword argument is ]passed to `PlanarEntry`. All other keyword arguments are
passed directly to `SciMLBase.ODEFunction`.

```julia
f = PlanarEntryFunction()
let u = randn(4), p = randn(7), t = NaN # time invariant
    f(u, p, t)
end
```
"""
function PlanarEntryFunction(; name=:PlanarEntry, kwargs...)
    defaults = (; jac=true)
    options = merge(defaults, kwargs)
    return ODEFunction{true,SciMLBase.FullSpecialize}(
        PlanarEntry(; name=name);
        options...
    )
end
