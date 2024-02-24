using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ModelingToolkit

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
model = PlanarEntrySystem()
```
"""
@memoize function PlanarEntrySystem(; name=:PlanarEntry)

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

model = complete(PlanarEntrySystem())
field = ODEFunction(model)

x = collect(0.1:0.1:0.6)
p = collect(0.1:0.1:0.7)
t = NaN

@show field(x, p, t)

open(joinpath(@__DIR__, "entry.jl"), "w") do file
    func = string(generate_function(model)[1])
    write(file, func)
end