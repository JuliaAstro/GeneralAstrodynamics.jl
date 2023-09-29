"""
A state vector for planar entry dynamics.
"""
mutable struct PlanarEntryState{F} <: AstrodynamicalState{F,4}
    γ::F
    v::F
    r::F
    θ::F

    PlanarEntryState{F}(::UndefInitializer) where {F} = new{F}()
    PlanarEntryState(::UndefInitializer) = PlanarEntryState{Float64}(undef)

    PlanarEntryState{F}(γ, v, r, θ) where {F} = new{F}(convert(F, γ), convert(F, v), convert(F, r), convert(F, θ))
    PlanarEntryState(γ, v, r, θ) = new{promote_type(typeof(γ), typeof(v), typeof(r), typeof(θ))}(γ, v, r, θ)
    PlanarEntryState(; γ=0, v=0, r=0, θ=0) = PlanarEntryState(γ, v, r, θ)
    PlanarEntryState{F}(; γ=0, v=0, r=0, θ=0) where {F} = PlanarEntryState{F}(γ, v, r, θ)
    PlanarEntryState(values::NamedTuple) =
        let (; γ, v, r, θ) = merge((; γ=0, v=0, r=0, θ=0), values)
            PlanarEntryState(γ, v, r, θ)
        end
    PlanarEntryState{F}(values::NamedTuple) where {F} =
        let (; γ, v, r, θ) = merge((; γ=0, v=0, r=0, θ=0), values)
            PlanarEntryState{F}(γ, v, r, θ)
        end
end

"""
A parameter vector for planar entry dynamics.
"""
struct PlanarEntryParameters{F} <: AstrodynamicalParameters{F,4}
    R::F
    P::F
    H::F
    m::F
    A::F
    C::F
    μ::F

    PlanarEntryParameters{F}(::UndefInitializer) where {F} = new{F}()
    PlanarEntryParameters(::UndefInitializer) = PlanarEntryParameters{Float64}(undef)

    PlanarEntryParameters{F}(R, P, H, m, A, C, μ) where {F} = new{F}(convert(F, γ), convert(F, v), convert(F, r), convert(F, θ))
    PlanarEntryParameters(R, P, H, m, A, C, μ) = new{promote_type(typeof(R), typeof(P), typeof(H), typeof(m), typeof(A), typeof(C), typeof(μ))}(R, P, H, m, A, C, μ)
    PlanarEntryParameters(; R=0, P=0, H=0, m=0, A=0, C=0, μ=0) = PlanarEntryParameters(R, P, H, m, A, C, μ)
    PlanarEntryParameters{F}(; R=0, P=0, H=0, m=0, A=0, C=0, μ=0) where {F} = PlanarEntryParameters{F}(R, P, H, m, A, C, μ)
    PlanarEntryParameters(values::NamedTuple) =
        let (; R, P, H, m, A, C, μ) = merge((; R=0, P=0, H=0, m=0, A=0, C=0, μ=0), values)
            PlanarEntryParameters(R, P, H, m, A, C, μ)
        end
    PlanarEntryParameters{F}(values::NamedTuple) where {F} =
        let (; R, P, H, m, A, C, μ) = merge((; R=0, P=0, H=0, m=0, A=0, C=0, μ=0), values)
            PlanarEntryParameters{F}(R, P, H, m, A, C, μ)
        end
end

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
@memoize function PlanarEntryFunction(; name=:PlanarEntry, kwargs...)
    defaults = (; jac=true)
    options = merge(defaults, kwargs)
    return ODEFunction{true,SciMLBase.FullSpecialize}(
        PlanarEntrySystem(; name=name);
        options...
    )
end
