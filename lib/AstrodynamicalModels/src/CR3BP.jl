#
# Circular Restricted Three-body Problem models
#

@doc CartesianState const CR3BState = CartesianState

"""
A parameter vector for CR3BP dynamics.
"""
Base.@kwdef struct CR3BParameters{F} <: AstrodynamicalParameters{F,1}
    μ::F

    CR3BParameters{F}(μ::Number) where {F} = new{F}(μ)
    CR3BParameters(μ::Number) = new{typeof(μ)}(μ)
    CR3BParameters{F}(μ::Tuple) where {F} = CR3BParameters{F}(μ...)
    CR3BParameters(μ::Tuple) = CR3BParameters(μ...)

end

system(::CR3BParameters, args...; kwargs...) = CR3BSystem(args...; kwargs...)
dynamics(::CR3BParameters, args...; kwargs...) = CR3BFunction(args...; kwargs...)
Base.@pure paradigm(::CR3BParameters) = "Circular Restricted Three Body Dynamics"


"""
A `ModelingToolkit.ODESystem` for the Circular Restricted Three-body Problem.

The order of the states follows: `[x, y, z, ẋ, ẏ, ż]`.

The order of the parameters follows: `[μ]`.

# Extended Help
The Circular Restricted Three-body Problem is a simplified dynamical model
describing one small body (spacecraft, etc.) and two celestial
bodies moving in a circle about their common center of mass.
This may seem like an arbitrary simplification, but this assumption
holds reasonably well for the Earth-Moon, Sun-Earth, and many other
systems in our solar system.

### Usage

```julia
model = CR3BSystem(; stm=true)
```
"""
@memoize function CR3BSystem(;
    stm = false,
    name = :CR3B,
    defaults = Pair{ModelingToolkit.Num,<:Number}[],
    kwargs...,
)
    @independent_variables t
    @parameters μ
    @variables x(t) y(t) z(t) ẋ(t) ẏ(t) ż(t)
    δ = Differential(t)
    r = [x, y, z]
    v = [ẋ, ẏ, ż]

    eqs = vcat(
        δ.(r) .~ v,
        δ(ẋ) ~
            x + 2ẏ - (μ * (x + μ - 1) * (sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3)) -
            ((x + μ) * (sqrt(y^2 + z^2 + (x + μ)^2)^-3) * (1 - μ)),
        δ(ẏ) ~
            y - (2ẋ) - (
                y * (
                    μ * (sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3) +
                    (sqrt(y^2 + z^2 + (x + μ)^2)^-3) * (1 - μ)
                )
            ),
        δ(ż) ~
            z * (
                -μ * (sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3) -
                ((sqrt(y^2 + z^2 + (x + μ)^2)^-3) * (1 - μ))
            ),
    )

    if stm
        @variables (Φ(t))[1:6, 1:6] [description = "state transition matrix estimate"]
        A = Symbolics.jacobian(map(el -> el.rhs, eqs), vcat(r, v))

        Φ = Symbolics.scalarize(Φ)
        LHS = δ.(Φ)
        RHS = A * Φ

        eqs = vcat(eqs, vec([LHS[i] ~ RHS[i] for i in eachindex(LHS)]))
    end

    if string(name) == "CR3B" && stm
        modelname = Symbol("CR3BWithSTM")
    else
        modelname = name
    end

    if stm
        defaults = vcat(defaults, vec(Φ .=> Float64.(I(6))))
        return ODESystem(
            eqs,
            t,
            vcat(r, v, vec(Φ)),
            [μ];
            name = modelname,
            defaults = defaults,
            kwargs...,
        )
    else
        return ODESystem(
            eqs,
            t,
            vcat(r, v),
            [μ];
            name = modelname,
            defaults = defaults,
            kwargs...,
        )
    end
end

"""
Returns an `ODEFunction` for CR3B dynamics.

The order of the states follows: `[μ]`.

The order of the parameters follows: `[μ]`.

# Extended Help

### Usage

The `stm`, and `name` keyword arguments
are passed to `CR3B`. All other keyword arguments are passed
directly to `SciMLBase.ODEFunction`.

```julia
f = CR3BFunction(; stm=false, jac=true)
let u = randn(6), p = randn(1), t = 0
    f(u, p, t)
end
```
"""
@memoize function CR3BFunction(; stm = false, name = :CR3B, kwargs...)
    defaults = (; jac = true)
    options = merge(defaults, kwargs)
    sys = complete(CR3BSystem(; stm = stm, name = name); split = false)
    return ODEFunction{true,SciMLBase.FullSpecialize}(
        sys;
        options...,
    )
end

"""
An `Orbit` which exists within CR3BP dynamics.
"""
const CR3BOrbit = Orbit{<:CR3BState,<:CR3BParameters}
AstrodynamicalModels.CR3BOrbit(state::AbstractVector, parameters::AbstractVector) =
    Orbit(CR3BState(state), CR3BParameters(parameters))
AstrodynamicalModels.CR3BOrbit(; state::AbstractVector, parameters::AbstractVector) =
    Orbit(CR3BState(state), CR3BParameters(parameters))

"""
Return an `ODEProblem` for the provided CR3B system.
"""
CR3BProblem(u0, tspan, p; kwargs...) = ODEProblem(CR3BFunction(), u0, tspan, p; kwargs...)
CR3BProblem(orbit::AstrodynamicalOrbit, tspan::Union{<:Tuple,<:AbstractArray}; kwargs...) =
    ODEProblem(
        CR3BFunction(),
        AstrodynamicalModels.state(orbit),
        tspan,
        AstrodynamicalModels.parameters(orbit);
        kwargs...,
    )
CR3BProblem(orbit::AstrodynamicalOrbit, Δt; kwargs...) =
    CR3BProblem(orbit, (zero(Δt), δt); kwargs...)
