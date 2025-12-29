#
# Restricted Two-body Problem models
#

@doc CartesianState const R2BState = CartesianState

"""
A parameter vector for R2BP dynamics.
"""
Base.@kwdef struct R2BParameters{F} <: AstrodynamicalParameters{F,1}
    μ::F

    R2BParameters(μ) = new{typeof(μ)}(μ)
    R2BParameters{F}(μ) where {F} = new{F}(μ)
    R2BParameters(p::R2BParameters) = R2BParameters(p.μ)
    R2BParameters{F}(p::R2BParameters) where {F} = new{F}(p.μ)
    R2BParameters{F}(μ::Tuple) where {F} = R2BParameters{F}(μ...)
    R2BParameters(μ::Tuple) = R2BParameters(μ...)
end


system(::R2BParameters, args...; kwargs...) = R2BSystem(args...; kwargs...)
dynamics(::R2BParameters, args...; kwargs...) = R2BFunction(args...; kwargs...)
Base.@pure paradigm(::R2BParameters) = "Restricted Two Body Dynamics"

"""
A `ModelingToolkit.System` for the Restricted Two-body Problem.

The order of the states follows: `[x, y, z, ẋ, ẏ, ż]`.

The order of the parameters follows: `[μ]`.

# Extended Help
The Restricted Two-body Problem is a simplified dynamical model
describing one small body (spacecraft, etc.) and one celestial
body. The gravity of the celestial body exhibits a force on the
small body. This model is commonly used as a simplification to
describe our solar systems' planets orbiting our sun, or a
spacecraft orbiting Earth.

### Usage

```julia
model = R2BSystem()
```
"""
@memoize function R2BSystem(;
    stm = false,
    name = :R2B,
    kwargs...,
)

    @parameters μ
    @variables x(t) y(t) z(t) ẋ(t) ẏ(t) ż(t)
    r = [x, y, z]
    v = [ẋ, ẏ, ż]

    eqs = vcat(D.(r) .~ v, D.(v) .~ -μ .* (r ./ norm(r)^3))

    if stm
        @variables Φ(t)[1:6, 1:6], [description = "state transition matrix estimate"]

        A = Symbolics.jacobian(map(el -> el.rhs, eqs), vcat(r, v))
        eqs = vcat(eqs, vec(Symbolics.scalarize(D(Φ) ~ A * Φ)))
    end

    if string(name) == "R2B" && stm
        modelname = Symbol("R2BWithSTM")
    else
        modelname = name
    end

    return System(eqs, t;
        name = modelname,
        kwargs...,
    )
end

"""
Returns an `ODEFunction` for R2B dynamics.

The order of the states follows: `[x, y, z, ẋ, ẏ, ż]`.

The order of the parameters follows: `[μ]`.

# Extended Help

### Usage

The `stm`, and `name` keyword arguments
are passed to `R2B`. All other keyword arguments are passed
directly to `SciMLBase.ODEFunction`.

```julia
f = R2BFunction(; stm=false, name=:R2B, jac=true)
let u = randn(6), p = randn(), t = 0
    sys = f.sys
    u0 = get_u0(sys, [:x, :y, :z, :ẋ, :ẏ, :ż] .=> u) # Or get_u0(sys, ModelingToolkit.unknowns(sys) .=> u)
    p = get_p(sys, [:μ => p]) # Or get_p(sys, ModelingToolkit.parameters(sys) .=> p)
    f(u0, p, t)
end
```
"""
@memoize function R2BFunction(; stm = false, name = :R2B, kwargs...)
    defaults = (; jac = true)
    options = merge(defaults, kwargs)
    sys = complete(R2BSystem(; stm = stm, name = name); split = true)
    return ODEFunction{true,SciMLBase.FullSpecialize}(
        sys;
        options...,
    )
end

"""
An `Orbit` which exists within R2BP dynamics.
"""
const R2BOrbit = Orbit{<:CartesianState,<:R2BParameters}

"""
Return an `ODEProblem` for the provided R2B system.
"""
R2BProblem(op, tspan; kwargs...) = ODEProblem(R2BFunction().sys, op, tspan; kwargs...)
R2BProblem(
    orbit::AstrodynamicalOrbit{<:CartesianState},
    tspan::Union{<:Tuple,<:AbstractArray};
    kwargs...,
) = ODEProblem(
    R2BFunction().sys,
    AstrodynamicalModels.op(orbit),
    tspan;
    kwargs...,
)
R2BProblem(orbit::AstrodynamicalOrbit{<:CartesianState}, Δt; kwargs...) =
    R2BProblem(orbit, (zero(Δt), Δt); kwargs...)

# TODO: Deprecate old methods? https://github.com/SciML/ModelingToolkit.jl/blob/master/NEWS.md#new-problem-and-constructors
#R2BProblem(u0, tspan, p; kwargs...) = ODEProblem(R2BFunction(), u0, tspan, p; kwargs...)
#R2BProblem(
#    orbit::AstrodynamicalOrbit{<:CartesianState},
#    tspan::Union{<:Tuple,<:AbstractArray};
#    kwargs...,
#) = ODEProblem(
#    R2BFunction(),
#    AstrodynamicalModels.state(orbit),
#    tspan,
#    AstrodynamicalModels.parameters(orbit);
#    kwargs...,
#)
#R2BProblem(orbit::AstrodynamicalOrbit{<:CartesianState}, Δt; kwargs...) =
#    R2BProblem(orbit, (zero(Δt), δt); kwargs...)
