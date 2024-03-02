#
# Restricted Two-body Problem models
#

@doc CartesianState
const R2BState = CartesianState

"""
A parameter vector for R2BP dynamics.
"""
Base.@kwdef struct R2BParameters{F} <: AstrodynamicalParameters{F,1}
    μ::F

    R2BParameters(μ) = new{typeof(μ)}(μ)
    R2BParameters{F}(μ) where {F} = new{F}(μ)
    R2BParameters(p::R2BParameters) = R2BParameters(p.μ)
    R2BParameters{F}(p::R2BParameters) where {F} = new{F}(p.μ)
end

Base.@pure paradigm(::R2BParameters) = "Restricted Two Body Dynamics"

"""
A `ModelingToolkit.ODESystem` for the Restricted Two-body Problem.

The order of the states follows: `[x, y, z, ẋ, ẏ, ż]`.

The order of the parameters follows: `[μ]`.

# Extended Help
The Restricted Two-body Problem is a simplified dynamical model
describing one small body (spacecraft, etc.) and one celestial
body. The gravity of the celestial body exhibits a force on the
small body. This model is commonly used as a simplification to
descibe our solar systems' planets orbiting our sun, or a
spacecraft orbiting Earth.

### Usage

```julia
model = R2BSystem()
```
"""
@memoize function R2BSystem(; stm=false, name=:R2B)

    @parameters t μ
    @variables x(t) y(t) z(t) ẋ(t) ẏ(t) ż(t)
    δ = Differential(t)
    r = [x, y, z]
    v = [ẋ, ẏ, ż]

    eqs = vcat(
        δ.(r) .~ v,
        δ.(v) .~ -μ .* (r ./ norm(r)^3)
    )

    if stm
        @variables (Φ(t))[1:6, 1:6] [description = "state transition matrix estimate"]
        Φ = Symbolics.scalarize(Φ)
        A = Symbolics.jacobian(map(el -> el.rhs, eqs), vcat(r, v))

        LHS = map(δ, Φ)
        RHS = map(simplify, A * Φ)

        eqs = vcat(eqs, [LHS[i] ~ RHS[i] for i in 1:length(LHS)])
    end

    if string(name) == "R2B" && stm
        modelname = Symbol("R2BWithSTM")
    else
        modelname = name
    end

    if stm
        return ODESystem(
            eqs, t, vcat(r, v, vec(Φ)), [μ];
            name=modelname,
            defaults=Dict(vec(Φ .=> I(6)))
        )
    else
        return ODESystem(
            eqs, t, vcat(r, v), [μ]; name=modelname
        )
    end
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
let u = randn(6), p = randn(1), t = 0
    f(u, p, t)
end
```
"""
@memoize function R2BFunction(; stm=false, name=:R2B, kwargs...)
    defaults = (; jac=true)
    options = merge(defaults, kwargs)
    sys = complete(R2BSystem(; stm=stm, name=name); split=false)
    return ODEFunction{true,SciMLBase.FullSpecialize}(
        sys, ModelingToolkit.unknowns(sys), ModelingToolkit.parameters(sys);
        options...
    )
end

"""
An `Orbit` which exists within R2BP dynamics.
"""
const R2BOrbit = Orbit{<:CartesianState,<:R2BParameters}

"""
Return an `ODEProblem` for the provided R2B system.
"""
R2BProblem(u0, tspan, p; kwargs...) = ODEProblem(R2BFunction(), u0, tspan, p; kwargs...)
R2BProblem(orbit::AstrodynamicalOrbit{<:CartesianState}, tspan::Union{<:Tuple,<:AbstractArray}; kwargs...) = ODEProblem(R2BFunction(), AstrodynamicalModels.state(orbit), tspan, AstrodynamicalModels.parameters(orbit); kwargs...)
R2BProblem(orbit::AstrodynamicalOrbit{<:CartesianState}, Δt; kwargs...) = R2BProblem(orbit, (zero(Δt), δt); kwargs...)