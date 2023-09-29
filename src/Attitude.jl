"""
A mutable state vector for attitude dynamics.
"""
mutable struct AttitudeState{F} <: AstrodynamicalState{F,7}
    q₁::F
    q₂::F
    q₃::F
    q₄::F
    ω₁::F
    ω₂::F
    ω₃::F


    AttitudeState{F}(::UndefInitializer) where {F} = new{F}()
    AttitudeState(::UndefInitializer) = AttitudeState{Float64}(undef)

    AttitudeState{F}(q₁, q₂, q₃, q₄, ω₁, ω₂, ω₃) where {F} = new{F}(q₁, q₂, q₃, q₄, ω₁, ω₂, ω₃)
    AttitudeState(q₁, q₂, q₃, q₄, ω₁, ω₂, ω₃) = new{promote_type(typeof(q₁), typeof(q₂), typeof(q₃), typeof(q₄), typeof(ω₁), typeof(ω₂), typeof(ω₃))}(q₁, q₂, q₃, q₄, ω₁, ω₂, ω₃)
    AttitudeState(; q₁=0.0, q₂=0.0, q₃=0.0, q₄=1.0, ω₁=0.0, ω₂=0.0, ω₃=0.0) = AttitudeState(q₁, q₂, q₃, q₄, ω₁, ω₂, ω₃)
    AttitudeState{F}(; q₁=zero(F), q₂=zero(F), q₃=zero(F), q₄=one(F), ω₁=zero(F), ω₂=zero(F), ω₃=zero(F)) where {F} = AttitudeState{F}(q₁, q₂, q₃, q₄, ω₁, ω₂, ω₃)
    AttitudeState{F}(values::NamedTuple) where {F} =
        let (; q₁, q₂, q₃, q₄, ω₁, ω₂, ω₃) = merge((; q₁=zero(F), q₂=zero(F), q₃=zero(F), q₄=one(F), ω₁=zero(F), ω₂=zero(F), ω₃=zero(F)), values)
            AttitudeState{F}(q₁, q₂, q₃, q₄, ω₁, ω₂, ω₃)
        end
    AttitudeState(values::NamedTuple) =
        let F = Float64, (; q₁, q₂, q₃, q₄, ω₁, ω₂, ω₃) = merge((; q₁=zero(F), q₂=zero(F), q₃=zero(F), q₄=one(F), ω₁=zero(F), ω₂=zero(F), ω₃=zero(F)), values)
            AttitudeState(q₁, q₂, q₃, q₄, ω₁, ω₂, ω₃)
        end
end

"""
A parameter vector for attitude dynamics.
"""
struct AttitudeParameters{F} <: AstrodynamicalParameters{F,15}
    J₁₁::F
    J₂₁::F
    J₃₁::F
    J₁₂::F
    J₂₂::F
    J₃₂::F
    J₁₃::F
    J₂₃::F
    J₃₃::F
    L₁::F
    L₂::F
    L₃::F
    f₁::F
    f₂::F
    f₃::F

    AttitudeParameters{F}(J₁₁, J₂₁, J₃₁, J₁₂, J₂₂, J₃₂, J₁₃, J₂₃, J₃₃, L₁, L₂, L₃, f₁, f₂, f₃) where {F} = new{F}(J₁₁, J₂₁, J₃₁, J₁₂, J₂₂, J₃₂, J₁₃, J₂₃, J₃₃, L₁, L₂, L₃, f₁, f₂, f₃)
    AttitudeParameters(J₁₁, J₂₁, J₃₁, J₁₂, J₂₂, J₃₂, J₁₃, J₂₃, J₃₃, L₁, L₂, L₃, f₁, f₂, f₃) = new{promote_type(typeof(J₁₁), typeof(J₂₁), typeof(J₃₁), typeof(J₁₂), typeof(J₂₂), typeof(J₃₂), typeof(J₁₃), typeof(J₂₃), typeof(J₃₃), typeof(L₁), typeof(L₂), typeof(L₃), typeof(f₁), typeof(f₂), typeof(f₃))}(J₁₁, J₂₁, J₃₁, J₁₂, J₂₂, J₃₂, J₁₃, J₂₃, J₃₃, L₁, L₂, L₃, f₁, f₂, f₂)
    AttitudeParameters(; J₁₁=one(F), J₂₁=zero(F), J₃₁=zero(F), J₁₂=zero(F), J₂₂=one(F), J₃₂=zero(F), J₁₃=zero(F), J₂₃=zero(F), J₃₃=one(F), L₁=zero(F), L₂=zero(F), L₃=zero(F), f₁=zero(F), f₂=zero(F), f₃=zero(F)) = AttitudeParameters(J₁₁, J₂₁, J₃₁, J₁₂, J₂₂, J₃₂, J₁₃, J₂₃, J₃₃, L₁, L₂, L₃, f₁, f₂, f₃)
    AttitudeParameters{F}(; J₁₁=one(F), J₂₁=zero(F), J₃₁=zero(F), J₁₂=zero(F), J₂₂=one(F), J₃₂=zero(F), J₁₃=zero(F), J₂₃=zero(F), J₃₃=one(F), L₁=zero(F), L₂=zero(F), L₃=zero(F), f₁=zero(F), f₂=zero(F), f₃=zero(F)) where {F} = AttitudeParameters{F}(J₁₁, J₂₁, J₃₁, J₁₂, J₂₂, J₃₂, J₁₃, J₂₃, J₃₃, L₁, L₂, L₃, f₁, f₂, f₃)
    AttitudeParameters(values::NamedTuple) =
        let (; J₁₁, J₂₁, J₃₁, J₁₂, J₂₂, J₃₂, J₁₃, J₂₃, J₃₃, L₁, L₂, L₃, f₁, f₂, f₃) = merge((; J₁₁=one(F), J₂₁=zero(F), J₃₁=zero(F), J₁₂=zero(F), J₂₂=one(F), J₃₂=zero(F), J₁₃=zero(F), J₂₃=zero(F), J₃₃=one(F), L₁=zero(F), L₂=zero(F), L₃=zero(F), f₁=zero(F), f₂=zero(F), f₃=zero(F)), values)
            AttitudeParameters(J₁₁, J₂₁, J₃₁, J₁₂, J₂₂, J₃₂, J₁₃, J₂₃, J₃₃, L₁, L₂, L₃, f₁, f₂, f₃)
        end
    AttitudeParameters{F}(values::NamedTuple) where {F} =
        let (; J₁₁, J₂₁, J₃₁, J₁₂, J₂₂, J₃₂, J₁₃, J₂₃, J₃₃, L₁, L₂, L₃, f₁, f₂, f₃) = merge((; J₁₁=one(F), J₂₁=zero(F), J₃₁=zero(F), J₁₂=zero(F), J₂₂=one(F), J₃₂=zero(F), J₁₃=zero(F), J₂₃=zero(F), J₃₃=one(F), L₁=zero(F), L₂=zero(F), L₃=zero(F), f₁=zero(F), f₂=zero(F), f₃=zero(F)), values)
            AttitudeParameters{F}(J₁₁, J₂₁, J₃₁, J₁₂, J₂₂, J₃₂, J₁₃, J₂₃, J₃₃, L₁, L₂, L₃, f₁, f₂, f₃)
        end
end

"""
A `ModelingToolkit.ODESystem` for atmospheric entry. Currently, only exponential atmosphere
models are provided! The output model is cached with `Memoize.jl`. Planet-specific
parameters default to Earth values.

The order of the states follows: `[q₁, q₂, q₃, q₄, ω₁, ω₂, ω₃]`.

The order of the parameters follows: `[]`

# Extended Help

This model describes how an object moves through an exponential atmosphere, above a
spherical planet.

## States

1. `q`: scalar-last attitude quaternion
2. `ω`: body rates (radians per second)

## Parameters

1. `J`: inertial matrix
2. `L`: TODO: what the hell is this?
3. `f`: torques on the vehicle body (Newton-meters)


### Usage

```julia
model = Attitude()
```
"""
@memoize function AttitudeSystem(; stm=false, name=:Attitude)

    @variables t
    @variables (q(t))[1:4] [description = "scalar-last attitude quaternion"]
    @variables (ω(t))[1:3] [description = "body rates in radians per second"]
    @parameters J[1:3, 1:3] [description = "moment of inertial matrix"]
    @parameters L[1:3] [description = "lever arm where input torque is applied"]
    @parameters f[1:3] [description = "input torque"]
    δ = Differential(t)

    q = collect(q)
    ω = collect(ω)
    J = collect(J)
    L = collect(L)
    f = collect(f)

    A = [
        0 ω[3] -ω[2] ω[1]
        -ω[3] 0 ω[1] ω[2]
        ω[2] -ω[1] 0 ω[3]
        -ω[1] -ω[2] -ω[3] 0
    ]

    ωx = [
        0 -ω[3] ω[2]
        ω[3] 0 -ω[1]
        -ω[2] ω[1] 0
    ]

    eqs = vcat(
        δ.(q) .~ (1 // 2) * (A * q),
        δ.(ω) .~ (-inv(J) * ωx * J * ω + inv(J) * L + f)
    )

    if stm
        @variables (Φ(t))[1:7, 1:7]
        Φ = Symbolics.scalarize(Φ)
        M = Symbolics.jacobian(map(el -> el.rhs, eqs), vcat(q, ω))

        LHS = map(δ, Φ)
        RHS = M * Φ

        eqs = vcat(eqs, [LHS[i] ~ RHS[i] for i in eachindex(LHS)])
    end

    if string(name) == "Attitude" && stm
        modelname = Symbol("AttitudeWithSTM")
    else
        modelname = name
    end

    if stm
        return ODESystem(eqs, t, vcat(q, ω, vec(Φ)), vcat(vec(J), L, f); name=name, defaults=Dict(vec(Φ .=> I(7))))
    else
        return ODESystem(eqs, t, vcat(q, ω), vcat(vec(J), L, f); name=name)
    end
end


"""
Returns an `ODEFunction` for spacecraft attitude dynamics.

# Extended Help

### Usage

The `stm` and `name` keyword arguments are passed to `Attitude`. All other keyword arguments
are passed directly to `SciMLBase.ODEFunction`.

```julia
f = AttitudeFunction()
let u = randn(7), p = randn(15), t = NaN # time invariant
    f(u, p, t)
end
```
"""
@memoize function AttitudeFunction(; stm=false, name=:Attitude, kwargs...)
    defaults = (; jac=true)
    options = merge(defaults, kwargs)
    return ODEFunction{true,SciMLBase.FullSpecialize}(
        AttitudeSystem(; stm=stm, name=name);
        options...
    )
end
