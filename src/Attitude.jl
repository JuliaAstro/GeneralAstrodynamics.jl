"""
A `ModelingToolkit.ODESystem` for atmospheric entry. Currently, only exponential atmosphere
models are provided! The output model is cached with `Memoize.jl`. Planet-specific
parameters default to Earth values.

The order of the states follows: `[γ, v, r, θ]`.

The order of the parameters follows: `[R, P, H, m, A, C, μ]`

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
function Attitude(; stm=false, name=:Attitude)

    @variables t
    @variables (q(t))[1:4] [description = "scalar-last attitude quaternion"]
    @variables (ω(t))[1:3] [description = "body rates in radians per second"]
    @parameters J[1:3, 1:3] [description = "moment of inertial matrix"]
    @parameters L[1:3] [description = "lever arm where input torque is applied"]
    @parameters u[1:3] [description = "input torque"]
    δ = Differential(t)

    q = collect(q)
    ω = collect(ω)
    J = collect(J)
    L = collect(L)
    u = collect(u)

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
        δ.(ω) .~ (-inv(J) * ωx * J * ω + inv(J) * L + u)
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
        return ODESystem(eqs, t, vcat(q, ω, vec(Φ)), vcat(vec(J), L, u); name=name, defaults=Dict(vec(Φ .=> I(7))))
    else
        return ODESystem(eqs, t, vcat(q, ω), vcat(vec(J), L, u); name=name)
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
function AttitudeFunction(; stm=false, name=:Attitude, kwargs...)
    defaults = (; jac=true)
    options = merge(defaults, kwargs)
    return ODEFunction{true,SciMLBase.FullSpecialize}(
        Attitude(; stm=stm, name=name);
        options...
    )
end
