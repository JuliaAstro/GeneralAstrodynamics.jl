#
# Circular Restricted Three-body Problem models
#

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
model = CR3BP(; stm=true)
```
"""
function CR3BP(; stm=false, name=:CR3BP)

    @parameters t μ
    @variables x(t) y(t) z(t) ẋ(t) ẏ(t) ż(t)
    δ = Differential(t)
    r = [x, y, z]
    v = [ẋ, ẏ, ż]

    eqs = vcat(
        δ.(r) .~ v,
        δ(ẋ) ~ x + 2ẏ - (μ * (x + μ - 1) * (sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3)) - ((x + μ) * (sqrt(y^2 + z^2 + (x + μ)^2)^-3) * (1 - μ)),
        δ(ẏ) ~ y - (2ẋ) - (y * (μ * (sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3) + (sqrt(y^2 + z^2 + (x + μ)^2)^-3) * (1 - μ))),
        δ(ż) ~ z * (-μ * (sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3) - ((sqrt(y^2 + z^2 + (x + μ)^2)^-3) * (1 - μ)))
    )

    if stm
        @variables (Φ(t))[1:6, 1:6] [description = "state transition matrix estimate"]
        Φ = Symbolics.scalarize(Φ)
        A = Symbolics.jacobian(map(el -> el.rhs, eqs), vcat(r, v))

        LHS = map(δ, Φ)
        RHS = map(simplify, A * Φ)

        eqs = vcat(eqs, [LHS[i] ~ RHS[i] for i in eachindex(LHS)])
    end

    if string(name) == "CR3BP" && stm
        modelname = Symbol("CR3BPWithSTM")
    else
        modelname = name
    end

    if stm
        return ODESystem(
            eqs, t, vcat(r, v, vec(Φ)), [μ];
            name=modelname,
            defaults=Dict(vec(Φ .=> map(Int, I(6))))
        )
    else
        return ODESystem(
            eqs, t, vcat(r, v), [μ]; name=modelname
        )
    end
end

"""
Returns an `ODEFunction` for CR3BP dynamics.

The order of the states follows: `[x, y, z, ẋ, ẏ, ż]`.

The order of the parameters follows: `[μ]`.

# Extended Help

### Usage

The `stm`, and `name` keyword arguments
are passed to `CR3BP`. All other keyword arguments are passed
directly to `SciMLBase.ODEFunction`.

```julia
f = CR3BPFunction(; stm=false, jac=true)
let u = randn(6), p = randn(1), t = 0
    f(u, p, t)
end
```
"""
function CR3BPFunction(; stm=false, name=:CR3BP, kwargs...)
    defaults = (; jac=true)
    options = merge(defaults, kwargs)
    return ODEFunction{true,SciMLBase.FullSpecialize}(
        CR3BP(; stm=stm, name=name);
        options...
    )
end
