#
# Restricted Two-body Problem models
#

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
model = R2BP()
```
"""
function R2BP(; stm=false, name=:R2BP)

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

    if string(name) == "R2BP" && stm
        modelname = Symbol("R2BPWithSTM")
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
Returns an `ODEFunction` for R2BP dynamics.

The order of the states follows: `[x, y, z, ẋ, ẏ, ż]`.

The order of the parameters follows: `[μ]`.

# Extended Help

### Usage

The `stm`, and `name` keyword arguments
are passed to `R2BP`. All other keyword arguments are passed
directly to `SciMLBase.ODEFunction`.

```julia
f = R2BPFunction(; stm=false, name=:R2BP, jac=true)
let u = randn(6), p = randn(1), t = 0
    f(u, p, t)
end
```
"""
function R2BPFunction(; stm=false, name=:R2BP, kwargs...)
    defaults = (; jac=true)
    options = merge(defaults, kwargs)
    return ODEFunction(
        R2BP(; stm=stm, name=name);
        options...
    )
end
