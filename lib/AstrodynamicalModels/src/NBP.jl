#
# N-body models
#

"""
A `ModelingToolkit.System` for the Newtonian N-body Problem.

The order of the states follows:
`[x₁, y₁, z₁, ..., xₙ, yₙ, zₙ, ẋ₁, ẏ₁, ż₁, ..., ẋₙ, ẏₙ, żₙ]`.

The order of the parameters follows:
`[G, m₁, m₂, ..., mₙ]`.

!!! warning
    Be careful about specifying `stm=true` for systems with
    `N ≥ 3`! If state transition matrix dynamics are enabled,
    you can calculate the total number of system states with
    `N*6 + (N*6)^2`. Note that this increases exponentially as
    `N` grows! For `N == 6`, unless you're using parallelization,
    your computer may run for several hours.

# Extended Help
The N-body problem is a model which describes how `N` bodies
will move with respect to a common origin. This problem
typically involves many bodies which act due to one force:
electromagentism, gravity, etc. This model applies most
closely to many celestial bodies moving due to gravity.
That's about right for a model in a package called
`AstrodynamicalModels`!

### Usage

```julia
# One model for ALL the planets in our solar system 😎
model = NBSystem(9)
```
"""
@memoize function NBSystem(
    N::Int;
    stm = false,
    name = :NBP,
    kwargs...,
)

    N > 0 || throw(ArgumentError("`N` must be a number greater than zero!"))
    T = N * 6 + (N * 6)^2
    @parameters G m[1:N]
    @variables x(t)[1:N] y(t)[1:N] z(t)[1:N] ẋ(t)[1:N] ẏ(t)[1:N] ż(t)[1:N]

    r = [[x[i], y[i], z[i]] for i = 1:N]
    v = [[ẋ[i], ẏ[i], ż[i]] for i = 1:N]

    poseqs = reduce(vcat, [D.(r[i]) .~ v[i] for i ∈ 1:N])

    veleqs = reduce(
        vcat,
        [
            D.(v[i]) .~ sum((
                ((m[i] * m[j]) / norm(r[j] .- r[i])^3 * (r[j] .- r[i])) .*
                (j == N ? G / m[i] : 1) for j ∈ 1:N if j ≢ i
            )) for i ∈ 1:N
        ],
    )

    eqs = [poseqs; veleqs]

    if stm
        if N ≥ 3
            @warn """
            You requested state transition matrix dynamics for $N bodies.
            This will result in $(T) states! That may take
            a long time! Consider setting `stm=false`, and using
            `ModelingToolkit.calculate_jacobian` instead.
            """
        end

        @variables Φ(t)[1:length(eqs), 1:length(eqs)], [description = "state transition matrix estimate"]
        A = jacobian(map(el -> el.rhs, eqs), [r...; v...])
        eqs = [eqs; vec(scalarize(D(Φ) ~ A * Φ))]
    end

    if string(name) == "NBP" && stm
        modelname = Symbol("NBPWithSTM")
    else
        modelname = name
    end

    return System(eqs, t;
        name = modelname,
        kwargs...,
    )
end

"""
Returns an `ODEFunction` for NBP dynamics.
The order of states and parameters in the
`ODEFunction` arguments are equivalent to the
order of states and parameters for the system
produced with `NBP(N)`. As a general rule,
the order of the states follows: `[x₁, y₁, z₁,
..., xₙ, yₙ, zₙ, ẋ₁, ẏ₁, ż₁, ..., ẋₙ, ẏₙ, żₙ]`.

!!! note
    Unlike `R2BP` and `CR3BP`, `jac` is set to
    `false` by default. The number of states for `NBP`
    systems can be very large for relatively small
    numbers of bodies (`N`). Enabling `jac=true`
    by default would cause unnecessarily long
    waiting times for this @memoize function to return for
    `N ≥ 3` or so. If `N=2` and `stm=true`,
    setting `jac=true` could still result in several
    minutes of calculations, depending on the computer
    you're using.

!!! warning
    Be careful about specifying `stm=true` for systems with
    `N ≥ 3`! If state transition matrix dynamics are enabled,
    you can calculate the total number of system states with
    `N*6 + (N*6)^2`. Note that this increases exponentially as
    `N` grows! For `N == 6`, unless you're using parallelization,
    your computer may run for several hours.

# Extended Help

### Usage

The `stm`, and `name` keyword arguments
are passed to `NBP`. All other keyword arguments are passed
directly to `SciMLBase.ODEFunction`.

```julia
f = NBFunction(3; stm=false, name=:NBP, jac=false, sparse=false)
let u = randn(3*6), p = randn(1 + 3), t = 0
    f(u, p, t)
end
```
"""
@memoize function NBFunction(N::Int; stm = false, name = :R2BP, kwargs...)
    defaults = (; jac = false)
    options = merge(defaults, kwargs)
    if N ≥ 2 && stm && options.jac
        @warn """
        With state transition matrix dynamics appended to the state
        vector, and Jacobian computations requested, your system
        will require $((((N*6) + (N*6)^2) + ((N*6) + (N*6)^2)^2)) scalar
        calculations! Consider setting `jac=false`, `stm=false`, or both.
        """
    end
    sys = complete(NBSystem(N; stm, name); split = true)
    return ODEFunction{true,SciMLBase.FullSpecialize}(
        sys;
        options...,
    )
end
