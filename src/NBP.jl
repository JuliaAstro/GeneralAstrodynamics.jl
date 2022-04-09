#
# N-body models
# 

"""
A `ModelingToolkit.ODESystem` for the Newtonian N-body Problem. 

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
model = NBP(9) 
```
"""
@memoize function NBP(N::Int; stm=false, structural_simplify=true, name=:NBP)

    N > 0 || throw(ArgumentError("`N` must be a number greater than zero!")) 
    @parameters t G m[1:N]
    @variables x[1:N](t) y[1:N](t) z[1:N](t) ẋ[1:N](t) ẏ[1:N](t) ż[1:N](t)
    δ = Differential(t)

    r = [[x[i], y[i], z[i]] for i in 1:N]
    v = [[ẋ[i], ẏ[i], ż[i]] for i in 1:N]

    poseqs = reduce(vcat, [ 
        δ.(r[i]) .~ v[i] 
        for i ∈ 1:N 
    ])

    veleqs = reduce(vcat, [
        δ.(v[i]) .~ sum((
            ((m[i] * m[j]) / norm(r[j] .- r[i])^3 * (r[j] .- r[i])) .* (j == N ? G / m[i] : 1)
            for j ∈ 1:N if j ≢ i
        ))  for i ∈ 1:N 
    ])

    eqs = vcat(poseqs, veleqs)

    if stm 
        if N ≥ 3
            @warn """
            You requested state transition matrix dynamics for $N bodies.
            This will result in $(N*6 + (N*6)^2) states! That may take 
            a long time! Consider setting `stm=false`, and using
            `ModelingToolkit.calculate_jacobian` instead.
            """
        end

        @variables Φ[1:length(eqs),1:length(eqs)](t)
        Φ = Symbolics.scalarize(Φ)
        A = Symbolics.jacobian(map(el -> el.rhs, eqs), vcat(r...,v...))
    
        LHS = map(δ, Φ)
        RHS = map(simplify, A * Φ)

        eqs = vcat(eqs, [LHS[i] ~ RHS[i] for i in 1:length(LHS)])
    end

    if string(name) == "NBP" && stm 
        modelname = Symbol("NBPWithSTM")
    else
        modelname = name
    end

    sys = ODESystem(
        eqs, t, stm  ? vcat(r...,v...,Φ...) : vcat(r...,v...), vcat(G, m...); 
        name = modelname
    )
    return structural_simplify ? ModelingToolkit.structural_simplify(sys) : sys
end

"""
Returns an `ODEFunction` for NBP dynamics. 
Results are cached with `Memoize.jl`.
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
    waiting times for this function to return for 
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

The `stm`, `structural_simplify`, and `name` keyword arguments 
are passed to `NBP`. All other keyword arguments are passed
directly to `SciMLBase.ODEFunction`.

```julia
f = NBPFunction(3; stm=false, structural_simplify=true, name=:NBP, jac=false, sparse=false)
let u = randn(3*6), p = randn(1 + 3), t = 0
    f(u, p, t)
end
```
"""
@memoize function NBPFunction(N::Int; stm=false, structural_simplify=true, name=:R2BP, kwargs...)
    defaults = (; jac=false)
    options  = merge(defaults, kwargs)
    if N ≥ 2 && stm && options.jac
        @warn """
        With state transition matrix dynamics appended to the state 
        vector, and Jacobian compuations reqeusted, your system 
        will require $((((N*6) + (N*6)^2) + ((N*6) + (N*6)^2)^2)) scalar 
        calculations! Consider setting `jac=false`, `stm=false`, or both.
        """
    end
    return ODEFunction(
        NBP(N; stm=stm, structural_simplify=structural_simplify, name=name);
        options...
    )
end