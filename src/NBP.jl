#
# N-body models
# 

"""
A `ModelingToolkit.ODESystem` for the Newtonian N-body Problem. 

# Extended Help
The N-body problem is a model which describes how `N` bodies 
will move with respect to a common origin. This problem 
typically involves many bodies which act due to one force:
electromagentism, gravity, etc. This model applies most 
closely to many celestial bodies moving due to gravity.
That's about right for a model in a package called
`AstrodynamicalModels`!
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
        if N > 5
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

# Extended Help

### Usage

The `stm`, `structural_simplify`, and `name` keyword arguments 
are passed to `NBP`. All other keyword arguments are passed
directly to `SciMLBase.ODEFunction`.

```julia
f = NBPFunction(; stm=false, structural_simplify=true, name=:NBP, jac=true, sparse=false)
```
"""
@memoize function NBPFunction(N::Int; stm=false, structural_simplify=true, name=:R2BP, kwargs...)
    defaults = (; jac=true, sparse=false)
    options  = merge(defaults, kwargs)
    return ODEFunction(
        NBP(N; stm=stm, structural_simplify=structural_simplify, name=name);
        options...
    )
end