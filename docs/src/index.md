# `AstrodynamicalSolvers.jl`

_Common solvers within orbital mechanics and astrodynamics._

## Installation

```julia
pkg> add AstrodynamicalSolvers
```

## Getting Started

This package currently provides periodic orbit, and manifold computations within 
Circular Restricted Three Body Problem dynamics.

### Periodic Orbits

This package contains differential correctors, and helpful wrapper functions, for 
finding periodic orbits within Circular Restricted Three Body Problem dynamics.

```@example usage
using AstrodynamicalSolvers
using AstrodynamicalModels
using OrdinaryDiffEq
using Plots

μ = 0.012150584395829193

planar = let
    u, T = halo(μ, 1) # lyapunov (planar) orbit
    u = [u.x, 0, 0, 0, u.ẏ, 0]
    problem = ODEProblem(CR3BFunction(), u, (0, T), (μ,))
    solution = solve(problem, Vern9(), reltol=1e-14, abstol=1e-14)
    plot(solution, idxs=(:x,:y,:z), title = "Lyapunov Orbit", label=:none, size=(1600,900), dpi=400, aspect_ratio=1)
end

extraplanar = let
    u, T = halo(μ, 2; amplitude=0.01) # halo (non-planar) orbit
    u = [u.x, 0, u.z, 0, u.ẏ, 0]
    problem = ODEProblem(CR3BFunction(), u, (0, T), (μ,))
    solution = solve(problem, Vern9(), reltol=1e-14, abstol=1e-14)
    plot(solution, idxs=(:x,:y,:z), title = "Halo Orbit", label=:none, size=(1600,900), dpi=400, aspect_ratio=1)
end

plot(planar, extraplanar, layout=(1,2))
```

### Manifold Computations

Manifold computations, provided by `AstrodynamicalCalculations.jl`, can perturb 
halo orbits onto their unstable or stable manifolds.

```@example
using AstrodynamicalSolvers
using AstrodynamicalCalculations
using AstrodynamicalModels
using OrdinaryDiffEq
using LinearAlgebra
using Plots

μ = 0.012150584395829193

unstable = let
    u, T = halo(μ, 1; amplitude=0.005)

    u = [u.x, 0, u.z, 0, u.ẏ, 0]
    Φ = monodromy(u, μ, T)

    ics = let
        problem = ODEProblem(CR3BFunction(stm=true), vcat(u, vec(I(6))), (0, T), (μ,))
        solution = solve(problem, Vern9(), reltol=1e-12, abstol=1e-12, saveat=(T / 10))

        solution.u
    end

    perturbations = [
        diverge(ic[1:6], reshape(ic[7:end], 6, 6), Φ; eps=-1e-7)
        for ic in ics
    ]

    problem = EnsembleProblem(
        ODEProblem(CR3BFunction(), u, (0.0, 2T), (μ,)),
        prob_func=(prob, i, repeat) -> remake(prob; u0=perturbations[i]),
    )

    solution = solve(problem, Vern9(), trajectories=length(perturbations), reltol=1e-14, abstol=1e-14)
end

stable = let
    u, T = halo(μ, 2; amplitude=0.005)

    u = [u.x, 0, u.z, 0, u.ẏ, 0]
    Φ = monodromy(u, μ, T)

    ics = let
        problem = ODEProblem(CR3BFunction(stm=true), vcat(u, vec(I(6))), (0, T), (μ,))
        solution = solve(problem, Vern9(), reltol=1e-12, abstol=1e-12, saveat=(T / 10))

        solution.u
    end
    
    perturbations = [
        converge(ic[1:6], reshape(ic[7:end], 6, 6), Φ; eps=1e-7)
        for ic in ics
    ]

    problem = EnsembleProblem(
        ODEProblem(CR3BFunction(), u, (0.0, -2.1T), (μ,)),
        prob_func=(prob, i, repeat) -> remake(prob; u0=perturbations[i]),
    )

    solution = solve(problem, Vern9(), trajectories=length(perturbations), reltol=1e-14, abstol=1e-14)
end

figure = plot(; 
    aspect_ratio = 1.0,
    background = :transparent,
    grid = true,
    title = "Unstable and Stable Invariant Manifolds",
    size = (1600,900),
    dpi = 400,
)

plot!(figure, unstable, idxs=(:x, :y), aspect_ratio=1, label=:none, palette=:blues)
plot!(figure, stable, idxs=(:x, :y), aspect_ratio=1, label=:none, palette=:blues)
scatter!(figure, [1-μ], [0], label="Moon", xlabel="X (Earth-Moon Distance)", ylabel="Y (Earth-Moon Distance)", marker=:x, color=:black, markersize=10,)

figure # hide
```
