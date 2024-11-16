"""
Wrappers around SciML differential equation solvers for fast and convenient 
orbit propagation.

# Extended Help

## Exports
$(EXPORTS)

## Imports
$(IMPORTS)
"""
module Propagation

using DocStringExtensions

@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(SIGNATURES)

                                         $(DOCSTRING)
                                         """

@template (TYPES, CONSTANTS) = """
                               $(TYPEDEF)

                               $(DOCSTRING)
                               """

using AstrodynamicalCalculations
using AstrodynamicalModels
using ModelingToolkit, OrdinaryDiffEqVerner, SciMLBase
using StaticArrays

export propagate, propagate!, monodromy, convergent_manifold, divergent_manifold

function SciMLBase.ODEProblem(
    orbit::AstrodynamicalModels.AstrodynamicalOrbit,
    Δt;
    stm = false,
    kwargs...,
)
    f = dynamics(orbit, stm = stm)
    u = AstrodynamicalModels.state(orbit)

    if stm
        u = Vector(vcat(u, vec(AstrodynamicalModels.CartesianSTM())))
    end

    p = AstrodynamicalModels.parameters(orbit)
    tspan = (Δt isa AbstractArray || Δt isa Tuple) ? Δt : (zero(Δt), Δt)

    return ODEProblem(f, u, tspan, p; kwargs...)
end

"""
Numerically integrate the orbit forward (or backward) in time, and return a new 
`AstrodynamicalOrbit` instance with identical parameters to the provided orbit.
"""
function propagate(
    orbit::AstrodynamicalModels.AstrodynamicalOrbit,
    Δt;
    stm = false,
    algorithm = Vern7(),
    reltol = 1e-12,
    abstol = 1e-12,
    kwargs...,
)
    problem = ODEProblem(orbit, Δt, stm = stm)
    return solve(problem, algorithm; reltol = reltol, abstol = abstol, kwargs...)
end

"""
Numerically integrate the orbit forward (or backward) in time, modifying the 
state vector in-place within the `AstrodynamicalOrbit` instance.
"""
function propagate!(
    orbit::AstrodynamicalModels.AstrodynamicalOrbit,
    Δt;
    stm = false,
    algorithm = Vern7(),
    reltol = 1e-12,
    abstol = 1e-12,
    kwargs...,
)
    overrides = (; save_everystep = false, save_start = false, save_end = true)
    options = (; kwargs..., overrides...)
    solution = propagate(
        orbit,
        Δt;
        stm = stm,
        algorithm = algorithm,
        reltol = reltol,
        abstol = abstol,
        options...,
    )
    AstrodynamicalModels.state(orbit) .= solution.u[end]
    return nothing
end

"""
Compute the monodromy matrix for any periodic orbit.
"""
function monodromy(
    orbit::AstrodynamicalModels.AstrodynamicalOrbit,
    Δt;
    algorithm = Vern7(),
    reltol = 1e-12,
    abstol = 1e-12,
    kwargs...,
)
    overrides = (; save_everystep = false, save_start = false, save_end = true)
    options = (; kwargs..., overrides...)
    solution = propagate(
        orbit,
        Δt;
        stm = true,
        algorithm = algorithm,
        reltol = reltol,
        abstol = abstol,
        options...,
    )

    N = length(AstrodynamicalModels.state(orbit))

    return AstrodynamicalModels.CartesianSTM(solution.u[end][begin+N:end])
end


"""
Solve for the monodromy matrix of the periodic orbit.
"""
function monodromy(
    u::AbstractVector,
    μ,
    T,
    f::Function;
    algorithm = Vern9(),
    reltol = 1e-12,
    abstol = 1e-12,
    save_everystep = false,
    kwargs...,
)
    problem = ODEProblem(
        f,
        [
            u[begin],
            u[begin+1],
            u[begin+2],
            u[begin+3],
            u[begin+4],
            u[begin+5],
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
        ],
        (zero(T), T),
        (μ,),
    )
    solution = solve(
        problem,
        algorithm;
        reltol = reltol,
        abstol = abstol,
        save_everystep = save_everystep,
        kwargs...,
    )

    if solution.u[begin][begin:begin+5] ≉ solution.u[end][begin:begin+5]
        @warn "The orbit does not appear to be periodic!"
    end

    return reshape((solution.u[end][begin+6:end]), 6, 6)
end

"""
Return a vector of orbits along the manifold which diverges from the provided 
halo orbit.
"""
function divergent_manifold(u, μ, Δt; eps = 1e-8, trajectories = nothing, kwargs...)
    orbit = Orbit(u, CR3BParameters(μ))
    Φ = monodromy(orbit, Δt)

    if isnothing(trajectories)
        ics = propagate(orbit, Δt; stm = true, kwargs...).u
    else
        if !isnothing(get(kwargs, :saveat, nothing))
            error(
                "you cannot provide both the `saveat` keyword argument and the `trajectories` keyword argument",
            )
        else
            ics =
                propagate(orbit, Δt; stm = true, saveat = Δt / trajectories, kwargs...).u[begin:begin+trajectories-1]
        end
    end

    states = [CartesianState(u[begin:begin+5]) for u in ics]
    stms = [CartesianSTM(u[begin+6:end]) for u in ics]

    for (u, ϕ) in zip(states, stms)
        diverge!(u, ϕ, Φ; eps = eps)
    end

    return map(
        u -> Orbit(CartesianState(u), AstrodynamicalModels.parameters(orbit)),
        states,
    )

end

"""
Return a vector of orbits along the manifold which converges to the provided 
halo orbit.
"""
function convergent_manifold(u, μ, Δt; eps = 1e-8, trajectories = nothing, kwargs...)
    orbit = Orbit(u, CR3BParameters(μ))
    Φ = monodromy(orbit, Δt)

    if isnothing(trajectories)
        ics = propagate(orbit, Δt; stm = true, kwargs...).u
    else
        if !isnothing(get(kwargs, :saveat, nothing))
            error(
                "you cannot provide both the `saveat` keyword argument and the `trajectories` keyword argument",
            )
        else
            ics =
                propagate(orbit, Δt; stm = true, saveat = Δt / trajectories, kwargs...).u[begin:begin+trajectories-1]
        end
    end

    states = [CartesianState(u[begin:begin+5]) for u in ics]
    stms = [CartesianSTM(u[begin+6:end]) for u in ics]

    for (u, ϕ) in zip(states, stms)
        converge!(u, ϕ, Φ; eps = eps)
    end

    return map(
        u -> Orbit(CartesianState(u), AstrodynamicalModels.parameters(orbit)),
        states,
    )

end

end