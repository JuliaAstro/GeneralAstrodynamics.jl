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

@template (
    FUNCTIONS,
    METHODS,
    MACROS,
) = """
    $(SIGNATURES)

    $(DOCSTRING)
    """

@template (TYPES, CONSTANTS) = """
                               $(TYPEDEF)

                               $(DOCSTRING)
                               """

using AstrodynamicalModels
using ModelingToolkit, OrdinaryDiffEq

export propagate, propagate!

function OrdinaryDiffEq.ODEProblem(orbit::AstrodynamicalModels.AstrodynamicalOrbit, Δt; kwargs...)
    f = dynamics(orbit)
    u = AstrodynamicalModels.state(orbit)
    p = AstrodynamicalModels.parameters(orbit)
    tspan = (Δt isa AbstractArray || Δt isa Tuple) ? Δt : (zero(Δt), Δt)

    return ODEProblem(f, u, tspan, p; kwargs...)
end

"""
Numerically integrate the orbit forward (or backward) in time, and return a new 
`AstrodynamicalOrbit` instance with identical parameters to the provided orbit.
"""
function propagate(orbit::AstrodynamicalModels.AstrodynamicalOrbit, Δt; algorithm=Vern7(), reltol=1e-12, abstol=1e-12, kwargs...)
    problem = ODEProblem(orbit, Δt)
    return solve(problem, algorithm; reltol=reltol, abstol=abstol, kwargs...)
end

"""
Numerically integrate the orbit forward (or backward) in time, modifying the 
state vector in-place within the `AstrodynamicalOrbit` instance.
"""
function propagate!(orbit::AstrodynamicalModels.AstrodynamicalOrbit, Δt; algorithm=Vern7(), reltol=1e-12, abstol=1e-12, kwargs...)
    overrides = (; save_everystep=false, save_start=false, save_end=true)
    options = (; kwargs..., overrides...)
    solution = propagate(orbit, Δt; algorithm=algorithm, reltol=reltol, abstol=abstol, options...)
    AstrodynamicalModels.state(orbit) .= solution.u[end]
    return nothing
end



end