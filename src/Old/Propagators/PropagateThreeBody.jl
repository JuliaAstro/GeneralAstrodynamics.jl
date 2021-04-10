#
#   PropagateThreeBody.jl
#
#   Includes functions and structures for propagating orbits 
#   within the circular restricted three-body problem.
#

"""
    threebody_tic

Currently not exported. Used for two-body numerical integration.
"""
function RestrictedThreeBodyTic!(∂u, u, p, t=0)
    ∂u.rₛ =  u.vₛ
    ThreeBody.accel!(∂u.vₛ, u.rₛ, u.vₛ, p.μ)
    return nothing
end

"""
Uses OrdinaryDiffEq solvers to propagate `sys` Δt into the future.
All keyword arguments are passed directly to OrdinaryDiffEq solvers.
"""
function propagate(r, v, μ, Δt; kwargs...)

    # Referencing:
    # [1] https://diffeq.sciml.ai/v4.0/tutorials/ode_example.html
    # [2] https://github.com/SciML/DifferentialEquations.jl/issues/393#issuecomment-658210231
    # [3] https://discourse.julialang.org/t/smart-kwargs-dispatch/14571/15

    # Set default kwargs (modified from [3])
    defaults = (;  reltol=1e-14, abstol=1e-14)
    options = merge(defaults, kwargs)

    # Initial conditions
    u₀ = ComponentArray((rₛ=r, vₛ=v))
    ts = (zero(Δt), Δt)
    p  = ComponentArray((μ=μ, x₁=-μ, x₂=1-μ))

    # Numerically integrate!
    sols = solve(ODEProblem(RestrictedThreeBodyTic!, u₀, ts, p); options...)

    # Return PropagationResult structure
    return sols
end
