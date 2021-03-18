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
    accel!(∂u.vₛ, u.rₛ, u.vₛ, p.μ)
    return nothing
end

"""
Uses OrdinaryDiffEq solvers to propagate `orbit` Δt into the future.
All keyword arguments are passed directly to OrdinaryDiffEq solvers.
"""
function propagate(sys::NondimensionalThreeBodyState, Δt::T = sys.Δt; kwargs...) where T<:Real

    # Referencing:
    # [1] https://diffeq.sciml.ai/v4.0/tutorials/ode_example.html
    # [2] https://github.com/SciML/DifferentialEquations.jl/issues/393#issuecomment-658210231
    # [3] https://discourse.julialang.org/t/smart-kwargs-dispatch/14571/15

    # Set default kwargs (modified from [3])
    defaults = (;  reltol=1e-14, abstol=1e-14)
    options = merge(defaults, kwargs)

    # Initial conditions
    u₀ = ComponentArray((rₛ=sys.r, vₛ=sys.v))
    ts = (zero(Δt), Δt)
    p  = ComponentArray((μ=sys.μ, x₁=-sys.μ, x₂=1-sys.μ))

    # Numerically integrate!
    sols = solve(ODEProblem(RestrictedThreeBodyTic!, u₀, ts, p), ode_alg; options...)

    # Return PropagationResult structure
    return Trajectory(
        map(step->NondimensionalThreeBodyState(step.rₛ, step.vₛ, sys.μ, sys.Δt, sys.DU, sys.DT), sols.u),
        sols.t,
        sols.retcode
    )

end

