#
#   propagator.jl
#
#   Includes functions and structures for propagating orbits 
#   within the two-body problem.
#

"""
Currently not exported. Used for ideal two-body numerical integration.
"""
function RestrictedTwoBodyTic!(∂u, u, p, t)
    ∂u.rᵢ =  u.vᵢ
    ∂u.vᵢ = -p.μ .* (u.rᵢ ./ norm(u.rᵢ,2)^3)
    return nothing
end

"""
Currently not exported. Used for ideal two-body numerical integration.
"""
function RestrictedBiasedTwoBodyTic!(∂u, u, p, t)
    ∂u.rᵢ =  u.vᵢ
    ∂u.vᵢ = (-p.μ .* (u.rᵢ ./ norm(u.rᵢ,2)^3)) .+ (normalize(u.vᵢ) .* p.T)
    return nothing
end

"""
Uses OrdinaryDiffEq solvers to propagate `orbit` Δt into the future.
All keyword arguments are passed directly to OrdinaryDiffEq solvers.

References:
* [1] https://diffeq.sciml.ai/v4.0/tutorials/ode_example.html
* [2] https://github.com/SciML/DifferentialEquations.jl/issues/393#issuecomment-658210231
* [3] https://discourse.julialang.org/t/smart-kwargs-dispatch/14571/15
"""
function propagate(orbit::CartesianOrbit{C, T}, 
                   Δt::Unitful.Time = period(orbit);
                   thrust::Unitful.Acceleration = 0.0u"N/kg",
                   kwargs...) where C where T

    # Set default kwargs
    defaults = (;  reltol=1e-14, abstol=1e-14)
    options = merge(defaults, kwargs)

    # Initial conditions
    r₀ = Array(ustrip.(u"m",   position_vector(orbit)))
    v₀ = Array(ustrip.(u"m/s", velocity_vector(orbit)))
    u₀ = ComponentArray((rᵢ=r₀, vᵢ=v₀))
    ts = T.(ustrip.(u"s", (zero(Δt), Δt)))
    f  = thrust == zero(thrust) ? RestrictedTwoBodyTic! : RestrictedBiasedTwoBodyTic!
    p  = let
        if thrust == zero(thrust)
            ComponentArray((μ=ustrip(u"m^3/s^2", orbit.body.μ)))
        else
            ComponentArray((μ=ustrip(u"m^3/s^2", orbit.body.μ), T=ustrip(u"N/kg", thrust)))
        end
    end

    # Integrate! 
    sols = solve(ODEProblem(f, u₀, ts, p); options...)

    # Return PropagationResult structure
    return Trajectory(
            map(x -> CartesianOrbit(u"m" * x.rᵢ, u"m/s" * x.vᵢ, orbit.body), sols.u),
            sols.t .* u"s",
            sols.retcode
    )

end

"""
Uses OrdinaryDiffEq solvers to propagate `orbit` Δt into the future.
All keyword arguments are passed directly to OrdinaryDiffEq solvers.
"""
function propagate(orbit::KeplerianOrbit, 
                   Δt::Unitful.Time=period(orbit), 
                   ode_alg::OrdinaryDiffEqAlgorithm=Tsit5(),
                   thrust::Unitful.Acceleration=0.0u"N/kg";
                   kwargs...)
    traj = propagate(CartesianOrbit(orbit), Δt, ode_alg, thrust, kwargs...)
    return Trajectory(KeplerianOrbit.(traj.step), traj.t, traj.status)
end