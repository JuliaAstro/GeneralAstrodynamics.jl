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
    ∂u.r =  u.v
    ∂u.v = -p.μ .* (u.r ./ norm(u.r,2)^3)
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

function SciMLBase.ODEProblem(orbit::CartesianOrbit, Δt::Real = ustrip(timeunit(orbit.state), period(orbit))) 
    u = CartesianState(orbit.r, orbit.v)
    t = (orbit.epoch, Δt)
    p = (μ = ustrip(lengthunit(orbit.state)^3 / timeunit(orbit.state)^2, mass_parameter(orbit.system)),)
    return ODEProblem(RestrictedTwoBodyTic!, u, t, p)
end
SciMLBase.ODEProblem(orbit::CartesianOrbit, Δt::Unitful.Time) = ODEProblem(orbit, ustrip(timeunit(orbit.state), Δt))

"""
Uses OrdinaryDiffEq solvers to propagate `orbit` Δt into the future.
All keyword arguments are passed directly to OrdinaryDiffEq solvers.

References:
* [1] https://diffeq.sciml.ai/v4.0/tutorials/ode_example.html
* [2] https://github.com/SciML/DifferentialEquations.jl/issues/393#issuecomment-658210231
* [3] https://discourse.julialang.org/t/smart-kwargs-dispatch/14571/15
"""
function propagate(orbit::CartesianOrbit, 
                   Δt = period(orbit);
                   thrust::Unitful.Acceleration = 0.0u"N/kg",
                   kwargs...)

    # Set default kwargs
    defaults = (;  reltol=1e-14, abstol=1e-14)
    options = merge(defaults, kwargs)

    # Initial conditions
    problem = ODEProblem(orbit, Δt)

    # Integrate! 
    sols = solve(problem; options...)

    # Return PropagationResult structure
    if sols.retcode != :success
        @warn "`DifferentialEquations` solvers returned code $(string(sols.retcode))."
    end
    return CartesianOrbit.(sols.t, CartesianState.(map(u->u.r, sols.u), map(u->u.v, sols.u); lengthunit = lengthunit(orbit), timeunit = timeunit(orbit)), (orbit.system,))

end
