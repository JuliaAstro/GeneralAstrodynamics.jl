#
# Restricted Two-body Problem Propagators
#

"""
Dynamics for ideal two-body numerical integration.
Compatable with `DifferentialEquations` solvers!
"""
function R2BPTic!(∂u, u, p, t)
    ∂u[1:3] .= u[4:6]
    ∂u[4:6] .= -p.μ .* (u[1:3] ./ norm(u[1:3])^3)
    return nothing
end

"""
Given a `CartesianOrbit`, return a `DifferentialEquations.ODEProblem` instance
that describes the Restricted Two-body Problem.
"""
function SciMLBase.ODEProblem(orbit::CartesianOrbit, Δt::Real = ustrip(timeunit(orbit.state), period(orbit))) 
    u = orbit.state.rv
    t = (orbit.state.t, Δt)
    p = (μ = ustrip(lengthunit(orbit.state)^3 / timeunit(orbit.state)^2, mass_parameter(orbit.system)),)
    return ODEProblem(R2BPTic!, u, t, p)
end

"""
Given a `CartesianOrbit`, return a `DifferentialEquations.ODEProblem` instance
that describes the Restricted Two-body Problem.
"""
SciMLBase.ODEProblem(orbit::CartesianOrbit, Δt::Unitful.Time) = ODEProblem(orbit, ustrip(timeunit(orbit.state), Δt))

"""
Uses `DifferentialEquations` solvers to propagate `orbit` Δt into the future.
All keyword arguments are passed directly to `DifferentialEquations.solve`.
"""
function propagate(orbit::CartesianOrbit, 
                   Δt = period(orbit);
                   kwargs...)

    # Set default kwargs
    defaults = (;  reltol=1e-14, abstol=1e-14)
    options = merge(defaults, kwargs)

    # Initial conditions
    problem = ODEProblem(orbit, Δt)

    # Integrate! 
    sols = solve(problem; options...)

    # Return PropagationResult structure
    if sols.retcode != :Success
        @warn "DifferentialEquations solvers returned code $(string(sols.retcode))."
    end

    return [
        CartesianOrbit(
            sols.u[i][1:3] * lengthunit(orbit.state), 
            sols.u[i][4:6] * velocityunit(orbit.state), 
            orbit.system, sols.t[i] * timeunit(orbit.state))
        for i ∈ 1:length(sols.t)
    ]
end

propagate(::KeplerianOrbit, args...; kwargs...) = throw(
    ArgumentError(
        "Keplerian orbit propagation isn't supported at this time! Please convert your orbit to a `CartesianOrbit` for propagation."
    )
)