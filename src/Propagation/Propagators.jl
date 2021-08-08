#
# Restricted Two-body problem propagation
#

"""
Create's an ODEProblem for a `R2BP` orbit.
"""
function SciMLBase.ODEProblem(orbit::R2BPOrbit, Δt::Number)
    u  = state(orbit)
    Δt = eltype(state(orbit))(Δt)
    ts = Δt isa Unitful.Length ? (zero(ustrip(timeunit(orbit), Δt)), ustrip(timeunit(orbit), Δt)) : (zero(Δt), Δt)
    ts = (min(ts...), max(ts...))
    p  = system(orbit)
    return ODEProblem(R2BPVectorField, u, ts, p)
end

"""
Create's an ODEProblem for a `R2BP` orbit.
"""
function SciMLBase.ODEProblem(orbit::R2BPOrbit, Δt::Real)
    u  = state(orbit)
    Δt = eltype(state(orbit))(Δt)
    ts = (zero(Δt), Δt)
    ts = (min(ts...), max(ts...))
    p  = system(orbit)
    return ODEProblem(CR3BPVectorField, u, ts, p)
end

"""
Propagates an orbit forward or backward in time.
Use `algorithm` to set the desired numerical integration
algorithm, e.g. `algorithm = Tsit5()`. All other `kwargs` 
are passed directly to `DifferentialEquations.solve`.
"""
function propagate(orbit::Orbit, Δt::Number; algorithm=nothing, kwargs...)
    problem = ODEProblem(orbit, Δt)

    if isnothing(algorithm)
        solution = solve(problem; kwargs...)
    else 
        solution =  solve(problem, algorithm; kwargs...)
    end

    return Trajectory{frame(orbit), typeof(epoch(orbit)), typeof(solution)}(epoch(orbit), solution)
end