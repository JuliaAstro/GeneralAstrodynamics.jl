#
# Restricted Two-body problem propagation
#

"""
Create's an ODEProblem for a `R2BP` orbit.
"""
function SciMLBase.ODEProblem(orbit::R2BPOrbit, Δt::Number)
    u  = state(orbit)
    Δt = eltype(state(orbit))(Δt)
    ts = Δt isa Unitful.Time ? (zero(ustrip(timeunit(orbit), Δt)), ustrip(timeunit(orbit), Δt)) : (zero(Δt), Δt)
    p  = system(orbit)
    return ODEProblem(R2BPFunction(), u, ts, p)
end

"""
Create's an ODEProblem for a `R2BP` orbit.
"""
function SciMLBase.ODEProblem(orbit::CR3BPOrbit, Δt::Real)
    u  = state(orbit)
    Δt = eltype(state(orbit))(Δt)
    ts = (zero(Δt), Δt)
    p  = system(orbit)
    f  = u isa CartesianStateWithSTM ? CR3BPFunction(; stm=true) : CR3BPFunction()
    return ODEProblem(f, u, ts, p)
end

"""
Propagates an orbit forward or backward in time.
Use `algorithm` to set the desired numerical integration
algorithm, e.g. `algorithm = Tsit5()`. All other `kwargs` 
are passed directly to `DifferentialEquations.solve`.
"""
function propagate(orbit::Orbit, Δt::Number; algorithm=nothing, kwargs...)
    defaults = (; reltol=1e-14, abstol=1e-14)
    options  = merge(defaults, kwargs)
    problem  = ODEProblem(orbit, Δt)

    if isnothing(algorithm)
        solution = solve(problem; options...)
    else 
        solution =  solve(problem, algorithm; options...)
    end

    return Trajectory{frame(orbit), typeof(solution.prob.u0), typeof(solution.prob.p), typeof(epoch(orbit)), typeof(solution)}(epoch(orbit), solution)
end