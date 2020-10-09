#
#   propagator.jl
#
#   Includes functions and structures for propagating orbits 
#   within the two-body problem.
#

"""
    TwobodyPropagationResult <: PropagationResult

Struct to hold two-body propagation results.
"""
struct TwobodyPropagationResult <: PropagationResult 

    t::AbstractVector{Unitful.Time{Float64}}
    step
    ode_solution::ODESolution

end

"""
    twobody_tic

Currently not exported. Used for two-body numerical integration.
"""
function twobody_tic(du, u, p, t)

    # Note the citation [2] above - this function was copied and
    # modified from [2]. I originally had a 6 state function,
    # but I had trouble with mixing units. The solution shown
    # in [2] allows the use of mixed units through the ComponentArrays
    # package.

    du.r̅ =  u.v̅
    du.v̅ = -p.μ * (u.r̅ ./ norm(u.r̅,2)^3)

end

"""
    propagate_twobody(orbit::Orbit, 
                   Δt::Unitful.Time=orbital_period(orbit), 
                   ode_alg::OrdinaryDiffEqAlgorithm=Tsit5(); 
                   kwargs...)

Uses OrdinaryDiffEq solvers to propagate `orbit` Δt into the future.
All keyword arguments are passed directly to OrdinaryDiffEq solvers.
"""
function propagate(orbit::Orbit, 
                   Δt::Unitful.Time=orbital_period(orbit), 
                   ode_alg::OrdinaryDiffEqAlgorithm=Tsit5(); 
                   kwargs...)

    # Referencing:
    # [1] https://diffeq.sciml.ai/v4.0/tutorials/ode_example.html
    # [2] https://github.com/SciML/DifferentialEquations.jl/issues/393#issuecomment-658210231
    # [3] https://discourse.julialang.org/t/smart-kwargs-dispatch/14571/15

    # Set default kwargs (modified from [3])
    defaults = (;  reltol=1e-14, abstol=1e-14)
    options = merge(defaults, kwargs)

    # Initial conditions
    r₀ = Array(ustrip.(u"m",orbit.r̅))
    v₀ = Array(ustrip.(u"m/s", orbit.v̅))

    # Define the problem (modified from [2])
    problem = ODEProblem(
    twobody_tic, 
    ComponentArray((r̅=r₀, v̅=v₀)), 
    ustrip.(u"s", (0.0u"s", Δt)), 
    ComponentArray((μ=ustrip(u"m^3 / s^2", orbit.body.μ))))

    # Solve the problem! 
    sols = solve(problem, ode_alg; options...)

    # Return PropagationResult structure
    return TwobodyPropagationResult(
        u"s" * sols.t,
        map(x -> Orbit(u"m" * x.r̅, u"m/s" * x.v̅, orbit.body), sols.u),
        sols
    )

end