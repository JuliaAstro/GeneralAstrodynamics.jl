#
#   propagator.jl
#
#   Includes functions and structures for propagating orbits 
#   within the two-body problem.
#

"""
    TwobodyPropagationResult <: PropagationResult

Wrapper for ODESolution, with optional units.
"""
struct TwobodyPropagationResult <: PropagationResult

    t::AbstractVector{Unitful.Time{Float64}}
    r̅::AbstractMatrix{Unitful.Length{Float64}}
    v̅::AbstractMatrix{Unitful.Velocity{Float64}}
    ode_solution::ODESolution

end

"""
    propagate_twobody(orbit::TwoBodyOrbit, 
                   Δt::Unitful.Time=orbital_period(orbit), 
                   ode_alg::OrdinaryDiffEqAlgorithm=Tsit5(); 
                   kwargs...)

Uses OrdinaryDiffEq solvers to propagate `orbit` Δt into the future.
All keyword arguments are passed directly to OrdinaryDiffEq solvers.
"""
function propagate_twobody(orbit::TwoBodyOrbit, 
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
    r₀ = Array(ustrip.(u"km",orbit.r̅))
    v₀ = Array(ustrip.(u"km/s", orbit.v̅))

    # Define the problem (modified from [2])
    problem = ODEProblem(
    twobody_tic, 
    ComponentArray((r̅=r₀, v̅=v₀)), 
    ustrip.(u"s", (0.0u"s", Δt)), 
    ComponentArray((μ=ustrip(u"km^3 / s^2", orbit.body.μ))))

    # Solve the problem! 
    sols = solve(problem, ode_alg; options...)

    # Return PropagationResult structure
    return TwobodyPropagationResult(
        u"s" * sols.t,
        u"km" * vcat(map(x->x.r̅', sols.u)...),
        u"km/s" * vcat(map(x->x.v̅', sols.u)...),
        sols
    )

end


function twobody_tic(du, u, p, t)

    # Note the citation [2] above - this function was copied and
    # modified from [2]. I originally had a 6 state function,
    # but I had trouble with mixing units. The solution shown
    # in [2] allows the use of mixed units through the ComponentArrays
    # package.

    du.r̅ =  u.v̅
    du.v̅ = -p.μ * (u.r̅ ./ norm(u.r̅,2)^3)

end