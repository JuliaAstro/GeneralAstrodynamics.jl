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
struct TwobodyPropagationResult{
        T<:Unitful.Time,
        O<:TwoBodySystem,
        VT<:AbstractVector{T},
        VO<:AbstractVector{O}
} <: PropagationResult 

    t::VT
    step::VO
    ode_solution::ODESolution

end

"""
    twobody_tic

Currently not exported. Used for two-body numerical integration.
"""
function twobody_tic!(∂u, u, p, t)

    ∂u.rᵢ =  u.vᵢ
    ∂u.vᵢ = (-p.μ .* (u.rᵢ ./ norm(u.rᵢ,2)^3)) .+ (normalize(u.vᵢ) .* p.T)

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
                   ode_alg::OrdinaryDiffEqAlgorithm=Tsit5(),
                   thrust::Unitful.Acceleration=0.0u"N/kg";
                   kwargs...)

    # Referencing:
    # [1] https://diffeq.sciml.ai/v4.0/tutorials/ode_example.html
    # [2] https://github.com/SciML/DifferentialEquations.jl/issues/393#issuecomment-658210231
    # [3] https://discourse.julialang.org/t/smart-kwargs-dispatch/14571/15

    # Set default kwargs (modified from [3])
    defaults = (;  reltol=1e-14, abstol=1e-14)
    options = merge(defaults, kwargs)

    # Initial conditions
    r₀ = Array(ustrip.(u"m",   orbit.rᵢ))
    v₀ = Array(ustrip.(u"m/s", orbit.vᵢ))

    # Define the problem (modified from [2])
    problem = ODEProblem(twobody_tic!, 
                         ComponentArray((rᵢ=r₀, vᵢ=v₀)), 
                         ustrip.(u"s", (0.0u"s", Δt)), 
                         ComponentArray((μ=ustrip(u"m^3 / s^2", orbit.body.μ), 
                                         T=ustrip(u"N/kg", thrust))))

    # Solve the problem! 
    sols = solve(problem, ode_alg; options...)

    # Return PropagationResult structure
    return TwobodyPropagationResult(
        u"s" .* sols.t,
        map(x -> Orbit(u"m" * x.rᵢ, u"m/s" * x.vᵢ, orbit.body), sols.u),
        sols
    )

end