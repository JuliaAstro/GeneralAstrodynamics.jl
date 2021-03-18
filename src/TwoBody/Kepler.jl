#
#   Kepler.jl
#
#   Solves Kepler's problem for TwoBody orbits.
#

"""
Solves Kepler's Problem for `orbit` and `Δtᵢ`.
"""
function kepler(orbit::O, Δtᵢ::T = period(orbit); tol=1e-6, max_iter=100) where O <: RestrictedTwoBodySystem where T<:Unitful.Time

    conic_section = conic(orbit)

    # Guess χ₀
    if conic_section == Circular || conic_section == Elliptical

        Δt = mod(Δtᵢ, period(orbit))
        χ₀ = √(orbit.body.μ) * Δt / semimajor_axis(orbit)

    elseif conic_section == Hyperbolic
        
        Δt = Δtᵢ
        χ₀ = sign(Δt) * √(-semimajor_axis(orbit)) * log(ℯ, (-2 * orbit.body.μ / semimajor_axis(orbit) * Δt) / 
                (radius_vector(orbit) ⋅ velocity_vector(orbit) + (sign(Δt) * √(-orbit.body.μ * semimajor_axis(orbit)) * (1 - norm(radius_vector(orbit)) / semimajor_axis(orbit)))))

    elseif conic_section == Parabolic

        Δt = Δtᵢ
        χ₀ = √(semi_parameter(orbit)) * tan(true_anomoly(orbit) / 2)

    else

        @warn "Kepler's problem failed to converge."
        return Orbit([NaN, NaN, NaN] * u"km", [NaN, NaN, NaN] * u"km/s", orbit.body)

    end

    # Iteratively solve for χ
    # TODO: Compare loop vs. recursion performance here.
    # There shouldn't be too large of a difference, since this tends
    # to converge with only a few iterations.
    χₙ, r, ψ, C₂, C₃ = χₖ(χ₀, Δt, radius_vector(orbit), velocity_vector(orbit), semimajor_axis(orbit), orbit.body.μ, tol=tol, max_iter=max_iter)

    # Convert to a Orbit
    f = 1 - χₙ^2 / norm(radius_vector(orbit)) * C₂
    ḟ = √(orbit.body.μ) / (norm(radius_vector(orbit)) * r) * χₙ * (ψ * C₃ - 1)
    g = Δt - (χₙ^3 / √(orbit.body.μ)) * C₃
    ġ = 1 - (χₙ^2 / r) * C₂

    return Orbit(f * radius_vector(orbit) + g * velocity_vector(orbit), ḟ * radius_vector(orbit) + ġ * velocity_vector(orbit), orbit.body)

end

function χₖ(χₙ, Δt, rᵢ₀, vᵢ₀, a, μ; iter=1, tol=1e-14, max_iter=100)
    
    r₀ = norm(rᵢ₀)
    ψ = upreferred(χₙ^2 / a)

    if ψ > tol
        C₂ = (1 - cos(√(ψ))) /  ψ
        C₃ = (√(ψ) - sin(√(ψ))) / √(ψ^3)
    elseif ψ < -tol
        println(√(-ψ))
        C₂ = (1 - cosh(√(-ψ))) / ψ
        C₃ = (sinh(√(-ψ)) - √(-ψ)) / √((-ψ)^3)
    else
        C₂ = 1.0 / 2.0
        C₃ = 1.0 / 6.0
    end

    r = χₙ^2 * C₂ + (rᵢ₀ ⋅ vᵢ₀) * χₙ / √(μ) * (1 - ψ*C₃) + r₀ * (1 - ψ * C₂)
    χₙ₊₁ = χₙ + ((√(μ) * Δt - χₙ^3 * C₃ - (rᵢ₀ ⋅ vᵢ₀) / √(μ) * χₙ^2 * C₂ - r₀ * χₙ * (1 - ψ * C₃)) / r)

    if iter > max_iter
        return NaN, NaN, NaN, NaN, NaN
    elseif abs(χₙ₊₁ - χₙ) < oneunit(χₙ) * tol
        return  χₙ, r, ψ, C₂, C₃
    else
        return χₖ(χₙ₊₁, Δt, rᵢ₀, vᵢ₀, a, μ; iter=iter+1, tol=tol, max_iter=max_iter)
    end
    
end