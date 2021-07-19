#
#   Kepler.jl
#
#   Solves Kepler's problem for TwoBody orbits.
#

"""
Solves Kepler's Problem for `orbit` and `Δtᵢ`.
"""
function kepler(orbit::RestrictedTwoBodyOrbit, Δtᵢ::Unitful.Time = period(orbit); tol=1e-6, max_iter=100) 

    conic_section = conic(orbit)

    # Guess χ₀
    if conic_section == Hyperbolic
        Δt = Δtᵢ
        χ₀ = sign(Δt) * √(-semimajor_axis(orbit)) * log(ℯ, (-2 * mass_parameter(orbit.system) / semimajor_axis(orbit) * Δt) / 
                (position_vector(orbit) ⋅ velocity_vector(orbit) + (sign(Δt) * √(-mass_parameter(orbit.system) * semimajor_axis(orbit)) * (1 - norm(position_vector(orbit)) / semimajor_axis(orbit)))))
    elseif conic_section == Parabolic
        Δt = Δtᵢ
        χ₀ = √(semi_parameter(orbit)) * tan(true_anomoly(orbit) / 2)
    else
        Δt = mod(Δtᵢ, period(orbit))
        χ₀ = √(mass_parameter(orbit.system)) * Δt / semimajor_axis(orbit)
    end

    # Iteratively solve for χ
    # TODO: Compare loop vs. recursion performance here.
    # There shouldn't be too large of a difference, since this tends
    # to converge with only a few iterations.
    χₙ, r, ψ, C₂, C₃ = χₖ(χ₀, Δt, position_vector(orbit), velocity_vector(orbit), semimajor_axis(orbit), mass_parameter(orbit.system), tol=tol, max_iter=max_iter)

    # Convert to a Orbit
    f = 1 - χₙ^2 / norm(position_vector(orbit)) * C₂
    ḟ = √(mass_parameter(orbit.system)) / (norm(position_vector(orbit)) * r) * χₙ * (ψ * C₃ - 1)
    g = Δt - (χₙ^3 / √(mass_parameter(orbit.system))) * C₃
    ġ = 1 - (χₙ^2 / r) * C₂

    return CartesianOrbit(f * position_vector(orbit) + g * velocity_vector(orbit), ḟ * position_vector(orbit) + ġ * velocity_vector(orbit), orbit.system, epoch(orbit.state))
end

"""
Implements Kepler's Prediction Problem for generic `r`, `v`.

Arguments:
* `r`: Initial spacecraft position. Any abstract vector with length 3 and type `<: Unitful.Length`
* `v`: Initial spacecraft velocity. Any abstract vector with length 3 and type `<: Unitful.Velocity`
* `μ`: Central body mass parameter. Any abstract vector with length 3 and type `<: Unitful.AbstractQuantity` with units compatabile with `km^3/s^2`
* `Δtᵢ`: Propagation time. Any scalar with type `<: Unitful.Time`
"""
function kepler(r, v, μ, Δtᵢ; tol=1e-6, max_iter=100)
    initial = RestrictedTwoBodyOrbit(r, v, μ)
    final   = kepler(final, Δtᵢ; tol=tol, max_iter=max_iter)
    return position_vector(r), velocity_vector(v)
end

kepler(orbit::RestrictedTwoBodyOrbit, Δtᵢ::Real; kwargs...) = kepler(orbit, timeunit(orbit) * Δtᵢ; kwargs...)

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
        @error "Failed to converge!"
        return χₙ, r, ψ, C₂, C₃
    elseif abs(χₙ₊₁ - χₙ) < oneunit(χₙ) * tol
        return  χₙ, r, ψ, C₂, C₃
    else
        return χₖ(χₙ₊₁, Δt, rᵢ₀, vᵢ₀, a, μ; iter=iter+1, tol=tol, max_iter=max_iter)
    end
    
end