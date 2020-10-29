#
# Solver for Lambert's problem
#
# References:
# [1] David, A. "Vallado. Fundamentals of Astrodynamics and Applications." (2013).
#

"""
    lambert(r̅₁, r̅₂, μ, Δt, trajectory=:short; tol=1e-6, max_iter=100)

Solves Lambert's problem through the use of univeral variables.
"""
function lambert(r̅₁, r̅₂, μ, Δt, trajectory=:short; tol=1e-6, max_iter=100)

    # Specify short way, or long way trajectory
    if trajectory == :short
        tₘ = 1
    elseif trajectory == :long
        tₘ = -1
    else
        throw(ArgumentError("`trajectory` must be set to `:short` or `:long`"))
    end

    r₁ = norm(r̅₁)
    r₂ = norm(r̅₂)

    cosΔν = (r̅₁⋅r̅₂) / (r₁*r₂)
    Δν    = asin(u"rad", tₘ * √(1 - (cosΔν)^2))

    A = upreferred(tₘ * √(r₂*r₁ * (1 + cosΔν)))

    if A ≈ 0
        throw(ErrorException("Can't calculate the orbit."))
    end

    ψₙ = 0.0
    c₂ = 1/2
    c₃ = 1/6

    ψ₊ = 4π^2
    ψ₋ = -4π
    yₙ = r₁ + r₂ + (A * (ψₙ*c₃ - 1) / √(c₂))

    Δtₙ = Δt + 1u"s"
    iter = 0

    while (iter < max_iter) && 
          (abs(Δtₙ - Δt) > (tol * oneunit(Δt))) || (A > 0u"m" && yₙ < 0u"m")

        yₙ = r₁ + r₂ + (A * (ψₙ*c₃ - 1) / √(c₂))
        χₙ = √(yₙ / c₂)
        Δtₙ = (χₙ^3 * c₃ + A*√(yₙ)) / √(μ)

        if Δtₙ < Δt
            ψ₋ = ψₙ
        else
            ψ₊ = ψₙ
        end

        ψₙ = (ψ₊ + ψ₋) / 2
        if ψₙ > tol
            C₂ = (1 - cos(√(ψₙ))) /  ψₙ
            C₃ = (√(ψₙ) - sin(√(ψₙ))) / √(ψₙ^3)
        elseif ψₙ < -tol
            C₂ = (1 - cosh(√(-ψₙ))) / ψₙ
            C₃ = (sinh(√(-ψₙ)) - √(-ψₙ)) / √((-ψₙ)^3)
        else
            C₂ = 1.0 / 2.0
            C₃ = 1.0 / 6.0
        end

        iter += 1

    end

    f = 1 - yₙ/r₁
    ġ = 1 - yₙ/r₂
    g = A * √(yₙ/μ)

    v̅₁ = (r̅₂ .- (f .* r̅₁)) ./ g
    v̅₂ = ((ġ .* r̅₂) .- r̅₁) ./ g

    return v̅₁, v̅₂
    
end