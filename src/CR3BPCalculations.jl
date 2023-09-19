"""
Common calculations within Circular Restricted Three Body Problem dynamics.

# Extended Help

## Imports

$(IMPORTS)

## Exports

$(EXPORTS)

"""
module CR3BPCalculations

using StaticArrays
using LinearAlgebra
using DocStringExtensions
using Roots

@template (
    FUNCTIONS,
    METHODS,
    MACROS,
) = """
    $(SIGNATURES)

    !!! warning "CR3BP Dynamics"
        This computation is valid for Circular Restricted Three Body Problem dynamics.

    $(DOCSTRING)
    """
export distance_scaling,
    time_scaling,
    reduced_mass,
    velocity_scaling,
    nondimensional,
    redimensioned,
    nondimensional_radii,
    distance_to_primary,
    distance_to_secondary,
    potential_energy,
    jacobi_constant,
    primary_synodic_position,
    secondary_synodic_position,
    lagrange_point,
    zero_velocity_curves,
    richardson_halo,
    perturbation,
    convergent_direction,
    divergent_direction

"""
The length scale factor used to nondimensionalize CR3BP states.
"""
Base.@pure distance_scaling(a) = a
Base.@pure distance_scaling(a, μ₁, μ₂) = a

"""
The time scale factor used to nondimensionalize CR3BP states.
"""
time_scaling(a, μ₁, μ₂) = 2π * √(a^3 / (μ₁ + μ₂))

"""
The velocity scale factor used to nondimensionalize CR3BP states.
"""
velocity_scaling(a, μ₁, μ₂) = distance_scaling(a) / time_scaling(a, μ₁, μ₂)

"""
Return the reduced mass.
"""
reduced_mass(μ₁, μ₂) = min(μ₁, μ₂) / (μ₁ + μ₂)

"""
Normalizes a CR3BP orbit in the rotating reference frame.
"""
function nondimensional(r::AbstractVector, v::AbstractVector, t, a, μ₁, μ₂)
    rₙ = r ./ a
    μₜ = μ₁ + μ₂
    Tₛ = 2π * √(a^3 / μₜ)
    vₙ = v ./ (a / Tₛ)
    tₙ = t / Tₛ
    μₙ = min(μ₁, μ₂) / μₜ

    return rₙ, vₙ, tₙ, μₙ, a, Tₛ
end

"""
Redimensionalizes a CR3BP orbit in the rotating reference frame.
"""
function redimensioned(rₙ::AbstractVector, vₙ::AbstractVector, tₙ, μ, a, Tₛ)
    r = rₙ .* a
    v = vₙ .* a ./ Tₛ
    t = tₙ * Tₛ

    sum_μs = a^3 / ((Tₛ / 2π)^2)
    μ₂ = μₙ * sum_μs
    μ₁ = sum_μ - μ₂

    return r, v, t, a, (μ₁, μ₂)
end

"""
Returns the spacecraft's nondimensional position w.r.t. body 1 (or 2).
"""
nondimensional_radii(r::AbstractVector, xᵢ) = sqrt((r[1] - xᵢ)^2 + r[2]^2 + r[3]^2)

"""
Returns synodic distance to primary body.
"""
distance_to_primary(r, μ) = r + μ

"""
Returns synodic distance to secondary body.
"""
distance_to_secondary(r, μ) = r - (one(μ) - μ)

"""
Returns the potential energy `U` in the Synodic frame with Normalized units.
"""
function potential_energy(r, μ)
    return (r[1]^2 + r[2]^2) / 2 +
           ((1 - μ) / nondimensional_radii(r, -μ)) +
           (μ / nondimensional_radii(r, 1 - μ))
end

"""
Returns the Jacobi Constant, `C` in the Synodic frame with Normalized units.
"""
jacobi_constant(r, v, μ) = 2 * potential_energy(r, μ) - (v ⋅ v)

"""
Given the Synodic frame vector, returns the vector in the barycentric-inertial reference frame.
"""
function synodic_to_inertial(rₛ::AbstractVector, t, ω)
    θ = ω * t
    O = zero(θ)
    l = one(θ)

    ᴵTₛ = SMatrix{3,3}(cos(θ), -sin(θ), O, sin(θ), cos(θ), O, O, O, l)

    return ᴵTₛ * rₛ
end

"""
Given a barycentric-inertial Cartesian state, returns the state in the synodic (rotating) reference frame.
"""
function inertial_to_synodic(rᵢ::AbstractVector, t, ω)
    θ = ω * t
    O = zero(θ)
    l = one(θ)

    ᴵTₛ = SMatrix{3,3}(cos(θ), -sin(θ), O, sin(θ), cos(θ), O, O, O, l)

    ˢTᵢ = inv(ᴵTₛ)

    return ˢTᵢ * rᵢ
end

"""
Position of primary body.
"""
function primary_synodic_position(μ)
    return SVector{3}(-μ, zero(μ), zero(μ))
end

"""
Position of primary body.
"""
function secondary_synodic_position(μ)
    return SVector{3}(one(μ) - μ, zero(μ), zero(μ))
end

"""
Returns the lagrange points for a CR3BP system.

__Arguments:__ 
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `L`: Langrange points requested, must be in range [1,5]

__Outputs:__
- Tuple of Lagrange points
- Throws `ArgumentError` if L is out of range [1,5]

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/)
"""
function lagrange_point(μ::Real, L=1:5)
    if !all(L[i] ∈ (1, 2, 3, 4, 5) for i = 1:length(L))
        throw(ArgumentError("Requested lagrange points must be in range [1,5]"))
    end

    expressions = SVector{3}(
        x -> x - (1 - μ) / (x + μ)^2 + μ / (x + μ - 1)^2,
        x -> x - (1 - μ) / (x + μ)^2 - μ / (x + μ - 1)^2,
        x -> x + (1 - μ) / (x + μ)^2 + μ / (x + μ + 1)^2,
    )

    return (
        map(f -> [find_zero(f, (-3, 3)), 0, 0], expressions)...,
        [(1 / 2) - μ, √(3) / 2, 0],
        [(1 / 2) - μ, -√(3) / 2, 0],
    )[L]
end

"""
Returns an analytical solution for a Halo orbit about `L`.

# Extended Help

## Arguments 
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `L`: Lagrange point to orbit (L1 or L2).
- `Z`: Desired non-dimensional Z-amplitude for Halo orbit.
- `hemisphere`: Specifies northern or southern Halo orbit.
- `ϕ`: Desired Halo orbit phase.
- `steps`: Number of non-dimensional timepoints in returned state.

## Outputs
- Synodic position vectors `r::Array{<:AbstractFloat}`
- Synodic velocity vectors `v::Array{<:Abstractfloat}`.
- Halo orbit period `Τ`.
- Throws `ArgumentError` if L is not `1` or `2`.

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/).
"""
function richardson_halo(μ, L::Int; Z=0.0, hemisphere=:northern, ϕ=0.0, steps=1)
    if L == 1
        point = first(lagrange_point(μ, 1))
        γ = abs(one(μ) - μ - point)
        n = collect(1:4) .* one(μ)
        c = @. (μ + (-one(1))^n * (one(μ) - μ)γ^(n + 1)) / (γ^3 * (one(μ) - γ^(n + 1)))
    elseif L == 2
        point = first(lagrange_point(μ, 2))
        γ = abs(point - one(μ) + μ)
        n = collect(1:4) .* one(μ)
        c = @. ((-one(μ))^n * μ + (-one(μ))^n * (one(μ) - μ)γ^(n + 1)) /
               (γ^3 * (one(μ) + γ^(n + 1)))
    else
        throw(ArgumentError("Only Halo orbits about L1 or L2 are supported."))
    end

    ωₚ = √((2 - c[2] + √((9c[2]^2 - 8c[2]))) / 2)
    k = (ωₚ^2 + 1 + 2c[2]) / (2ωₚ)

    d₁ = (3ωₚ^2 / k) * (k * (6ωₚ^2 - 1) - 2ωₚ)
    d₂ = (8ωₚ^2 / k) * (k * (11ωₚ^2 - 1) - 2ωₚ)
    a₂₁ = (3c[3] * (k^2 - 2)) / (4(1 + 2c[2]))
    a₂₂ = (3c[3]) / (4(1 + 2c[2]))
    a₂₃ = (-3c[3]ωₚ / (4k * d₁)) * (3k^3 * ωₚ - 6k * (k - ωₚ) + 4)
    a₂₄ = (-3c[3]ωₚ / (4k * d₁)) * (2 + 3k * ωₚ)
    b₂₁ = (-3c[3]ωₚ / (2d₁)) * (3k * ωₚ - 4)
    b₂₂ = -3c[3] * ωₚ / d₁
    d₂₁ = -c[3] / (2ωₚ^2)
    a₃₁ =
        (-9ωₚ / (4d₂)) * (4c[3] * (k * a₂₃ - b₂₁) + k * c[4] * (4 + k^2)) +
        ((9ωₚ^2 + 1 - c[2]) / (2d₂)) * (3c[3] * (2a₂₃ - k * b₂₁) + c[4] * (2 + 3k^2))
    a₃₂ =
        (-9ωₚ / (4d₂)) * (4c[3] * (3k * a₂₄ - b₂₂) + k * c[4]) -
        (3 / (2d₂)) * (9ωₚ^2 + 1 - c[2]) * (c[3] * (k * b₂₂ + d₂₁ - 2a₂₄) - c[4])
    b₃₁ =
        (3 / (8d₂)) * 8ωₚ * (3c[3] * (k * b₂₁ - 2a₂₃) - c[4] * (2 + 3k^2)) +
        (3 / (8d₂)) * (9ωₚ^2 + 1 + 2c[2]) * (4c[3] * (k * a₂₃ - b₂₁) + k * c[4] * (4 + k^2))
    b₃₂ =
        (9ωₚ / d₂) * (c[3] * (k * b₂₂ + d₂₁ - 2a₂₄) - c[4]) +
        (3(9ωₚ^2 + 1 + 2c[2]) / (8d₂) * (4c[3] * (k * a₂₄ - b₂₂) + k * c[4]))
    d₃₁ = (3 / (64ωₚ^2)) * (4c[3] * a₂₄ + c[4])
    d₃₂ = (3 / (64 + ωₚ^2)) * (4c[3] * (a₂₃ - d₂₁) + c[4] * (4 + k^2))

    s₁ =
        (1 / (2ωₚ * (ωₚ * (1 + k^2) - 2k))) * (
            3c[3] / 2 * (2a₂₁ * (k^2 - 2) - a₂₃ * (k^2 + 2) - 2k * b₂₁) -
            (3c[4] / 8) * (3k^4 - 8k^2 + 8)
        )
    s₂ =
        (1 / (2ωₚ * (ωₚ * (1 + k^2) - 2k))) * (
            3c[3] / 2 * (2a₂₂ * (k^2 - 2) + a₂₄ * (k^2 + 2) + 2k * b₂₂ + 5d₂₁) +
            (3c[4] / 8) * (12 - k^2)
        )
    l₁ = (-3c[3] / 2) * (2a₂₁ + a₂₃ + 5d₂₁) - (3c[4] / 8) * (12 - k^2) + 2ωₚ^2 * s₁
    l₂ = (3c[3] / 2) * (a₂₄ - 2a₂₂) + (9c[4] / 8) + 2ωₚ^2 * s₂
    Δ = ωₚ^2 - c[2]

    Aᵧ = Z / γ
    Aₓ = √((-l₂ * Aᵧ^2 - Δ) / l₁)

    ν = 1 + s₁ * Aₓ^2 + s₂ * Aᵧ^2
    Τ = 2π / (ωₚ * ν)
    τ =
        ν .*
        (steps > 1 ? range(0, stop=Τ, length=steps) : range(0, stop=Τ, length=1000))

    if hemisphere == :northern
        m = 1.0
    elseif hemisphere == :southern
        m = 3.0
    else
        throw(ArgumentError("`hemisphere` must be `:northern` or `:southern`."))
    end

    δₘ = 2 - m
    τ₁ = @. ωₚ * τ + ϕ

    x = @. γ * (
        a₂₁ * Aₓ^2 + a₂₂ * Aᵧ^2 - Aₓ * cos(τ₁) +
        (a₂₃ * Aₓ^2 - a₂₄ * Aᵧ^2) * cos(2τ₁) +
        (a₃₁ * Aₓ^3 - a₃₂ * Aₓ * Aᵧ^2) * cos(3τ₁)
    ) + 1 - μ - (L == 1 ? γ : -γ)
    y = @. γ * (
        k * Aₓ * sin(τ₁) +
        (b₂₁ * Aₓ^2 - b₂₂ * Aᵧ^2) * sin(2τ₁) +
        (b₃₁ * Aₓ^3 - b₃₂ * Aₓ * Aᵧ^2) * sin(3τ₁)
    )
    z = @. γ * (
        δₘ * Aᵧ * cos(τ₁) +
        δₘ * d₂₁ * Aₓ * Aᵧ * (cos(2τ₁) - 3) +
        δₘ * (d₃₂ * Aᵧ * Aₓ^2 - d₃₁ * Aᵧ^3) * cos(3τ₁)
    )

    ẋ = @. γ * (
        ωₚ * ν * Aₓ * sin(τ₁) - 2ωₚ * ν * (a₂₃ * Aₓ^2 - a₂₄ * Aᵧ^2) * sin(2τ₁) -
        3ωₚ * ν * (a₃₁ * Aₓ^3 - a₃₂ * Aₓ * Aᵧ^2) * sin(3τ₁)
    )
    ẏ = @. γ * (
        ωₚ * ν * k * Aₓ * cos(τ₁) +
        2ωₚ * ν * (b₂₁ * Aₓ^2 - b₂₂ * Aᵧ^2) * cos(2τ₁) +
        3ωₚ * ν * (b₃₁ * Aₓ^3 - b₃₂ * Aₓ * Aᵧ^2) * cos(3τ₁)
    )
    ż = @. γ * (
        -ωₚ * ν * δₘ * Aᵧ * sin(τ₁) - 2ωₚ * ν * δₘ * d₂₁ * Aₓ * Aᵧ * sin(2τ₁) -
        3ωₚ * ν * δₘ * (d₃₂ * Aᵧ * Aₓ^2 - d₃₁ * Aᵧ^2) * sin(3τ₁)
    )

    return hcat(x, y, z)[1:steps, :], hcat(ẋ, ẏ, ż)[1:steps, :], Τ
end

"""
Returns a `Vector` of `Matrix` values.
Each `Matrix` contains a 3-column nondimensional
position trajectory in the Synodic frame which
represents a Zero Velocity Curve.
"""
function zero_velocity_curves(
    r::AbstractVector,
    v::AbstractVector,
    μ::Real;
    nondimensional_range=range(-2; stop=2, length=1000)
)
    x = nondimensional_range
    y = nondimensional_range

    z = [
        2 * potential_energy([xs, ys, 0.0], μ) - jacobi_constant(r, v, μ) for xs in x,
        ys in y
    ]

    zvc_contour = contours(x, y, z, [0.0])

    curves = Vector{Matrix{Float64}}()

    for zvc_level in Contour.levels(zvc_contour)
        x = Float64[]
        y = Float64[]
        for zvc_line in Contour.lines(zvc_level)
            xs, ys = coordinates(zvc_line)
            if length(x) == length(y) == 0
                x = xs
                y = ys
            else
                x = vcat(x, xs)
                y = vcat(y, ys)
            end
        end

        if length(curves) == 0
            curves = [hcat(x, y)]
        else
            curves = vcat(curves, [hcat(x, y)])
        end
    end

    return curves
end

"""
Calculates the eigenvector associated with the stable manifold of a Monodromy matrix.
"""
function convergent_direction(stm::AbstractMatrix)
    evals, evecs = eigen(stm)
    evals = filter(e -> isreal(e) && isposdef(e), evals) .|> real
    evecs =
        filter(
            x -> !isempty(x),
            map(vec -> filter(x -> all(isreal.(vec)), vec), eachcol(evecs)),
        ) .|> real

    imin = findmin(evals)[2]
    imax = findmax(evals)[2]

    if (evals[imin] * evals[imax]) ≉ 1
        @warn "The dynamics appear to be ill-formed; the minimum and maximum real eigenvalues should be multiplicative inverses of one another. Product equals $(evals[imin] * evals[imax]), not 1."
    end

    return evecs[imin]
end

"""
Calculates the direction associated with the unstable manifold of a Monodromy matrix.
"""
function divergent_direction(stm::AbstractMatrix)
    evals, evecs = eigen(stm)
    evals = filter(e -> isreal(e) && isposdef(e), evals) .|> real
    evecs =
        filter(
            x -> !isempty(x),
            map(vec -> filter(x -> all(isreal.(vec)), vec), eachcol(evecs)),
        ) .|> real

    imin = findmin(evals)[2]
    imax = findmax(evals)[2]

    if (evals[imin] * evals[imax]) ≉ 1
        @warn "The dynamics appear to be ill-formed; the minimum and maximum real eigenvalues should be multiplicative inverses of one another. Product equals $(evals[imin] * evals[imax]), not 1."
    end

    return evecs[imax]
end

"""
Perturbs a Cartesian state along a halo orbit onto the provided direction of the provided manifold.
"""
function perturbation(stm::AbstractMatrix, direction::AbstractVector; eps=1e-8)
    return eps * normalize(stm * normalize(direction))
end

end