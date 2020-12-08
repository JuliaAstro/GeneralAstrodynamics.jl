#
# ThreeBodyCalculations.jl
# Calculations for the Circular Restricted
# Three Body Problem.
#

"""
Returns time scale factor, `Tₛ`.
"""
time_scale_factor(a, μ₁, μ₂) = orbital_period(a, μ₁+μ₂)

"""
Returns nondimensional length unit, `DU`.
"""
nondimensionalize_length(rᵢ, a) = upreferred.(rᵢ ./ a)

"""
Returns nondimensional velocity unit, `DU/DT`.
"""
nondimensionalize_velocity(vᵢ, a, Tₛ) = upreferred.(vᵢ ./ (a / Tₛ))

"""
Returns nondimensional time unit, `DT`.
"""
nondimensionalize_time(t, a, μ₁, μ₂) = t / time_scale_factor(a, μ₁, μ₂)

"""
Returns nondimensional mass parameter, `μ`.
"""
nondimensionalize_mass_parameter(μ₁, μ₂) = min(μ₁,μ₂) / (μ₁+μ₂)

"""
Returns nondimensional form of (`Unitful`) scalar posiion.
"""
nondimensionalize(rᵢ::R, a::A) where {
        R<:Length, A<:Length
    } = nondimensionalize_length(rᵢ, a)

"""
Returns nondimensional form of (`Unitful`) position vector.
"""
nondimensionalize(rᵢ::R, a::A) where {
        U<:Length, R<:AbstractVector{U}, A<:Length
    } = nondimensionalize_length(rᵢ, a)

"""
Returns nondimensional form of (`Unitful`) scalar velocity.
"""
nondimensionalize(vᵢ::V, a::A, Tₛ::T) where {
        V<:Velocity, A<:Length, T<:Time
    } = nondimensionalize_velocity(vᵢ, a, Tₛ)

"""
Returns nondimensional form of (`Unitful`) velocity vector.
"""
nondimensionalize(vᵢ::V, a::A, Tₛ::T) where {
        U<:Velocity, V<:AbstractVector{U}, A<:Length, T<:Time
    } = nondimensionalize_velocity(vᵢ, a, Tₛ)

"""
Returns nondimensional form of (`Unitful`) velocity vector.
"""
nondimensionalize(vᵢ::V, a::A, μ₁::U1, μ₂::U2) where {
        U<:Velocity, V<:AbstractVector{U}, A<:Length, U1<:MassParameter, U2<:MassParameter
    } = nondimensionalize_velocity(vᵢ, a, time_scale_factor(a, μ₁, μ₂))

"""
Returns nondimensional form of (`Unitful`) time duration.
"""
nondimensionalize(t::T1, Tₛ::T2) where {
        T1<:Time, T2<:Time
    } = t / Tₛ

"""
Returns nondimensional form of (`Unitful`) time duration.
"""
nondimensionalize(t::T1, a::A, μ₁::U1, μ₂::U2) where {
        T1<:Time, A<:Length, U1<:MassParameter, U2<:MassParameter
    } = nondimensionalize(t, time_scale_factor(a, μ₁, μ₂))

"""
Returns nondimensional form of (`Unitful`) graviational parameters.
"""
nondimensionalize(μ₁::U1, μ₂::U2) where {
        U1<:MassParameter, U2<:MassParameter
    } = min(μ₁, μ₂) / (μ₁+μ₂)

"""
Returns nondimensional Circular Restricted Three-body State.
"""
function nondimensionalize(r₃::R, v₃::V, Δt::T, μ₁::U1, μ₂::U2, a::A) where {
        RT<:Length, R<:AbstractVector{RT},
        VT<:Velocity, V<:AbstractVector{VT},
        T<:Time, U1<:MassParameter, U2<:MassParameter,
        A<:Length
    }

    Tₛ = time_scale_factor(a, μ₁, μ₂)
    return nondimensionalize(r₃, a), 
           nondimensionalize(v₃, a, Tₛ),
           nondimensionalize(Δt, Tₛ),
           nondimensionalize(μ₁, μ₂)

end

"""
Returns dimensional length units.
"""
redimensionalize_length(rᵢ, a) = upreferred(rᵢ .* a)

"""
Returns dimensional velocity units.
"""
redimensionalize_velocity(vᵢ, a, Tₛ) = vᵢ .* (a / Tₛ)

"""
Returns dimensional time unit.
"""
redimensionalize_time(t, a, μ₁, μ₂) = t * time_scale_factor(a, μ₁, μ₂)

"""
Returns dimensional (inertial) form of (`Unitful`) scalar posiion.
"""
redimensionalize(rᵢ::R, a::A) where {
        R<:Length, A<:Length
    } = redimensionalize_length(rᵢ, a)

"""
Returns dimensional (inertial) form of (`Unitful`) position vector.
"""
redimensionalize(rᵢ::R, a::A) where {
        U<:Length, R<:AbstractVector{U}, A<:Length
    } = redimensionalize_length(rᵢ, a)


"""
Returns dimensional (inertial) form of (`Unitful`) scalar velocity.
"""
redimensionalize(vᵢ::V, a::A) where {
        V<:Velocity, A<:Length
    } = redimensionalize_length(vᵢ, a)


"""
Returns dimensional (inertial) form of (`Unitful`) velocity vector.
"""
redimensionalize(vᵢ::V, a::A, Tₛ::T) where {
        U<:Velocity, V<:AbstractVector{U}, A<:Length, T<:Time
    } = redimensionalize_velocity(vᵢ, a, Tₛ)

"""
Returns dimensional (inertial) form of (`Unitful`) velocity vector.
"""
redimensionalize(vᵢ::V, a::A, μ₁::U1, μ₂::U2) where {
        U<:Velocity, V<:AbstractVector{U}, A<:Length, U1<:MassParameter, U2<:MassParameter
    } = redimensionalize_velocity(vᵢ, a, time_scale_factor(a, μ₁, μ₂))

"""
Returns dimensional (inertial) form of (`Unitful`) time duration.
"""
redimensionalize(t::T1, Tₛ::T2) where {
        T1<:Time, T2<:Time
    } = t / Tₛ

"""
Returns dimensional (inertial) form of (`Unitful`) time duration.
"""
redimensionalize(t::T1, a::A, μ₁::U1, μ₂::U2) where {
        T1<:Time, A<:Length, U1<:MassParameter, U2<:MassParameter
    } = redimensionalize(t, time_scale_factor(a, μ₁, μ₂))

"""
Returns the spacecraft's nondimensional position w.r.t. body 1 (or 2).
"""
nondimensional_radius(r, xᵢ=0) = √( (r[1]-xᵢ)^2 + r[2]^2 + r[3]^2 )

"""
Returns the potential energy `U`.
"""
potential_energy(r, μ) = (r[1]^2 + r[2]^2)/2 + ((1-μ)/nondimensional_radius(r,-μ)) + (μ/nondimensional_radius(r,1-μ))

"""
Returns the Jacobi Constant `C`.
"""
jacobi_constant(r, v, μ) = 2*potential_energy(r, μ) - (v⋅v)

"""
Given the Synodic frame vector, returns the vector in the inertial reference frame.
"""
function inertial(vecₛ, t, ω=1.0u"rad"/unit(t))

    θ = ω*t
    ᴵTₛ = @SMatrix [
        cos(θ) sin(θ) 0
       -sin(θ) cos(θ) 0
        0      0      1
    ]

    return  ᴵTₛ * vecₛ

end

"""
Returns the position and velocity vectors in the synodic (rotating) reference frame.
"""
synodic(rᵢ, vᵢ, a, Tₛ) =  nondimensionalize(rᵢ, a), nondimensionalize(vᵢ, a, Tₛ)

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
function lagrange(μ, L=1:5)
    
    if !all(L[i] ∈ (1,2,3,4,5) for i ∈ 1:length(L))
        throw(ArgumentError("Requested lagrange points must be in range [1,5]"))
    end

	expressions = @SVector [
		x -> x - (1-μ)/(x+μ)^2 + μ/(x+μ-1)^2,
		x -> x - (1-μ)/(x+μ)^2 - μ/(x+μ-1)^2,
		x -> x + (1-μ)/(x+μ)^2 + μ/(x+μ+1)^2
	]
	
	return  (map(f->[find_zero(f, (-3,3)), 0, 0], expressions)..., 
			[(1/2) - μ, √(3)/2, 0], [(1/2) - μ, -√(3)/2, 0])[L]
	
end

"""
Returns the partial derivative matrix of potential `U`.

__Arguments:__
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `r`: Non-dimensional position vector for the spacecraft.

__Outputs:__
- Partial derivative matrix of potential `U`.

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/)
"""
function Jᵤ(μ, r, x₀=0.0, y₀=0.0, z₀=0.0)

    x₁ = -μ
    x₂ = 1-μ
    r₁ = nondimensional_radius(r, x₁)
    r₂ = nondimensional_radius(r, x₂)
    x  = @views r[1]
    y  = @views r[2]
    z  = @views r[3]

    #=
    A = (1-μ)/r₁^3 + μ/r₂^3
    B = 3( (1-μ)/r₁^5 + μ/r₂^5 )
    C = 3( (1-μ)/r₁^5 * (x₀ - x₁) + μ/r₂^5 * (x₀ - x₂) )
    
    Uxx = 1 - A + 3(1-μ)*( (x₀-x₁)^2 / r₁^5 ) + 3μ * (x₀ - x₂)^2 / r₂^5
    Uyy = 1 - A + B*y₀^2
    Uzz = -A + B*z₀^2
    Uxy = C*y₀
    Uxz = C * z₀
    Uyz = B * y₀ * z₀
    =#

    Uxx = (1 + (1-μ)*(-r₁^2 + 3(x+μ)^2)) / r₁^5 + (μ * (-r₂^2 + 3(x - 1 + μ)^2)) / r₂^5
    Uyy = 1 + ((1-μ)*(-r₁^2 + 3y^2) / r₁^5) + (μ * (-r₂^2 + 3y^2) / r₂^5)
    Uzz = ((1-μ) * (-r₁^2 + 3z^2) / r₁^5) + (μ * (-r₂^2 + 3z^2) / r₂^5)
    Uxy = (3y  * (1-μ) * (x+μ) / r₁^5) + (μ * (x - 1 + μ) / r₂^5)
    Uxz = (3z * (1-μ)  * (x+μ) / r₁^5)  + (μ * (x - 1 + μ) / r₂^5)
    Uyz = (3y*z*(1-μ) / r₁^5) + μ / r₂^5

    return [
        Uxx Uxy Uxz;
        Uxy Uyy Uyz;
        Uxz Uyz Uzz
    ]

end

"""
Returns the derivative mapping of CR3BP state transition matrix, `F`.

__Arguments:__
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `r`: Non-dimensional position vector for the spacecraft.

__Outputs:__
- Linear mapping from Φ to Φ̇, `F`.

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/)
"""
function state_transition_dynamics(μ, r, L)

    return Matrix(vcat(
        hcat(zeros((3,3)), I(3)),
        hcat(Jᵤ(μ, r, L...), [0 2 0; -2 0 0; 0 0 0])
    ))

end

"""
Returns an analytical solution for a Halo orbit about `L`.

__Arguments:__ 
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `Zₐ`: Desired non-dimensional Z-amplitude for Halo orbit.
- `ϕ`: Desired Halo orbit phase.
- `steps`: Number of non-dimensional timepoints in returned state.
- `L`: Lagrange point to orbit (L1 or L2).
- `hemisphere`: Specifies northern or southern Halo orbit.

__Outputs:__
- Synodic position vector `r::Array{<:AbstractFloat}`
- Synodic velocity vector `v::Array{<:Abstractfloat}`.
- Halo orbit period `Τ`.
- Throws `ArgumentError` if L is not `:L1` or `:L2`.

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/).
"""
function halo_analytic(μ::T1; Zₐ::T2=0.0, ϕ::T3=0.0, steps::T4=1,
                       L::Symbol=:L1, hemisphere::Symbol=:northern) where {
        T1 <: AbstractFloat, 
        T2 <: AbstractFloat, 
        T3 <: Number,
        T4 <: Integer}

    if L == :L1
        point = first(lagrange(μ, 1))
        γ = abs(1 - μ - point)
        n = collect(1:4)
        c = @. (μ + (-1)^n * (1-μ)γ^(n+1)) / (γ^3 * (1 - γ^(n+1)))
    elseif L == :L2
        point = first(lagrange(μ, 2))
        γ = abs(point - 1 + μ)
        n = collect(1:4)
        c = @. ((-1)^n * μ + (-1)^n * (1-μ)γ^(n+1)) / (γ^3 * (1 + γ^(n+1)))
    else
        throw(ArgumentError("Only Halo orbits about L1 or L2 are supported."))
    end

    ωₚ  = √(2 - c[2] + √((9c[2]^2 - 8c[2])/(2)))
    k   = (ωₚ^2 + 1 + 2c[2]) / (2ωₚ)

    d₁  = (3ωₚ^2 / k) * (k*(6ωₚ^2 - 1) - 2ωₚ)
    d₂  = (8ωₚ^2 / k) * (k*((11ωₚ^2 - 1) - 2ωₚ))
    a₂₁ = (3c[3] * (k^2 - 2)) / (4(1 + 2c[2]))
    a₂₂ = (3c[3]) / (4(1 + 2c[2]))
    a₂₃ = (-3c[3]ωₚ / (4k*d₁)) * (3k^3 * ωₚ - 6k*(k-ωₚ) + 4)
    a₂₄ = (-3c[3]ωₚ / (4k*d₁)) * (2 + 3k*ωₚ)
    b₂₁ = (-3c[3]ωₚ / (2d₁)) * (3k*ωₚ - 4)
    b₂₂ = -3c[3]*ωₚ / d₁
    d₂₁ = -c[3] / (2ωₚ^2)
    a₃₁ = (-9ωₚ / (4d₂)) * (4c[3] * (k*a₂₃-b₂₁) + k*c[4]*(4+k^2)) + 
          ((9ωₚ^2 + 1 - c[2]) / (2d₂)) * (3c[3]*(2a₂₃-k*b₂₁) + c[4]*(2+3k^2))
    a₃₂ = (-9ωₚ / (4d₂)) * (4c[3] * (3k*a₂₄-b₂₂) + k*c[4]) -
          (3 / (2d₂)) * (9ωₚ^2 + 1 - c[2]) * (c[3]*(k*b₂₂+d₂₁-2a₂₄) - c[4])
    b₃₁ = (3 / (8d₂)) * 8ωₚ * (3c[3] * (k*b₂₁ - 2a₂₃) - c[4]*(2+3k^2)) + 
          (3/(8d₂)) * (9ωₚ^2 + 1 + 2c[2]) * (4c[3]*(k*a₂₃-b₂₁) + k*c[4]*(4+k^2))
    b₃₂ = (9ωₚ/d₂)*(c[3]*(k*b₂₂+d₂₁-2a₂₄)-c[4]) + 
          (3(9ωₚ^2 + 1 + 2c[2]) / (8d₂) * (4c[3]*(k*a₂₄-b₂₂)+k*c[4]))
    d₃₁ = (3 / (64ωₚ^2)) * (4c[3]*a₂₄ + c[4])
    d₃₂ = (3 / (64 + ωₚ^2)) * (4c[3]*(a₂₃ - d₂₁) + c[4]*(4+k^2))

    s₁  = (1 / (2ωₚ*(ωₚ*(1+k^2) - 2k))) * 
          (3c[3]/2 * (2a₂₁*(k^2 - 2) - a₂₃*(k^2 + 2) - 2k*b₂₁) - 
            (3c[4]/8) * (3k^4 - 8k^2 + 8))
    s₂  = (1 / (2ωₚ*(ωₚ*(1+k^2) - 2k))) * 
          (3c[3]/2 * (2a₂₂*(k^2-2) + a₂₄*(k^2 + 2) + 2k*b₂₂ + 5d₂₁) + 
            (3c[4]/8) * (12 - k^2))
    l₁  = (-3c[3] / 2) * (2a₂₁ + a₂₃ + 5d₂₁) - (3c[4]/8)*(12 - k^2) + 2ωₚ^2 * s₁
    l₂  = (3c[3]/2) * (a₂₄ - 2a₂₂) + (9c[4]/8) + 2ωₚ^2 * s₂
    Δ   = ωₚ^2 - c[2]

    Aᵧ  = Zₐ * γ
    Aₓ  = √((-l₂*Aᵧ^2 - Δ) / l₁)

    ν   = 1 + s₁*Aₓ^2 + s₂*Aᵧ^2
    Τ   = 2π / (ωₚ*ν)
    τ   = ν .* (steps > 1 ? range(0, stop=Τ, length=steps) : range(0, stop=Τ, length=1000))

    if hemisphere == :northern
        m = 1.0
    elseif hemisphere == :southern
        m = 3.0
    else
        throw(ArgumentError("`hemisphere` must be `:northern` or `:southern`."))
    end

    δₘ = 2 - m
    τ₁ = @. ωₚ*τ + ϕ

    x = @. γ * (a₂₁*Aₓ^2 + a₂₂*Aᵧ^2 - Aₓ*cos(τ₁) + (a₂₃*Aₓ^2 - 
                    a₂₄*Aᵧ^2)*cos(2τ₁) + (a₃₁*Aₓ^3 - a₃₂*Aₓ*Aᵧ^2)*cos(3τ₁)) + point
    y = @. γ * (k*Aₓ*sin(τ₁) + (b₂₁*Aₓ^2 - b₂₂*Aᵧ^2)*sin(2τ₁) + 
                    (b₃₁*Aₓ^3 - b₃₂*Aₓ*Aᵧ^2)*sin(3τ₁))
    z = @. γ * (δₘ*Aᵧ*cos(τ₁) + δₘ*d₂₁*Aₓ*Aᵧ*(cos(2τ₁)-3) + 
                    δₘ*(d₃₂*Aᵧ*Aₓ^2 - d₃₁*Aᵧ^3)*cos(3τ₁))

    ẋ = @. γ * (ωₚ*ν*Aₓ*sin(τ₁) - 2ωₚ*ν*(a₂₃*Aₓ^2 - a₂₄*Aᵧ^2)*sin(2τ₁) - 
                    3ωₚ*ν*(a₃₁*Aₓ^3 - a₃₂*Aₓ*Aᵧ^2)*sin(3τ₁))
    ẏ = @. γ * (ωₚ*ν*k*Aₓ*cos(τ₁) + 2ωₚ*ν*(b₂₁*Aₓ^2 - b₂₂*Aᵧ^2)*cos(2τ₁) + 
                    3ωₚ*ν*(b₃₁*Aₓ^3 - b₃₂*Aₓ*Aᵧ^2)*cos(3τ₁))
    ż = @. γ * (-ωₚ*ν*δₘ*Aᵧ*sin(τ₁) - 2ωₚ*ν*δₘ*d₂₁*Aₓ*Aᵧ*sin(2τ₁) - 
                    3ωₚ*ν*δₘ*(d₃₂*Aᵧ*Aₓ^2 - d₃₁*Aᵧ^2)*sin(3τ₁))

    return hcat(x, y, z)[1:steps, :], hcat(ẋ, ẏ, ż)[1:steps, :], Τ

end

"""
Returns a numerical solution for a Halo orbit about `L`.

__Arguments:__ 
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `Zₐ`: Desired non-dimensional Z-amplitude for Halo orbit.
- `ϕ`: Desired Halo orbit phase.
- `L`: Lagrange point to orbit (L1 or L2).
- `hemisphere`: Specifies northern or southern Halo orbit.

__Outputs:__
- Tuple of initial states: `(r, v)` where `r::Vector{<:AbstractFloat}`, `v::Vector{<:Abstractfloat}`.
- Throws `ArgumentError` if L is not `:L1` or `:L2`

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/).
"""
function halo(μ::T1; Zₐ::T2=0.0, ϕ::T3=0.0,
              L::Symbol=:L1, hemisphere::Symbol=:northern,
              tolerance::T4=1e-8, max_iter::T5=10,
              reltol=1e-14, abstol=1e-14) where {
        T1 <: AbstractFloat, 
        T2 <: AbstractFloat, 
        T3 <: Number,
        T4 <: AbstractFloat,
        T5 <: Integer
    }

    μ = promote_type(typeof(μ), typeof(BigFloat(μ)))(μ)
    r₀, v₀, Τ₀ = halo_analytic(μ; Zₐ=Zₐ, ϕ=ϕ, L=L, hemisphere=hemisphere, steps=1000)
    r₀ = r₀[1,:]
    v₀ = v₀[1,:]
    iter = 0
    Τ = Τ₀
    
    if  L == :L1
        point = lagrange(μ, 1)
    else
        point = lagrange(μ, 2)
    end

    @dowhile ((abs(δẋ) ≥ tolerance || abs(δż) ≥ tolerance) && iter < max_iter) begin

        problem = ODEProblem(
            halo_numerical_tic!,
            ComponentArray(rₛ  = r₀,
                           vₛ  = v₀,
                           Φ₁  = [1.0, 0, 0, 0, 0, 0],
                           Φ₂  = [0, 1.0, 0, 0, 0, 0],
                           Φ₃  = [0, 0, 1.0, 0, 0, 0],
                           Φ₄  = [0, 0, 0, 1.0, 0, 0],
                           Φ₅  = [0, 0, 0, 0, 1.0, 0],
                           Φ₆  = [0, 0, 0, 0, 0, 1.0]),
            (0.0, Inf),
            ComponentArray(μ   =  μ, L = point)
        )    

        condition(u, t, integrator) = u.rₛ[2]
        affect!(integrator) = terminate!(integrator)
        halt = ContinuousCallback(condition, affect!)
        integrator = init(problem, Vern9(); callback=halt, reltol=reltol, abstol=abstol)

        solve!(integrator)

        δẋ, δż, r₀, v₀, Τ = reset_halo(
            r₀, v₀, 
            transpose(hcat(integrator.u.Φ₁, integrator.u.Φ₂, integrator.u.Φ₃, integrator.u.Φ₄, integrator.u.Φ₅, integrator.u.Φ₆)), 
            integrator.u.rₛ, integrator.u.vₛ, integrator.t, μ; 
            tol=tolerance
        )

        iter += 1

    end

    if iter == max_iter
        @warn string(
            "Maximum iterations reached: ",
            "Halo orbit may have not converged to desired tolerance of ",
            tolerance
        )
    end

    return r₀, v₀, Τ

end

"""
Returns dynamics tic for combined Halo iterative solver state vector.

__Arguments:__ 
- `∂u`: Derivative of state `u`.
- `u`: State vector: `[x, y, z, ẋ, ẏ, ż, Φ₁, Φ₂, Φ₃, Φ₄, Φ₅, Φ₆]`.
- `p`: Parameters (contains non-dimensional mass parameter `μ`, positions `x₁`, `x₂`, and configuration).
- `t`: Time in seconds.

__Outputs:__
- None (sets derivative `∂u` in-place).

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/).
"""
function halo_numerical_tic!(∂u, u, p, t)

    # Cartesian state
    ∂u.rₛ   =  u.vₛ
    ∂u.vₛ   =  accel(u.rₛ, u.vₛ, p.μ)

    # State transition matrix
    ∂Φ  = state_transition_dynamics(p.μ, u.rₛ, p.L) * Matrix(transpose(hcat(u.Φ₁, u.Φ₂, u.Φ₃, u.Φ₄, u.Φ₅, u.Φ₆)))
    ∂u.Φ₁ = copy(∂Φ[1,:])[:]
    ∂u.Φ₂ = copy(∂Φ[2,:])[:]
    ∂u.Φ₃ = copy(∂Φ[3,:])[:]
    ∂u.Φ₄ = copy(∂Φ[4,:])[:]
    ∂u.Φ₅ = copy(∂Φ[5,:])[:]
    ∂u.Φ₆ = copy(∂Φ[6,:])[:]
    
end

"""
Returns non-dimensional acceleration for CR3BP state.
"""
function accel(rₛ, vₛ, μ)

    x₁ = -μ
    x₂ = 1-μ
    r₁ = nondimensional_radius(rₛ, x₁)
    r₂ = nondimensional_radius(rₛ, x₂)

    return [
         2vₛ[2] + rₛ[1] - (1-μ)*(rₛ[1] - x₁) / r₁^3 - μ*(rₛ[1] - x₂) / r₂^3,
        -2vₛ[1] + rₛ[2] - ((1-μ) / r₁^3 + (μ / r₂^3)) * rₛ[2],
        -((1-μ) / r₁^3 + (μ / r₂^3)) * rₛ[3]
    ]
    
end

"""
Resets the initial conditions for iterative numerical Halo orbit solver.

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/).
- [SciML Documentation](https://diffeq.sciml.ai/stable/features/callback_functions/#callbacks)
"""
function reset_halo(r₀, v₀, Φ, rₛ, vₛ, t, μ; tol=1e-12)

    δẋ = -vₛ[1]
    δż = -vₛ[3]

    ∂vₛ = accel(rₛ, vₛ, μ)

    F = [Φ[4,1] Φ[4,5]; Φ[6,1] Φ[6,5]] - ((1/vₛ[2]) * [∂vₛ[1]; ∂vₛ[3]] * [Φ[2,1] Φ[2,5]])

    δx₀, δẏ₀ = inv(F) * [δẋ; δż] 

    Τₙ = t * 2

    @show [δẋ δż rₛ[2] δx₀ δẏ₀]
    r₀ₙ = r₀ .- [δx₀, 0, 0]
    v₀ₙ = v₀ .- [0, δẏ₀, 0]

    return δẋ, δż, r₀ₙ, v₀ₙ, Τₙ

end
