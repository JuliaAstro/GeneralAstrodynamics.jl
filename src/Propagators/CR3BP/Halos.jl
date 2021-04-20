#
# Iterative Halo orbit and manifold solvers
#


"""
Returns the Monodromy Matrix for a Halo orbit.
"""
function monodromy(orbit::CircularRestrictedThreeBodyOrbit, T; reltol = 1e-14, abstol = 1e-14, atol = 1e-8)
    orb = (NormalizedSynodicSTMCR3BPOrbit ∘ normalize ∘ synodic)(orbit)
    problem = ODEProblem(orb, T)

    integrator = init(problem, Vern9(); reltol=reltol, abstol=abstol)
    solve!(integrator)

    SMatrix{6,6}(integrator.u[7:end]...) |> transpose |> Matrix
end

"""
Returns true if a `RestrictedThreeBodySystem` is numerically periodic.
"""
function isperiodic(orbit::CircularRestrictedThreeBodyOrbit, T; reltol = 1e-14, abstol = 1e-14, atol = 1e-8)
    orb = (normalize ∘ synodic)(orbit)
    problem = ODEProblem(orbit, T)

    integrator = init(problem, Vern9(); reltol=reltol, abstol=abstol)
    solve!(integrator)
        
    return all((isapprox.(orbit.state.r, integrator.u.r; atol = atol)..., isapprox(orbit.state.v, integrator.u.v; atol = atol)...))
end


"""
Returns a numerical solution for a Halo orbit about `L`.

__Arguments:__ 
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `Az`: Desired non-dimensional Z-amplitude for Halo orbit.
- `ϕ`: Desired Halo orbit phase.
- `L`: Lagrange point to orbit (L1 or L2).
- `hemisphere`: Specifies northern or southern Halo orbit.

__Outputs:__
- Tuple of initial states: `(r, v)` where `r::Vector{<:AbstractFloat}`, `v::Vector{<:Abstractfloat}`.
- Throws `ArgumentError` if L is not `:L1` or `:L2`

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/).
"""
function halo(μ; Az=0.0, L=1, hemisphere=:northern,
              tolerance=1e-8, max_iter=20,
              reltol=1e-14, abstol=1e-14,
              nan_on_fail = true)

    r₀, v₀, Τ = analyticalhalo(μ; Az=Az, ϕ=0.0, L=L, hemisphere=hemisphere)
    r₀ = r₀[1,:]
    v₀ = v₀[1,:]
    τ  = Τ/2

    Φ  = Matrix{promote_type(eltype(r₀), eltype(v₀), typeof(τ))}(undef, 6, 6)

    for i ∈ 1:max_iter

        problem = ODEProblem(
            CR3BPSTMTic!,
            ComponentVector(vcat(r₀, v₀, [row for row ∈ eachrow(I(6))]...),
                            Axis(r=1:3, v=4:6, Φ₁=7:12, Φ₂=13:18,
                                 Φ₃=19:24, Φ₄=25:30, Φ₅=31:36, Φ₆=37:42)),
            (0.0, τ),
            (μ = μ,)
        )    

        integrator = init(problem, Vern9(); reltol=reltol, abstol=abstol)
        solve!(integrator)

        rₛ = integrator.u.r
        vₛ = integrator.u.v

        Φ = hcat(integrator.u.Φ₁, integrator.u.Φ₂, integrator.u.Φ₃, integrator.u.Φ₄, integrator.u.Φ₅, integrator.u.Φ₆) |> transpose

        ∂vₛ = accel(rₛ, vₛ, μ)

        if Az ≉ 0
        F = @SMatrix [
            Φ[4,1] Φ[4,5] ∂vₛ[1];
            Φ[6,1] Φ[6,5] ∂vₛ[3];
            Φ[2,1] Φ[2,5]  vₛ[2]
        ]

        TERM1 = @SMatrix [r₀[1]; v₀[2]; τ] 
        TERM2 = - inv(F) * @SMatrix [vₛ[1]; vₛ[3]; rₛ[2]] 
        xᵪ = TERM1 + TERM2

        r₀[1] = xᵪ[1]
        v₀[2] = xᵪ[2]
        τ     = xᵪ[3]
        else
        F = @SMatrix [
            Φ[4,3] Φ[4,5] ∂vₛ[1];
            Φ[6,3] Φ[6,5] ∂vₛ[3];
            Φ[2,3] Φ[2,5]  vₛ[2]
        ]

        TERM1 = @SMatrix [r₀[3]; v₀[2]; τ] 
        TERM2 = - inv(F) * @SMatrix [vₛ[1]; vₛ[3]; rₛ[2]] 
        xᵪ = TERM1 + TERM2

        r₀[3] = xᵪ[1]
        v₀[2] = xᵪ[2]
        τ     = xᵪ[3]
        end

        if abs(integrator.u.v[1]) ≤ tolerance && abs(integrator.u.v[3]) ≤ tolerance
            break;
        elseif i == max_iter
            @warn "Desired tolerance was not reached, and iterations have hit the maximum number of iterations: $max_iter."
            return [NaN, NaN, NaN], [NaN, NaN, NaN], NaN
        end

    end

    return r₀, v₀, 2τ

end

"""
A `halo` wrapper! Returns a `CircularRestrictedThreeBodyOrbit`.
Returns a tuple: `halo_orbit, halo_period`.
"""
function halo(sys::CircularRestrictedThreeBodySystem; kwargs...)
    r,v,T = halo(normalized_mass_parameter(sys); kwargs...)
    orbit = CircularRestrictedThreeBodyOrbit(redimensionalize_length.(r, normalized_length_unit(sys)), 
                                             redimensionalize_velocity.(v, normalized_length_unit(sys), normalized_time_unit(sys)),
                                             sys) |> normalize
    return orbit, T
end

"""
Calculates the eigenvector associated with the __stable manifold__
of a Monodromy matrix.
"""
function stable_eigenvector(monodromy::AbstractMatrix)
    evals, evecs = eigen(monodromy)
    evals = filter(isreal, evals)
    evecs = filter(x->!isempty(x), map(vec->filter(x->all(isreal.(vec)), vec), eachcol(evecs)))

    @assert length(evals) == length(evecs) == 2 "There should only be one real eigenvalue pair."

    if min(real.(evals)...) == evals[1]
        return real.(evecs[1])
    else
        return real.(evecs[2])
    end
end

"""
Calculates the eigenvector associated with the __unstable manifold__
of a Monodromy matrix.
"""
function unstable_eigenvector(monodromy::AbstractMatrix)
    evals, evecs = eigen(monodromy)
    evals = filter(isreal, evals)
    evecs = filter(x->!isempty(x), map(vec->filter(x->all(isreal.(vec)), vec), eachcol(evecs)))

    @assert length(evals) == length(evecs) == 2 "There should only be one real eigenvalue pair."

    if max(real.(evals)...) == evals[1]
        return real.(evecs[1])
    else
        return real.(evecs[2])
    end
end

"""
Given a __periodic__ `CircularRestrictedThreeBodyOrbit`,
returns a nearby `CircularRestrictedThreeBodyOrbit` on the 
__stable manifold__ of the initial state.
"""
function manifold(orbit::NormalizedSynodicCR3BPOrbit, V::AbstractVector; eps = 1e-4)
    r = @views position_vector(orbit) + eps * V[1:3]
    v = @views velocity_vector(orbit) + eps * V[4:6]

    cart = CartesianState(r, v, orbit.state.t, Synodic;
                          lengthunit = NormalizedLengthUnit(), 
                          timeunit = NormalizedTimeUnit())

    return CircularRestrictedThreeBodyOrbit(cart, orbit.system)
end

"""
Given a __periodic__ `CircularRestrictedThreeBodyOrbit`,
returns a collection of `Trajectory` instances representing
the stable manifold near orbit.
"""
function stable_manifold(orbit::NormalizedSynodicCR3BPOrbit, T::Real; 
                         duration = T, reltol = 1e-16, 
                         abstol = 1e-16, eps = 1e-6)

    @assert duration > 0 "The provided duration cannot be zero or negative."
    @assert isperiodic(orbit, T) "The provided orbit is not periodic!"

    M = monodromy(orbit, T)
    V = stable_eigenvector(M)

    return [
        propagate(manifold(step, V; eps = eps), -duration; 
                  reltol = reltol, abstol = abstol)
        for step ∈ propagate(orbit, T; reltol = reltol, abstol = abstol)
    ]
end

"""
Given a __periodic__ `CircularRestrictedThreeBodyOrbit`,
returns a collection of `Trajectory` instances representing
the unstable manifold near orbit.
"""
function unstable_manifold(orbit::NormalizedSynodicCR3BPOrbit, T::Real; 
                           duration = 2T, reltol = 1e-16, 
                           abstol = 1e-16, eps = 1e-6)

    @assert duration > 0 "The provided duration cannot be zero or negative."
    @assert isperiodic(orbit, T) "The provided orbit is not periodic!"

    M = monodromy(orbit, T)
    V = unstable_eigenvector(M)

    return [
        propagate(manifold(step, V; eps = eps), duration; 
                  reltol = reltol, abstol = abstol)
        for step ∈ propagate(orbit, T; reltol = reltol, abstol = abstol)
    ]
end