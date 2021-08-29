#
# Iterative Halo orbit and manifold solvers
#

"""
Returns the Monodromy Matrix for a Halo orbit.
"""
function monodromy(orbit::CR3BPOrbit, T; verify=true, reltol = 1e-14, abstol = 1e-14)
    initial = Orbit(CartesianStateWithSTM(state(orbit)), system(orbit), epoch(orbit); frame=frame(orbit))
    final   = propagate(initial, T; reltol=reltol, abstol=abstol, save_everystep=false)[end]
    
    if verify && !(state(initial)[1:6] ≈ final[1:6])
        throw(ErrorException("The provided `orbit` is not periodic!"))
    end

    SMatrix{6,6}(final[7:end]...) |> Matrix
end

"""
Returns true if a `RestrictedThreeBodySystem` is numerically periodic.
"""
function isperiodic(orbit::CR3BPOrbit, T; reltol = 1e-14, abstol = 1e-14)
    final = propagate(orbit, T; reltol=reltol, abstol=abstol, save_everystep=false)[end]
    return state(orbit)[1:6] ≈ final[1:6] 
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

The iterative scheme was pulled from directly from literature
and sample code, including Rund 2018,
and Dr. Mireles' lecture notes and EarthSunHaloOrbit_NewtonMewhod.m 
file available on their website.
Specifically, the half-period iterative scheme (the `F` matrix
in the source code, and corresponding "next guess" computation) 
was ported __directly__ from Dr. Mireles' public code and notes, which are
available online. 
- [Dr. Mireles Notes](http://cosweb1.fau.edu/~jmirelesjames/hw5Notes.pdf)
- [Dr. Mireles Code](http://cosweb1.fau.edu/~jmirelesjames/notes.html)
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/).
"""
function halo(μ; Az=0.0, L=1, hemisphere=:northern,
              tolerance=1e-8, max_iter=100,
              reltol=1e-14, abstol=1e-14,
              nan_on_fail = true, disable_warnings = false)

    r₀, v₀, Τ = analyticalhalo(μ; Az=Az, ϕ=0.0, L=L, hemisphere=hemisphere)
    r₀ = r₀[1,:]
    v₀ = v₀[1,:]
    τ  = Τ/2
    δu = Vector{typeof(τ)}(undef, 6)
    Id = Matrix(I(6))

    for i ∈ 1:max_iter

        problem = ODEProblem(
            AstrodynamicalModels.CR3BPWithSTMVectorField,
            vcat(r₀, v₀, Id...),
            (0.0, τ),
            (μ,)
        )    

        retcode, rₛ, vₛ, Φ = let 
            sols  = solve(problem, Vern9(); reltol=reltol, abstol=abstol)
            final = sols.u[end] 
            sols.retcode, final[1:3], final[4:6], SMatrix{6,6}(final[7:end]...)
        end

        AstrodynamicalModels.CR3BPVectorField(δu, vcat(rₛ, vₛ), (μ,), NaN)
        ∂vₛ = δu[4:6]
        
        # All code in this `if, else, end` block is ported from
        # Dr. Mireles' MATLAB code, which is available on his 
        # website: http://cosweb1.fau.edu/~jmirelesjames/notes.html.
        # Lecture notes, which describe this algorithm further,
        # are available for reference as well:
        # http://cosweb1.fau.edu/~jmirelesjames/hw5Notes.pdf
        if Az ≉ 0
            F = @SMatrix [
                Φ[4,1] Φ[4,5] ∂vₛ[1];
                Φ[6,1] Φ[6,5] ∂vₛ[3];
                Φ[2,1] Φ[2,5]  vₛ[2]
            ]

            xᵪ = SVector(r₀[1], v₀[2], τ) - inv(F) * SVector(vₛ[1], vₛ[3], rₛ[2])

            r₀[1] = xᵪ[1]
            v₀[2] = xᵪ[2]
            τ     = xᵪ[3]
        else
            F = @SMatrix [
                Φ[4,3] Φ[4,5] ∂vₛ[1];
                Φ[6,3] Φ[6,5] ∂vₛ[3];
                Φ[2,3] Φ[2,5]  vₛ[2]
            ]

            xᵪ = SVector(r₀[3], v₀[2], τ) - inv(F) * SVector(vₛ[1], vₛ[3], rₛ[2])

            r₀[3] = xᵪ[1]
            v₀[2] = xᵪ[2]
            τ     = xᵪ[3]
        end

        if abs(vₛ[1]) ≤ tolerance && abs(vₛ[3]) ≤ tolerance
            break;
        elseif τ > 5 * one(τ)
            disable_warnings || @warn "Unreasonably large halo period, $τ, ending iterations."
            !nan_on_fail || return [NaN, NaN, NaN], [NaN, NaN, NaN], NaN
            break
        elseif i == max_iter
            disable_warnings || @warn "Desired tolerance was not reached, and iterations have hit the maximum number of iterations: $max_iter."
            !nan_on_fail || return [NaN, NaN, NaN], [NaN, NaN, NaN], NaN
            break
        elseif retcode != :Success
            !disable_warnings || @warn "Integrator returned $(string(retcode))."
            !nan_on_fail || return [NaN, NaN, NaN], [NaN, NaN, NaN], NaN    
            break        
        end
    end

    return r₀, v₀, 2τ

end

"""
A `halo` wrapper! Returns a `CircularRestrictedThreeBodyOrbit`.
Returns a tuple: `halo_orbit, halo_period`.
"""
function halo(sys::CR3BPParameters, epoch=TAIEpoch(now()); Az=0.0, kwargs...)
    
    if Az isa Unitful.Length
        Az = Az / lengthunit(sys)
    end

    r,v,T = halo(massparameter(sys); Az = Az, kwargs...)
    cart  = CartesianState(vcat(r,v); lengthunit=lengthunit(sys), timeunit=timeunit(sys), angularunit=angularunit(sys))
    orbit = Orbit(cart, sys)
    return orbit, T
end

"""
Given a __periodic__ `CircularRestrictedThreeBodyOrbit`,
returns a nearby `CircularRestrictedThreeBodyOrbit` on the 
__stable manifold__ of the initial state.
"""
function manifold_perturbation(orbit::CR3BPOrbit, V::AbstractVector; eps = 1e-8)
    V = normalize(V)
    r = state(orbit).r + eps * V[1:3]
    v = state(orbit).v + eps * V[4:6]

    cart = CartesianState(vcat(r,v); lengthunit=lengthunit(orbit), timeunit=timeunit(orbit), angularunit=angularunit(orbit))
    
    return Orbit(cart, system(orbit), epoch(orbit))
end

#=
"""
Given a __periodic__ `CircularRestrictedThreeBodyOrbit`,
returns a collection of `Trajectory` instances representing
the stable manifold near orbit.
"""
function stable_manifold(orbit::NormalizedSynodicCR3BPOrbit, T::Real; 
                         duration = T, reltol = 1e-14, 
                         abstol = 1e-14, eps = 1e-8,
                         saveat = 1e-2, num_trajectories = 10, parallel = false)

    @assert duration > 0 "The provided duration cannot be zero or negative. Duration provided was $duration"
    @assert isperiodic(orbit, T) "The provided orbit is not periodic!"

    M = monodromy(orbit, T)
    V = stable_eigenvector(M)

    orbit = NormalizedSynodicSTMCR3BPOrbit(orbit)
    traj = propagate(orbit, T; reltol = reltol, abstol = abstol, saveat = saveat)
    if length(traj) ÷ num_trajectories < 1
        ind = 1:length(traj)
    else
        ind  = 1 : length(traj) ÷ num_trajectories : length(traj) - (length(traj) % num_trajectories)
    end
    
    f = step -> reverse(propagate(manifold(step, step.state.stm * V; eps = eps), -duration; 
                                  reltol = reltol, abstol = abstol))
    return parallel ? pmap(f, traj[ind]) : map(f, traj[ind])
end

"""
Given a __periodic__ `CircularRestrictedThreeBodyOrbit`,
returns a collection of `Trajectory` instances representing
the unstable manifold near orbit.
"""
function unstable_manifold(orbit::NormalizedSynodicCR3BPOrbit, T::Real; 
                           duration = T, reltol = 1e-14, 
                           abstol = 1e-14, eps = 1e-8,
                           num_trajectories = 10,
                           saveat = 1e-2, parallel = false)

    @assert duration > zero(duration) "The provided duration cannot be zero or negative."
    @assert isperiodic(orbit, T) "The provided orbit is not periodic!"

    orbit = NormalizedSynodicSTMCR3BPOrbit(orbit)
    M = monodromy(orbit, T)
    V = unstable_eigenvector(M)

    traj = propagate(orbit, T; reltol = reltol, abstol = abstol, saveat = saveat)
    if length(traj) ÷ num_trajectories < 1
        ind = 1:length(traj)
    else
        ind  = 1 : length(traj) ÷ num_trajectories : length(traj) - (length(traj) % num_trajectories)
    end
    
    f = step -> propagate(manifold(step, step.state.stm * V; eps = eps), duration; 
                          reltol = reltol, abstol = abstol, saveat = saveat)
    return parallel ? pmap(f, traj[ind]) : map(f, traj[ind])
end
=#