#
# Invariant manifolds in the CR3BP
#

"""
An alias for some abstract `EnsembleSolution` 
with `States` state vector and parameter 
vector types.
"""
const AbstractOrbitalEnsembleSolution = SciMLBase.AbstractEnsembleSolution{T,N,<:AbstractVector{U}} where {T,N,U<:Trajectory}

"""
An wrapper for a `SciMLBase.ODESolution` with a `GeneralAstrodynamics.States.AbstractState` 
state vector type. This represents an object's `Trajectory` in space!
"""
struct Manifold{FR, S, P, E, O<:AbstractOrbitalEnsembleSolution}
    solution::O
end

"""
Calculates the eigenvector associated with the __stable manifold__
of a Monodromy matrix.
"""
function stable_eigenvector(monodromy::AbstractMatrix; verify=true)
    evals, evecs = eigen(monodromy)
    evals = filter(isreal, evals) .|> real
    evecs = filter(x->!isempty(x), map(vec->filter(x->all(isreal.(vec)), vec), eachcol(evecs))) .|> real

    imin = findmin(evals)[2]
    imax = findmax(evals)[2]

    !verify || @assert (evals[imin] * evals[imax]) ≈ one(eltype(evals)) "Min and max eigenvalue should be multiplicative inverses. Invalid Monodromy matrix. Product equals $(evals[imin] * evals[imax]), not $(one(eltype(evals)))."
    return evecs[imin]
end

"""
Calculates the eigenvector associated with the __unstable manifold__
of a Monodromy matrix.
"""
function unstable_eigenvector(monodromy::AbstractMatrix; verify=true)
    evals, evecs = eigen(monodromy)
    evals = filter(isreal, evals) .|> real
    evecs = filter(x->!isempty(x), map(vec->filter(x->all(isreal.(vec)), vec), eachcol(evecs))) .|> real

    imin = findmin(evals)[2]
    imax = findmax(evals)[2]

    !verify || @assert (evals[imin] * evals[imax]) ≈ one(eltype(evals)) "Min and max eigenvalue should be multiplicative inverses. Invalid Monodromy matrix. Product equals $(evals[imin] * evals[imax]), not $(one(eltype(evals)))."
    return evecs[imax]
end

"""
Perturbs a `CartesianStateWithSTM` in the direction of an 
eigenvector `V`.
"""
function perturb(state::CartesianStateWithSTM, V::AbstractVector; verify=true, eps = 1e-8)
    if verify && length(V) != 6
        throw(ArgumentError("Perturbation vector `V` must have length 6. Vector provided has length $(length(V))"))
    end

    V = normalize(States.get_stm(state) * normalize(V))
    r = States.get_r(state) + eps * V[1:3]
    v = States.get_v(state) + eps * V[4:6]

    return CartesianState(vcat(r,v); lengthunit=lengthunit(state), timeunit=timeunit(state), angularunit=angularunit(state))
end

"""
Returns an orbit perturbed in the direction of the local linearization right-multiplied 
by the provided eigenvector `V`. Only available for `Trajectory` instances with 
`CartesianStateWithSTM` state types.
"""
function perturb(traj::Trajectory{FR, S}, t, V::AbstractVector; eps=1e-8) where {FR, S<:CartesianStateWithSTM}
    state = perturb(state(traj, t), V; eps=eps)
    return Orbit(state, system(traj), epoch(traj,t))
end

"""
Returns an `EnsembleProblem` which represents perturbations 
off of a Halo orbit onto a stable or unstable manifold.
Use `kwarg` `direction=Val{:stable}` or `direction=Val{:unstable}`
to specify whether you want to solve for the stable or unstable 
invariant manifold about the provided Halo orbit.
All `kwargs` arguments are passed directly to `DifferentialEquations`
solvers.
"""
function SciMLBase.EnsembleProblem(halo::Trajectory{FR, <:CartesianStateWithSTM} where FR, duration::Number=(solution(halo).t[end] - solution(halo).t[1]); direction=Val{:stable}, distributed=Val{false}, Trajectory=Val{false}, verify=true, eps=1e-8, kwargs...)

    if verify && halo[1][1:6] ≉ halo[end][1:6]
        throw(ArgumentError("The Halo orbit `Trajectory` provided is not periodic!"))
    end

    if direction ∉ (Val{:stable}, Val{:unstable})
        throw(ArgumentError("Invalid direction provided: $(direction)"))
    end
    
    T₀      = solution(halo).t[1]
    T       = solution(halo).t[end]
    dur     = duration isa Unitful.Time ? duration / timeunit(halo) : duration
    Φ       = monodromy(Orbit(halo, 0), T; verify=verify)
    V       = direction == Val{:stable} ? stable_eigenvector(Φ; verify=verify) : unstable_eigenvector(Φ; verify=verify)
    start   = CartesianState(state(halo, 0))
    nominal = ODEProblem(Orbit(start, system(halo), initialepoch(halo)), dur; kwargs...)

    prob_func   = distributed == Val{true} ? DistributedManifoldIteration(halo, V, dur; verify=verify, eps=eps) : ManifoldIteration(halo, V, dur; verify=verify, eps=eps)
    output_func = Trajectory  == Val{true} ? ManifoldTrajectoryOutput(halo) : ManifoldSolutionOutput(halo)

    return EnsembleProblem(
        nominal; prob_func = prob_func, output_func = output_func, safetycopy = false
    )
end


"""
Returns an `EnsembleProblem` which represents perturbations 
off of a Halo orbit onto a stable or unstable manifold.
Use `kwarg` `direction=Val{:stable}` or `direction=Val{:unstable}`
to specify whether you want to solve for the stable or unstable 
invariant manifold about the provided Halo orbit.
All `kwargs` arguments are passed directly to `DifferentialEquations`
solvers.
"""
function SciMLBase.EnsembleProblem(orbit::Orbit, period::Number, duration::Number=period; kwargs...)
    cart  = state(orbit) isa CartesianStateWithSTM ? CartesianStateWithSTM(CartesianState(state(orbit))) : CartesianStateWithSTM(state(orbit))
    start = Orbit(cart, system(orbit), epoch(orbit); frame=frame(orbit))
    T     = period isa Unitful.Time ? period / timeunit(orbit) : period
    traj  = propagate(start, T)
    return EnsembleProblem(traj, duration; kwargs...)
end


"""
Distributed perturbation for `EnsembleProblem` itration.
"""
function DistributedManifoldIteration(traj::Trajectory, V::AbstractVector, dur::Real; verify=true, eps=1e-8)
    T₀ = solution(traj).t[1]
    T  = solution(traj).t[end]
    return @everywhere function f(prob, i, repeat)
        t = rand() * T
        remake(prob, u0 = perturb(state(traj, t), V; verify=verify, eps=eps), tspan = (T₀+t, T₀+t+dur))
    end
end

"""
Non-distributed perturbation for `EnsembleProblem` itration.
"""
function ManifoldIteration(traj::Trajectory, V::AbstractVector, dur::Real; verify=true, eps=1e-8)
    T₀ = solution(traj).t[1]
    T  = solution(traj).t[end]
    return function f(prob, i, repeat)
        t = rand() * T
        remake(prob, u0 = perturb(state(traj, t), V; verify=verify, eps=eps), tspan = (T₀+t, T₀+t+dur))
    end
end

"""
Maps a solution iteration within an `EnsembleProblem` solver to a `Trajectory`.
"""
function ManifoldTrajectoryOutput(traj::Trajectory)
    return function f(sol,i)
        e = epoch(traj, sol.t[1])
        return (
            Trajectory{frame(traj), typeof(sol.prob.u0), typeof(sol.prob.p), typeof(e), typeof(sol)}(e, sol),
            false
        )
    end
end

"""
Maps a solution iteration within an `EnsembleProblem` solver to a `Trajectory`.
"""
function ManifoldSolutionOutput(traj::Trajectory)
    return function OutputSolution(sol,i)
        return (sol, false)
    end
end

"""
Perturbs a periodic orbit's `Trajectory` in the direction of the stable or
unstable eigenvector of its `monodromy` matrix to form a stable
or unstable `manifold`.
"""
function manifold(orbit::Orbit, period::Number, duration::Number=period; kwargs...)
    start = Orbit(CartesianStateWithSTM(state(orbit)), system(orbit), epoch(orbit); frame=frame(orbit))
    return manifold(propagate(start, period), duration; kwargs...)
end

"""
Perturbs a periodic orbit `traj` in the direction of the stable or
unstable eigenvector of its `monodromy` matrix to form a stable
or unstable `manifold`.
"""
function manifold(traj::Trajectory{FR, S, P, E}, duration::Number=(solution(traj).t[end] - solution(traj).t[1]); verify=true, algorithm=nothing, ensemble_algorithm=nothing, trajectories=100, kwargs...) where {FR, S<:CartesianStateWithSTM, P, E}
    problem  = EnsembleProblem(traj, duration; Trajectory=Val{true}, verify=verify)
    if isnothing(algorithm)
        isnothing(ensemble_algorithm) || @warn "Argument `algorithm` is nothing: `ensemble_algorithm` keyword argument will be ignored."
        solutions = solve(problem; trajectories=trajectories, kwargs...)
    else
        if isnothing(ensemble_algorithm)
            solutions = solve(problem, algorithm; trajectories=trajectories, kwargs...)
        else
            solutions = solve(problem, algorithm, ensemble_algorithm; trajectories=trajectories, kwargs...)
        end
    end

    return Manifold{FR, typeof(initialstate(first(solutions.u))), P, E, typeof(solutions)}(solutions)
end