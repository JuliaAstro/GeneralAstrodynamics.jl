#
# Invariant manifolds in the CR3BP
#

function manifold_start(halo::Trajectory, ϵ::Real; manifold=:unstable)
    if !(system(halo) isa CR3BPParameters)
        throw(ArgumentError("The orbit must be a CR3BP orbit!"))
    end

    if !(halo[0] ≈ halo[end])
        throw(ArgumentError("Provided orbit must be periodic!"))
    end

    if manifold ∉ (:unstable, :stable)
        throw(ArgumentError("Keyword argument `manifold` must be set to `:stable` or `:unstable`."))
    end

    T = solution(halo).t[end]

    if halo[0] isa CartesianStateWithSTM && stm(halo[0]) ≈ I
        Φ = deepcopy(stm(halo[end]))
    else
        Φ = let
            initial = let state = deepcopy(state(orbit))
                statevector = CartesianState(
                    initial[1:6]; 
                    lengthunit  = lengthunit(initial), 
                    timeunit    = timeunit(initial), 
                    angularunit = angularunit(initial)
                ) |> CartesianStateWithSTM

                Orbit(statevector, system(halo), initialepoch(halo))
            end
            traj = propagate(initial, T)
            
            deepcopy(stm(traj[end]))
        end
    end

    if manifold == :unstable
        V = unstable_eigenvector(Φ)
    else
        V = stable_eigenvector(Φ)
    end

end