#
# Calculations related to CR3BP transfers
#

function closest_approach(man::Manifold{<:NormalizedSynodicCR3BPOrbit}, 
                          body_position::AbstractVector = (secondary_synodic_position ∘ first ∘ first)(man)) where O <: NormalizedSynodicCR3BPOrbit

    @assert length(body_position) == 3 "The position vector provided must have 3 elements."
    distance(orb) = norm(position_vector(orb) - body_position)
    return findmin(distance, map(tuple->tuple[2], findmin.(distance, man)))[2]
    
end

function optimal_approach(man::Manifold{<:NormalizedSynodicCR3BPOrbit},
                          body_position::AbstractVector = (secondary_synodic_position ∘ first ∘ first)(man),
                          cost_function::Function = (r,t)->√(t^2)) where O <: NormalizedSynodicCR3BPOrbit

    @assert length(body_position) == 3 "The position vector provided must have 3 elements."
    distance(orb)  = norm(position_vector(orb) - body_position)
    minimizer(orb) = cost_function(distance(orb), epoch(orb))
    max_distance = nondimensionalize_length(20_000u"km", 
                                            normalized_length_unit(
                                                (first ∘ first)(man).system))
    orbits = filter(orbit->distance(orbit) ≤ max_distance, vcat(man...))
    getindex(orbits, findmin(minimizer.(orbits))[2])

end

