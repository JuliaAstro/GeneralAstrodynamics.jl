#
# Calculations related to CR3BP transfers
#

function closest_approach(man::Manifold{<:NormalizedSynodicCR3BPOrbit}, 
                          body_position::AbstractVector  = (secondary_synodic_position ∘ first ∘ first)(man),
                          max_velocity::Unitful.Velocity = 10u"km/s",
                          max_position::Unitful.Length   = 50_000u"km") where O <: NormalizedSynodicCR3BPOrbit

    @assert length(body_position) == 3 "The position vector provided must have 3 elements."
    orbit         = man |> first |> first
    max_velocity  = nondimensionalize_velocity(max_velocity, normalized_length_unit(orbit.system), normalized_time_unit(orbit.system))
    max_position  = nondimensionalize_length(max_position, normalized_length_unit(orbit))

    distance(orb)    = norm(position_vector(orb) - body_position)
    isrealistic(orb) = norm(position_vector(orb)) < max_position && norm(velocity_vector(orb)) < max_velocity

    orbits = filter(isrealistic, vcat(man...))

    !isempty(orbits) || throw(ArgumentError("No state within the manifold meets the `max_velocity` and `max_position` criteria! Use the keyword arguments to loosen the constraints."))
    return findmin(distance.(orbits))[2]
    
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

