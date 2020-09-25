#
# Handles NBody problem states.
#

"""
    BodyState

Stores the state of each body in the NBody problem.
"""
struct BodyState
   
    r̅::SVector{3, Unitful.Length}
    v̅::SVector{3, Unitful.Velocity}
    m::Unitful.Mass

end

function BodyState(r̅::AbstractArray{Unitful.Length}, v̅::AbstractArray{Unitful.Velocity}, m::Unitful.Mass)
    
    return BodyState(SVector{3}(r̅), SVector{3}(v̅), m)

end

"""
    System

Describes a system of n `BodyState`'s.
"""
struct System

    body::Array{BodyState}

end