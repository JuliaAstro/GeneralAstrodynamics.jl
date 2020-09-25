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

"""
    System

Describes a system of n `BodyState`'s.
"""
struct System

    body::Array{BodyState}

end