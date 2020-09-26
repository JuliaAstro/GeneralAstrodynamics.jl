#
# Handles NBody problem states.
#

"""
    MultibodyState

Stores the state of each body in the NBody problem.
"""
struct Body <: AbstractBody
   
    r̅::SVector{3, Unitful.Length{Float64}}
    v̅::SVector{3, Unitful.Velocity{Float64}}
    m::Unitful.Mass{Float64}

end

function MultibodyState(r̅::AbstractArray{Unitful.Length}, v̅::AbstractArray{Unitful.Velocity}, m::Unitful.Mass)
    
    return MultibodyState(SVector{3}(Float64.(r̅)), SVector{3}(Float64.(v̅)), Float64(m))

end

"""
    System

Describes a system of `n` `MultibodyStates`'s.
"""
struct MultibodySystem <: OrbitalSystem

    body::Array{MultibodyState}

end