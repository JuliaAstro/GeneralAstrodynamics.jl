#
# Handles NBody problem states.
#

"""
Stores the state of each body in the NBody problem.
"""
struct Body{F<:AbstractFloat} <: AbstractBody
   
    r̅::SVector{3, Unitful.Length{F}}
    v̅::SVector{3, Unitful.Velocity{F}}
    m::Unitful.Mass{F}

    function Body(r̅, v̅, m) where F <: AbstractFloat
    
        return new{F}(SVector{3}(F.(r̅)), SVector{3}(F.(v̅)), F(m))
    
    end

    Body(r̅, v̅, m) = new{Float64}(SVector{3}(Float64.(r̅)), SVector{3}(Float64.(v̅)), Float64(m))

end



"""
Describes a system of `n` `MultibodyStates`'s.
"""
struct MultibodySystem <: OrbitalSystem

    body::Array{Body}

end