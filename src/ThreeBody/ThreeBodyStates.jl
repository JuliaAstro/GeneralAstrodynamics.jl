#
# Handles CR3BP problem states.
#

"""
Abstract type for restricted three-body systems.
"""
abstract type RestrictedThreeBodySystem <: AbstractOrbitalSystem end


"""
    `const CR3BPFrames = Union{Synodic, Bodycentric}`

Coordinate frames that are valid within the Circular Restricted Three-body Problem.
"""
const CR3BPFrames = Union{Synodic, Bodycentric}

"""
    `const NormalizedCartesianState{F, FR<:CR3BPFrames} = CartesianState{F, NormalizedLengthUnit, NormalizedTimeUnit, FR}`

CR3BP systems often use normalized units for calculations.
This `NormalizedCartesianState` represents a Cartesian 
state with no dimensioned units.
"""
const NormalizedCartesianState{F, FR<:CR3BPFrames} = CartesianState{F, NormalizedLengthUnit, NormalizedTimeUnit, FR}

function Base.show(io::IO, state::NormalizedCartesianState{F,FR}) where {F,FR} 
    println(io, "  ", "Normalized ", string(FR), " Cartesian State:")
    println("")
    println(io, "    t:  ", state.t)
    println(io, "    r = ", [state.r[1] state.r[2] state.r[3]])
    println(io, "    v = ", [state.v[1] state.v[2] state.v[3]])
end

Base.show(io::IO, ::MIME"text/plain", state::NormalizedCartesianState{F,FR}) where {F,FR} = show(io, state)

"""
    `NormalizedCartesianState(r::AbstractVector{<:Real}, v::AbstractVector{<:Real}, frame::FR = Synodic, epoch::Real = 0) where FR <: CR3BPFrames = CartesianState(r, v, frame, epoch)`

Constructor for a `NormalizedCartesianState`. No `Unitful` arguments allowed!
"""
NormalizedCartesianState(r::AbstractVector{<:Real}, v::AbstractVector{<:Real}, frame::FR = Synodic, epoch::Real = 0) where FR <: CR3BPFrames = CartesianState(r, v, frame, epoch)

"""
    `MinimalCircularRestrictedThreeBodySystem(μ::Real) <: AbstractSystem{F, NormalizedLengthUnit, NormalizedTimeUnit}`

Sometimes you want to do CR3BP calculations with only one system parameter,
the nondimensional mass parameter `μ`. That's what this structure is for!
""" 
mutable struct MinimalCircularRestrictedThreeBodySystem{F} <: AbstractSystem{F, NormalizedLengthUnit, NormalizedTimeUnit}
    μ::F
    name::String

    function MinimalCircularRestrictedThreeBodySystem(μ::Real; name = "") 
        @assert μ ≤ 1//2 "Nondimensional mass parameter must be less than 0.5 by definition."
        F = typeof(μ)
        return new{F}(F(μ), name)
    end
end

function normalized_mass_parameter(state::MinimalCircularRestrictedThreeBodySystem)
    return state.μ
end

function mass_parameter(::MinimalCircularRestrictedThreeBodySystem)
    throw(ArgumentError(
        "Only a normalized mass parameter can be computed from a " * 
        "MinimalCircularRestrictedThreeBodySystem. Try " * 
        "`normalized_mass_parameter(sys)` instead."
    ))
end

function Base.show(io::IO, sys::MinimalCircularRestrictedThreeBodySystem{F}) where {F} 
    println(io, "  Normalized Circular Restricted Three-body System", isempty(sys.name) ? "\n" : " ($(sys.name))\n")
    println(io, "    μ:    ", string(normalized_mass_parameter(sys)))
end

Base.show(io::IO, ::MIME"text/plain", sys::MinimalCircularRestrictedThreeBodySystem{F}) where {F} = show(io, sys)

"""
    `mutable struct CircularRestrictedThreeBodySystem{F, LU, TU} <: AbstractSystem{F, LU, TU}`

This structure represents _all_ required parameters that are used to describe 
a Circular Restricted Three-body System: normalized length unit `DU`, 
normalized time unit `DT`, and a `Tuple` of two mass parameters, `μ`.
"""
mutable struct CircularRestrictedThreeBodySystem{F, LU, TU} <: AbstractSystem{F, LU, TU} 
    DU::F
    DT::F
    μ::Tuple{F,F}
    name::String

    function CircularRestrictedThreeBodySystem(μ::Tuple{Real, Real}, DU::Real, DT::Real = time_scale_factor(DU, μ...), name = ""; lengthunit = u"km", timeunit = u"s")
        F = promote_type(typeof.(μ)..., typeof(DU), typeof(DT))
        if !(F <: AbstractFloat)
            @warn "Promoted type $(string(F)) is not a float. Defaulting to Float64."
            F = Float64
        end
        return new{F, typeof(lengthunit), typeof(timeunit)}(F(DU), F(DT), F.(μ), name)
    end

    function CircularRestrictedThreeBodySystem(μ::Tuple{MassParameter,MassParameter}, DU::Unitful.Length, DT::Unitful.Time = time_scale_factor(DU, μ...), name = "")
        lengthunit = unit(DU)
        timeunit = unit(DT)
        return CircularRestrictedThreeBodySystem(ustrip.(lengthunit^3 / timeunit^2, μ), 
                                                 ustrip(lengthunit, DU), ustrip(timeunit, DT), name; 
                                                 lengthunit = lengthunit, timeunit = timeunit)
    end
end

normalized_distance_unit(sys::CircularRestrictedThreeBodySystem) = sys.DU * lengthunit(sys)
normalized_time_unit(sys::CircularRestrictedThreeBodySystem) = sys.DT * timeunit(sys)

function Base.show(io::IO, sys::CircularRestrictedThreeBodySystem{F,LU,TU}) where {F,LU,TU} 
    println(io, "  Circular Restricted Three-body System", isempty(sys.name) ? "\n" : " ($(sys.name))\n")
    println(io, "        Length Unit:    ", sys.DU, " ", string(lengthunit(sys)))
    println(io, "          Time Unit:    ", sys.DT, " ", string(timeunit(sys)))
    println(io, "    Mass Parameters:    ", "(", max(sys.μ...), ", ", min(sys.μ...), ") ", string(massparameterunit(sys)))
end

Base.show(io::IO, ::MIME"text/plain", sys::CircularRestrictedThreeBodySystem{F,LU,TU}) where {F,LU,TU} = show(io, sys)

normalized_mass_parameter(sys::CircularRestrictedThreeBodySystem) = min(sys.μ) / reduce(+, sys.μ)
mass_parameters(sys::CircularRestrictedThreeBodySystem) = sys.μ .* massparameterunit(sys)
primary_mass_parameter(sys::CircularRestrictedThreeBodySystem) = max(sys.μ...) * massparameterunit(sys)
secondary_mass_parameter(sys::CircularRestrictedThreeBodySystem) = min(sys.μ...) * massparameterunit(sys)

mutable struct CircularRestrictedThreeBodyOrbit{
        F, LU, TU,
        C<:CartesianState{F, <:Unitful.LengthFreeUnits, <:Unitful.TimeFreeUnits, <:CR3BPFrames},
        S<:Union{MinimalCircularRestrictedThreeBodySystem{F}, CircularRestrictedThreeBodySystem{F,LU,TU}}
    } <: AbstractOrbit{F, LU, TU}

    state::C
    system::S

    #=
        Several CR3BP constructors are below. We need to handle the following cases...
           
        Fully described system:
        1) Unitless state vectors, full & unitless system
        2) Unitful state vectors, full & unitful system

        Minimal system:
        3) Unitless state vectors, incomplete & unitless system

    =#

    function CircularRestrictedThreeBodyOrbit(
            r::AbstractVector{<:Real}, v::AbstractVector{<:Real}, 
            μ::Tuple{<:Real, <:Real}, DU::Real, DT::Real, t::Real = 0; 
            frame = Synodic, lengthunit = u"km", timeunit = u"s")     

        F = promote_type(eltype(r), eltype(v), typeof(μ)..., typeof(DU), typeof(DT), typeof(t))
        if !(F <: AbstractFloat)
            @warn "Non-float promoted type $(string(F)) provided. Defaulting to Float64."
            F = Float64
        end

        state = CartesianState(F.(r), F.(v), F(t), frame; lengthunit = lengthunit, timeunit = timeunit)
        system = CircularRestrictedThreeBodySystem(
            F.(ustrip.(lengthunit(state)^3 / timeunit(state)^2, μ)), 
            F(ustrip(lengthunit(state), DU)), 
            F(ustrip(timeunit(state), DT)))

        return new{F, typeof(lengthunit(state)), typeof(timeunit(state)), typeof(state), typeof(system)}(state, system)

    end

    function CircularRestrictedThreeBodyOrbit(
            r::AbstractVector{Unitful.Length}, v::AbstractVector{Unitful.Length}, 
            μ::Tuple{MassParameter, MassParameter}, 
            DU::Unitful.Length, DT::Unitful.Time, t::Unitful.Time = 0u"s"; 
            frame = Synodic)        

        F = promote_type(eltype(ustrip.(r)), eltype(ustrip.(v)), typeof(ustrip.(μ))..., typeof(ustrip(DU)), typeof(ustrip(DT)), typeof(ustrip(t)))
        if !(F <: AbstractFloat)
            @warn "Non-float promoted type $(string(F)) provided. Defaulting to Float64."
            F = Float64
        end
        state = CartesianState(F.(r), F.(v), F(t), frame)
        system = CircularRestrictedThreeBodySystem(
            F.(ustrip.(lengthunit(state)^3 / timeunit(state)^2, μ)), 
            F(ustrip(lengthunit(state), DU)), 
            F(ustrip(timeunit(state), DT)))

        return new{F, typeof(lengthunit(state)), typeof(timeunit(state)), typeof(state), typeof(system)}(state, system)

    end


    function CircularRestrictedThreeBodyOrbit(r::AbstractVector{<:Real}, v::AbstractVector{<:Real}, μ::Real, t::Real = 0)
        F = promote_type(eltype(r), eltype(v), typeof(μ), typeof(t))
        if !(F <: AbstractFloat)
            @warn "Non-float promoted type $(string(F)) provided. Defaulting to Float64."
            F = Float64
        end
        state = CartesianState(F.(r),F.(v),F(t),Synodic; lengthunit = NormalizedLengthUnit(), timeunit = NormalizedTimeUnit())
        system = MinimalCircularRestrictedThreeBodySystem(F(μ))
        return new{F, NormalizedLengthUnit, NormalizedTimeUnit, typeof(state), typeof(system)}(state, system)
    end

    function CircularRestrictedThreeBodyOrbit(r::AbstractVector{<:Real}, v::AbstractVector{<:Real}, system::CircularRestrictedThreeBodySystem, t::Real = 0; frame = Synodic)  
        state = CartesianState(r, v, frame, t; lengthunit = lengthunit(system), timeunit = timeunit(system))   
        F = promote_type(eltype(state), eltype(system))
        L = lengthunit(system) |> typeof
        T = timeunit(system)   |> typeof

        state  = convert(CartesianState{F, L, T}, state)
        system = convert(CircularRestrictedThreeBodySystem{F, L, T}, system)

        return new{F, L, T, typeof(state), typeof(system)}(state, system)   
    end

    function CircularRestrictedThreeBodyOrbit(r::AbstractVector{Unitful.Length}, v::AbstractVector{Unitful.Velocity}, system::CircularRestrictedThreeBodySystem, t::Unitful.Time; frame = Synodic)
        state = CartesianState(r, v, frame, t)   
        F = promote_type(eltype(state), eltype(system))
        L = lengthunit(system) |> typeof
        T = timeunit(system)   |> typeof

        state  = convert(CartesianState{F, L, T}, state)
        system = convert(CircularRestrictedThreeBodySystem{F, L, T}, system)

        return new{F, L, T, typeof(state), typeof(system)}(state, system)   
    end

    function CircularRestrictedThreeBodyOrbit(r::AbstractVector{<:Real}, v::AbstractVector{<:Real}, system::MinimalCircularRestrictedThreeBodySystem, t::Real = 0; frame = Synodic)        
        F = promote_type(eltype(r), eltype(v), eltype(system), typeof(t))
        state = NormalizedCartesianState(F.(r), F.(v), F(t), frame)
        system = MinimalCircularRestrictedThreeBodySystem(F(normalized_mass_parameter(system)))
        return new{F, typeof(NormalizedLengthUnit), typeof(NormalizedTimeUnit), typeof(state), typeof(system)}(state, system)
    end

end

normalized_distance_unit(orb::CircularRestrictedThreeBodyOrbit) = normalized_distance_unit(orb.system)
normalized_time_unit(orb::CircularRestrictedThreeBodyOrbit) = normalized_time_unit(orb.system)
normalized_mass_parameter(orb::CircularRestrictedThreeBodyOrbit) = normalized_mass_parameter(orb.system)
mass_parameters(orb::CircularRestrictedThreeBodyOrbit) = mass_parameters(orb.system)
primary_mass_parameter(orb::CircularRestrictedThreeBodyOrbit) = primary_mass_parameter(orb.system)
secondary_mass_parameter(orb::CircularRestrictedThreeBodyOrbit) = secondary_mass_parameter(orb.system)

function Base.show(io::IO, orb::CircularRestrictedThreeBodyOrbit{F,LU,TU}) where {F,LU,TU} 
    println(io, "Circular Restricted Three-body Orbit\n")
    println(io, orb.state)
    println(io, orb.system)
end

Base.show(io::IO, ::MIME"text/plain", orb::CircularRestrictedThreeBodyOrbit{F,LU,TU}) where {F,LU,TU} = show(io, orb)
