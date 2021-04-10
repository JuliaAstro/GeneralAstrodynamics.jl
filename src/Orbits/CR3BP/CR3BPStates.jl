#
# Handles CR3BP problem states.
#

"""
    `const CR3BPFrames = Union{Synodic, Inertial}`

Coordinate frames that are valid within the Circular Restricted Three-body Problem.
"""
const CR3BPFrames = Union{Synodic, Inertial}

"""
    `const NormalizedCartesianState{F, FR<:CR3BPFrames} = CartesianState{F, NormalizedLengthUnit, NormalizedTimeUnit, FR}`

CR3BP systems often use normalized units for calculations.
This `NormalizedCartesianState` represents a Cartesian 
state with no dimensioned units.
"""
const NormalizedCartesianState{F, FR<:CR3BPFrames} = CartesianState{F, NormalizedLengthUnit, NormalizedTimeUnit, FR}

"""
Print `NormalizedCartesianState` instances to `io`.
"""
function Base.show(io::IO, state::NormalizedCartesianState{F,FR}) where {F,FR} 
    println(io, "  ", "Normalized ", string(FR), " Cartesian State:")
    println("")
    println(io, "    t:  ", state.t)
    println(io, "    r = ", [state.r[1] state.r[2] state.r[3]])
    println(io, "    v = ", [state.v[1] state.v[2] state.v[3]])
end

"""
Print `NormalizedCartesianState` instances to `io`.
"""
Base.show(io::IO, ::MIME"text/plain", state::NormalizedCartesianState{F,FR}) where {F,FR} = show(io, state)

"""
A structure which contains the linearization of CR3BP dynamics at timepoint `t`.
"""
mutable struct SynodicCartesianSTMState{F, LU, TU} <: AbstractState{F, LU, TU, Synodic}
    state::CartesianState{F, LU, TU, Synodic}
    stm::MMatrix{6,6,F}

    function SynodicCartesianSTMState(r::AbstractVector{<:Real}, v::AbstractVector{<:Real}, t::Real = 0, stm::AbstractMatrix{<:Real} = I(6); 
                                      lengthunit = NormalizedLengthUnit(), timeunit = NormalizedTimeUnit())
        state = CartesianState(r, v, t, Synodic; lengthunit = lengthunit, timeunit = timeunit)
        F = eltype(state)
        LU = typeof(Orbits.lengthunit(state))
        TU = typeof(Orbits.timeunit(state))

        return new{F, LU, TU}(state, SMatrix{6,6,eltype(state)}(stm))
    end

    function SynodicCartesianSTMState(r::AbstractVector{Unitful.Length}, v::AbstractVector{Unitful.Time}, t::Unitful.Time = 0u"s", stm::AbstractMatrix{<:Real} = I(6))
        state = CartesianState(r, v, t, Synodic)
        return new{eltype(state), typeof(lengthunit(state)), typeof(timeunit(state))}(state, SMatrix{6,6,eltype(state)}(stm))
    end

end


# Overrides `Unitful` unit conversions for `AbstractUnitfulStructure` instances.
(u::Unitful.LengthFreeUnits)(state::SynodicCartesianSTMState{F,LU,TU}) where {F,LU,TU} = convert(SynodicCartesianSTMState{F, typeof(u), TU}, state)


# Overrides `Unitful` unit conversions for `AbstractUnitfulStructure` instances.
(u::Unitful.TimeFreeUnits)(state::SynodicCartesianSTMState{F,LU,TU}) where {F,LU,TU} = convert(SynodicCartesianSTMState{F, LU, typeof(u)}, state)

"""
Print `SynodicCartesianSTMState` instances to `io`.
"""
function Base.show(io::IO, state::SynodicCartesianSTMState{F,LU, TU}) where {F,LU,TU} 
    println(io, state.state,"\n")
    println(io, "  Local Linearization (Φ):\n")
    for row ∈ eachrow(state.stm)
        print(io, "    ")
        for el ∈ row
            print(io, el, " ")
        end
        println(io,"")
    end
end

"""
Print `SynodicCartesianSTMState` instances to `io`.
"""
Base.show(io::IO, ::MIME"text/plain", state::SynodicCartesianSTMState{F,LU,TU}) where {F,LU,TU} = show(io, state)


"""
Returns the timepoint associated with a `CartesianState`.
"""
epoch(cart::SynodicCartesianSTMState) = epoch(cart.state)

"""
Returns the `Unitful` position vector of the `Cartesianstate`.
"""
position_vector(cart::SynodicCartesianSTMState) = position_vector(cart.state)

"""
Returns the `Unitful` velocity vector of the `CartesianState`.
"""
velocity_vector(cart::SynodicCartesianSTMState) = velocity_vector(cart.state)

"""
Returns the `Unitful` scalar position of the `CartesianState`.
"""
scalar_position(cart::SynodicCartesianSTMState) = scalar_position(cart.state)

"""
Returns the `Unitful` scalar velocity of the `CartesianState`.
"""
scalar_velocity(cart::SynodicCartesianSTMState) = scalar_velocity(cart.state)

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

    function MinimalCircularRestrictedThreeBodySystem(μ::Real, name = "") 
        @assert μ ≤ 1//2 "Nondimensional mass parameter must be less than 0.5 by definition."
        F = typeof(μ)
        return new{F}(F(μ), name)
    end
end

"""
Returns the normalized (nondimensional) mass parameter for the CR3BP system.
"""
function normalized_mass_parameter(state::MinimalCircularRestrictedThreeBodySystem)
    return state.μ
end

function mass_parameter(::MinimalCircularRestrictedThreeBodySystem)
    throw(ArgumentError("Only a normalized mass parameter can be computed from a MinimalCircularRestrictedThreeBodySystem. Try `normalized_mass_parameter(sys)` instead."))
end

"""
Print `MinimalCircularRestrictedThreeBodySystem` instances to `io`.
"""
function Base.show(io::IO, sys::MinimalCircularRestrictedThreeBodySystem{F}) where {F} 
    println(io, "  Normalized Circular Restricted Three-body System", isempty(sys.name) ? "\n" : " ($(sys.name))\n")
    println(io, "    μ:    ", string(normalized_mass_parameter(sys)))
end

"""
Print `MinimalCircularRestrictedThreeBodySystem` instances to `io`.
"""
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

"""
Convert between `eltype`, `lengthunit`, and `timeunit` types for CR3BP systems.
"""
function Base.convert(::Type{CircularRestrictedThreeBodySystem{F, LU, TU}}, sys::CircularRestrictedThreeBodySystem) where {F, LU, TU}
    return CircularRestrictedThreeBodySystem(
        F.(ustrip.(LU()^3 / TU()^2, mass_parameters(sys))),
        F.(ustrip.(LU(), normalized_length_unit(sys))),
        F.(ustrip.(TU(), normalized_time_unit(sys))),
        sys.name; lengthunit = LU(), timeunit = TU()
    )
end

"""
Returns the normalized length unit for a CR3BP system.
"""
normalized_length_unit(sys::CircularRestrictedThreeBodySystem) = sys.DU * lengthunit(sys)

"""
Returns the normalized time unit for a CR3BP system.
"""
normalized_time_unit(sys::CircularRestrictedThreeBodySystem) = sys.DT * timeunit(sys)

"""
Print `CircularRestrictedThreeBodySystem` instances to `io`.
"""
function Base.show(io::IO, sys::CircularRestrictedThreeBodySystem{F,LU,TU}) where {F,LU,TU} 
    println(io, "  Circular Restricted Three-body System", isempty(sys.name) ? "\n" : " ($(sys.name))\n")
    println(io, "        Length Unit:    ", sys.DU, " ", string(lengthunit(sys)))
    println(io, "          Time Unit:    ", sys.DT, " ", string(timeunit(sys)))
    println(io, "    Mass Parameters:    ", "(", max(sys.μ...), ", ", min(sys.μ...), ") ", string(massparameterunit(sys)))
end

"""
Print `CircularRestrictedThreeBodySystem` instances to `io`.
"""
Base.show(io::IO, ::MIME"text/plain", sys::CircularRestrictedThreeBodySystem{F,LU,TU}) where {F,LU,TU} = show(io, sys)

"""
Returns the normalized (nondimensional) mass parameter for a CR3BP system.
"""
normalized_mass_parameter(sys::CircularRestrictedThreeBodySystem) = min(sys.μ...) / reduce(+, sys.μ)

"""
Returns the dimensioned (not normalized) mass parameters for a CR3BP system.
"""
mass_parameters(sys::CircularRestrictedThreeBodySystem) = sys.μ .* massparameterunit(sys)

"""
Returns the primary body's mass parameter within a CR3BP system.
"""
primary_mass_parameter(sys::CircularRestrictedThreeBodySystem) = max(sys.μ...) * massparameterunit(sys)

"""
Returns the secondary body's mass parameter within a CR3BP system.
"""
secondary_mass_parameter(sys::CircularRestrictedThreeBodySystem) = min(sys.μ...) * massparameterunit(sys)

"""
An alias for `CartesianState` instances with `Synodic` coordinate frames.
"""
const SynodicCartesianState{F, LU, TU} = CartesianState{F, LU, TU, Synodic}

"""
An alias for `CarteisianState` instances with `Inertial` coordinate frames.
"""
const InertialCartesianState{F, LU, TU} = CartesianState{F, LU, TU, Inertial}

"""
A `Union` of all valid `CR3BP` state types.
"""
const CR3BPState{F, LU, TU, FR<:CR3BPFrames} = Union{SynodicCartesianSTMState{F, LU, TU}, CartesianState{F, LU, TU, FR}}

"""
Contains __all__ relevant state and system information within the CR3BP.
"""
mutable struct CircularRestrictedThreeBodyOrbit{
        F, LU, TU, 
        C<:CR3BPState{F, <:Unitful.LengthFreeUnits, <:Unitful.TimeFreeUnits, <:CR3BPFrames},
        S<:Union{MinimalCircularRestrictedThreeBodySystem{F}, CircularRestrictedThreeBodySystem{F,LU,TU}}
    } <: AbstractOrbit{F, LU, TU}

    state::C
    system::S

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
            DU::Unitful.Length, t::Unitful.Time = 0u"s"; 
            frame = Synodic, name = "")        

        F = promote_type(eltype(ustrip.(r)), eltype(ustrip.(v)), typeof(ustrip.(μ))..., typeof(ustrip(DU)), typeof(ustrip(t)))
        if !(F <: AbstractFloat)
            @warn "Non-float promoted type $(string(F)) provided. Defaulting to Float64."
            F = Float64
        end
        state = CartesianState(F.(r), F.(v), F(t), frame)
        DT = time_scale_factor(DU, μ...)
        system = CircularRestrictedThreeBodySystem(
            F.(ustrip.(lengthunit(state)^3 / timeunit(state)^2, μ)), 
            F(ustrip(lengthunit(state), DU)), 
            F(ustrip(timeunit(state), DT)), 
            name)

        return new{F, typeof(lengthunit(state)), typeof(timeunit(state)), typeof(state), typeof(system)}(state, system)
    end


    function CircularRestrictedThreeBodyOrbit(r::AbstractVector{<:Real}, v::AbstractVector{<:Real}, μ::Real, t::Real = 0; frame = Synodic, name = "")
        F = promote_type(eltype(r), eltype(v), typeof(μ), typeof(t))
        if !(F <: AbstractFloat)
            @warn "Non-float promoted type $(string(F)) provided. Defaulting to Float64."
            F = Float64
        end
        state = CartesianState(F.(r),F.(v),F(t),frame; lengthunit = NormalizedLengthUnit(), timeunit = NormalizedTimeUnit())
        system = MinimalCircularRestrictedThreeBodySystem(F(μ), name)
        return new{F, NormalizedLengthUnit, NormalizedTimeUnit, typeof(state), typeof(system)}(state, system)
    end

    function CircularRestrictedThreeBodyOrbit(r::AbstractVector{<:Real}, v::AbstractVector{<:Real}, system::CircularRestrictedThreeBodySystem, t::Real = 0; frame = Synodic)  
        state = CartesianState(r, v, t, frame; lengthunit = Orbits.lengthunit(system), timeunit = Orbits.timeunit(system))   
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

    function CircularRestrictedThreeBodyOrbit(state::CartesianState, sys::CircularRestrictedThreeBodySystem)
        F = promote_type(eltype(state), eltype(sys))
        state = convert(CartesianState{F, typeof(lengthunit(state)), typeof(timeunit(state))}, state)
        return new{F, typeof(lengthunit(sys)), typeof(timeunit(sys)), typeof(state), typeof(sys)}(state, sys)
    end

    function CircularRestrictedThreeBodyOrbit(state::SynodicCartesianSTMState, sys::CircularRestrictedThreeBodySystem)
        F = promote_type(eltype(state), eltype(sys))
        state = convert(SynodicCartesianSTMState{F, typeof(lengthunit(state)), typeof(timeunit(state))}, state)
        return new{F, typeof(lengthunit(sys)), typeof(timeunit(sys)), typeof(state), typeof(sys)}(state, sys)
    end
end

"""
Convert between `eltype`, `lengthunit`, and `timeunit` types for CR3BP orbits.
"""
function Base.convert(::Type{CircularRestrictedThreeBodyOrbit{F, LU, TU}}, orb::CircularRestrictedThreeBodyOrbit) where {F, LU, TU}
    return CircularRestrictedThreeBodyOrbit(
        F.(ustrip.(LU(), position_vector(orb.state))),
        F.(ustrip.(LU()/TU(), velocity_vector(orb.state))),
        convert(CircularRestrictedThreeBodySystem{F, LU, TU}, orb.sys),
        F(ustrip(TU(), epoch(orb)));
        frame = coordinateframe(orb.state)
    )
end

"""
An alias for `CircularRestrictedThreeBodyOrbit` instances with `NormalizedCartesianState` states.
"""
const NormalizedCR3BPOrbit{F,LU,TU,ST<:NormalizedCartesianState{F,<:AbstractFrame},SR} = CircularRestrictedThreeBodyOrbit{F, LU, TU, ST, SR}

"""
An alias for `CircularRestrictedThreeBodyOrbit` instances with `NormalizedCartesianState` states in the `Synodic` frame.
"""
const NormalizedSynodicCR3BPOrbit{F,LU,TU,SR} = NormalizedCR3BPOrbit{F, LU, TU, NormalizedCartesianState{F,Synodic},SR}

"""
An alias for `CircularRestrictedThreeBodyOrbit`.
"""
const CR3BPOrbit = CircularRestrictedThreeBodyOrbit

"""
an alias for `MinimalCircularRestrictedThreeBodySystem`.
"""
const CR3BPSystem = CircularRestrictedThreeBodySystem

"""
Returns the normalized (nondimensional) length unit associated with a CR3BP orbit.
"""
normalized_length_unit(orb::CircularRestrictedThreeBodyOrbit) = normalized_distance_unit(orb.system)

"""
Returns the normalized (nondimensional) time unit associated with a CR3BP orbit.
"""
normalized_time_unit(orb::CircularRestrictedThreeBodyOrbit) = normalized_time_unit(orb.system)

"""
Returns the normalized (nondimensional) mass parameter associated with a CR3BP orbit.
"""
normalized_mass_parameter(orb::CircularRestrictedThreeBodyOrbit) = normalized_mass_parameter(orb.system)

"""
Returns the dimensioned (not normalized) mass parameters associated with a CR3BP orbit.
"""
mass_parameters(orb::CircularRestrictedThreeBodyOrbit) = mass_parameters(orb.system)

"""
Returns the mass parameter of the primary body within a CR3BP orbit.
"""
primary_mass_parameter(orb::CircularRestrictedThreeBodyOrbit) = primary_mass_parameter(orb.system)

"""
Returns the mass parameter of the primary body within a CR3BP orbit.
"""
secondary_mass_parameter(orb::CircularRestrictedThreeBodyOrbit) = secondary_mass_parameter(orb.system)

"""
Prints a `CircularRestrictedThreeBodyOrbit` to `io`.
"""
function Base.show(io::IO, orb::CircularRestrictedThreeBodyOrbit{F,LU,TU}) where {F,LU,TU} 
    println(io, "Circular Restricted Three-body Orbit\n")
    println(io, orb.state)
    println(io, orb.system)
end

"""
Prints a `CircularRestrictedThreeBodyOrbit` to `io`.
"""
Base.show(io::IO, ::MIME"text/plain", orb::CircularRestrictedThreeBodyOrbit{F,LU,TU}) where {F,LU,TU} = show(io, orb)

"""
Prints a `Trajectory` instance to `io`.
"""
Base.show(io::IO, traj::Trajectory{<:CircularRestrictedThreeBodyOrbit}) = println(io, "Circular Restricted Three-body trajectory with ", length(traj), " steps")
