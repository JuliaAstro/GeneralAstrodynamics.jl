#
# Common data types and associated functions
#

# A custom `Unitful` dimension & unit for mass parameters.
Unitful.@derived_dimension MassParameter Unitful.𝐋^3/Unitful.𝐓^2

"""
Absract type for orbital coordinate frames.
"""
abstract type AbstractFrame end

# TODO // NOTE a future Julia patch will clarify how types and associated values are parameterized
"""
Abstract type for structures parameterized (in part) by `Unitful.Unit` types.
"""
abstract type AbstractUnitfulStructure{F <: AbstractFloat, LU <: Unitful.LengthFreeUnits, TU <: Unitful.TimeFreeUnits} end

"""
Abstract type for all structures describing orbital states (typically Cartesian or Keplerian).
"""
abstract type AbstractState{F, LU, TU, FR<:AbstractFrame} <: AbstractUnitfulStructure{F, LU, TU} end

"""
Abstract type for all structures describing orbital systems (R2BP systems, CR3BP systems, NBP systems, etc.)
"""
abstract type AbstractSystem{F, LU, TU} <: AbstractUnitfulStructure{F, LU, TU} end

"""
Abstract type for all structures which hold __both__ states and systems (aka the whole `Orbit`!)
"""
abstract type AbstractOrbit{F, LU, TU} <: AbstractUnitfulStructure{F, LU, TU} end

"""
An inertial coordinate frame, centered at _some_ location (typically the center of a body, or a barycenter).
"""
struct Inertial <: AbstractFrame end

"""
A Synodic coordinate frame. Located at the barycenter of a CR3BP system, and rotates with the two Celestial Bodies.
"""
struct Synodic <: AbstractFrame end

"""
A Perifocal coordinate frame. An `Inertial` frame which is 2D (orbit falls completely within the XY axis).
"""
struct Perifocal <: AbstractFrame end

"""
Cartesian state which describes a spacecraft or body's position and velocity with respect to _something_,
as described the the frame `FR`.
"""
mutable struct CartesianState{F, LU, TU, FR} <: AbstractState{F, LU, TU, FR} 
    t::F
    r::SubArray{F, 1, MVector{6, F}, Tuple{UnitRange{Int64}}, true}
    v::SubArray{F, 1, MVector{6, F}, Tuple{UnitRange{Int64}}, true}
    rv::MVector{6,F}

    function CartesianState(r::AbstractVector{R}, v::AbstractVector{V}, epoch::E=0, frame=Bodycentric; 
                            lengthunit=u"km", timeunit=u"s") where {R<:Real, V<:Real, E<:Real}
        F  = promote_type(R, V, E)
        if !(F <: AbstractFloat)
            @warn "Type provided ($(string(F))) is not a float: defaulting to Float64."
            F = Float64
        end
        @assert length(r) == length(v) == 3 "Both arguments must have length equal to 3!"
        rv = MVector{6, F}(r[1], r[2], r[3], v[1], v[2], v[3])
        pos = @views rv[1:3]
        vel = @views rv[4:6]
        return new{F, typeof(lengthunit), typeof(timeunit), frame}(F(epoch), pos, vel, rv)
    end

    function CartesianState(r::AbstractVector{R}, v::AbstractVector{V}, 
                            epoch::Unitful.Time = 0u"s", frame=Bodycentric) where {R<:Unitful.Length, V<:Unitful.Velocity}

        # Get ready to commit a crime... we need the Time unit provided in velocity quantity V

        # V is a quantity with units Length / Time.  Let's reverse that so we have Time instead of Time^(-1)
        TL = inv(unit(V))

        # Now we need to find which index (1 or 2) contains the Time unit
        timeaxis = findfirst(T -> T isa Unitful.Dimension{:Time}, collect(typeof(typeof(TL).parameters[2]).parameters[1]))
        
        # Now that we have the proper index, let's select the time unit
        timeunit = Unitful.FreeUnits{(typeof(TL).parameters[1][timeaxis],), Unitful.𝐓, nothing}()

        # This is easy...
        lengthunit = unit(R)

        # Phew!
        return CartesianState(ustrip.(lengthunit, r), ustrip.(lengthunit/timeunit, v), 
                              ustrip(timeunit, epoch), frame; lengthunit = lengthunit, timeunit = timeunit)
    end

    function CartesianState(cart::CartesianState{F,LU,TU,FR}) where {F,LU,TU,FR}
        return CartesianState(cart.r,  cart.v, cart.t, FR; lengthunit = LU(), timeunit = TU())
    end

    function CartesianState(arr::StaticVector{6})
        return CartesianState(arr[1:3], arr[4:6])
    end
end

"""
Convert `CartesianState` values between different `eltype` and `lengthunit` and `timeunit` types.
"""
function Base.convert(::Type{CartesianState{F,LU,TU}}, cart::CartesianState) where {F,LU,TU}
    r = F.(ustrip.(LU(), position_vector(cart)))
    v = F.(ustrip.(LU()/TU(), velocity_vector(cart)))
    t = F.(ustrip.(TU(), epoch(cart)))
    return CartesianState(r, v, t, coordinateframe(cart); lengthunit = LU(), timeunit = TU())
end

"""
Convert `CartesianState` values between the same `eltype` and `lengthunit` and `timeunit` types.
"""
function Base.convert(::Type{CartesianState{F,LU,TU}}, cart::CartesianState{F,LU,TU}) where {F,LU,TU}
    return cart
end

# Overrides `Unitful` unit conversions for `AbstractUnitfulStructure` instances.
(u::Unitful.LengthFreeUnits)(state::CartesianState{F,LU,TU}) where {F,LU,TU} = convert(CartesianState{F, typeof(u), TU}, state)

# Overrides `Unitful` unit conversions for `AbstractUnitfulStructure` instances.
(u::Unitful.TimeFreeUnits)(state::CartesianState{F,LU,TU}) where {F,LU,TU} = convert(CartesianState{F, LU, typeof(u)}, state)

"""
Print a `CartesianState` to `io`.
"""
function Base.show(io::IO, state::CartesianState{F,LU,TU,FR}) where {F,LU,TU,FR} 
    println(io, "  ", string(FR), " Cartesian State:")
    println("")
    println(io, "    t:  ", state.t, " ", string(TU()))
    println(io, "    r = ", [state.r[1] state.r[2] state.r[3]], " ", string(LU()))
    println(io, "    v = ", [state.v[1] state.v[2] state.v[3]], " ", string(LU()/TU()))
end

"""
Print a `CartesianState` to `io`.
"""
Base.show(io::IO, ::MIME"text/plain", state::CartesianState{F,LU,TU,FR}) where {F,LU,TU,FR} = show(io, state)

"""
A normalized length unit, typically used in normalized CR3B problems.
"""
const NormalizedLengthUnit = Unitful.FreeUnits{(), Unitful.𝐋, nothing}

"""
A normalized time unit, typically used in normalized CR3B problems.
"""
const NormalizedTimeUnit = Unitful.FreeUnits{(), Unitful.𝐓, nothing}

"""
Returns the timepoint associated with a `CartesianState`.
"""
epoch(cart::CartesianState) = cart.t * timeunit(cart)

"""
Returns the first paremeter (commonly labeled as `F`) of a `AbstractUnitfulStructure`.
"""
Base.eltype(::AbstractUnitfulStructure{F}) where F = F

"""
Returns the `Unitful.Length` unit associated with a structure.
"""
lengthunit(::AbstractUnitfulStructure{F, LU, TU}) where {F,LU,TU} = LU()

"""
Returns the `Unitful.Time` unit associated with a structure.
"""
timeunit(::AbstractUnitfulStructure{F, LU, TU}) where {F, LU, TU} = TU()

"""
Returns the `Unitful.Velocity` unit associated with a structure.
"""
velocityunit(state::AbstractUnitfulStructure) = lengthunit(state) / timeunit(state)

"""
Returns the `MassParameter` unit associated with a structure.
"""
massparameterunit(state::AbstractUnitfulStructure) = lengthunit(state)^3 / timeunit(state)^2

"""
Returns the coordinate frame associated with a state.
"""
coordinateframe(::AbstractState{F, LU, TU, FR}) where {F, LU, TU, FR} = FR

"""
Returns the `Unitful` position vector of the `Cartesianstate`.
"""
position_vector(state::CartesianState{F,LU,TU}) where {F, LU, TU} = state.r * lengthunit(state)

"""
Returns the `Unitful` velocity vector of the `CartesianState`.
"""
velocity_vector(state::CartesianState{F,LU,TU}) where {F, LU, TU} = state.v * velocityunit(state)

"""
Returns the `Unitful` scalar position of the `CartesianState`.
"""
scalar_position(state::CartesianState{F,LU,TU}) where {F, LU, TU} = norm(position_vector(state))

"""
Returns the `Unitful` scalar velocity of the `CartesianState`.
"""
scalar_velocity(state::CartesianState{F,LU,TU}) where {F, LU, TU} = norm(velocity_vector(state))