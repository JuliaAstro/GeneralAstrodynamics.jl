#
# Provides default frames, and the means for users to create their own.
#


"""
$(TYPEDEF)

A `supertype` for all coordinate frames used in space.
"""
abstract type OrbitalFrame end

"""
$(TYPEDEF)

A `supertype` for all inertial coordinate frames in space.
"""
abstract type Inertial <: OrbitalFrame end

"""
$(TYPEDEF)

A `supertype` for all bodycentric-inertial coordinate frames in space.
"""
abstract type BodycentricInertial <: Inertial end

"""
$(TYPEDEF)

A `supertype` for all barycentric-ineretial coordinate frames in space.
"""
abstract type BarycentricInertial <: Inertial end

"""
$(TYPEDEF)

A `supertype` for all rotating coordinate frames in space.
"""
abstract type Rotating <: OrbitalFrame end

"""
$(TYPEDEF)

A `supertype` for all bodycentric-rotating coordinate frames in space.
"""
abstract type BodycentricRotating <: Rotating end

"""
$(TYPEDEF)

A `supertype` for all barycentric-rotating coordinate frames in space.
"""
abstract type BarycentricRotating <: Rotating end


"""
$(SIGNATURES)

Returns `true` or `false` to indicate whether the provided frame is `Inertial`.
"""
isinertial(::Type{T}) where T <: OrbitalFrame = T <: Inertial 

"""
$(SIGNATURES)

Returns `true` or `false` to indicate whether the provided frame is `Rotating`.
"""
isrotating(::Type{T}) where T <: OrbitalFrame = T <: Rotating

"""
$(SIGNATURES)

Returns `true` or `false` to indicate whether the provided frame is bodycentric.
"""
isbodycentric(::Type{T}) where T <: OrbitalFrame = T <: Union{<:BodycentricInertial, BodycentricRotating}

"""
$(SIGNATURES)

Returns `true` or `false` to indicate whether the provided frame is barycentric.
"""
isbarycentric(::Type{T}) where T <: OrbitalFrame = T <: Union{<:BarycentricInertial, BarycentricRotating}

"""
$(SIGNATURES)

A simple convenience macro for defining new coordinate frames in space.
Creates a new `abstract` sub-type of the final argument, `property`.

# Extended Help

**Usage**

`@frame <name>`
`@frame <name> is <property>`

Here, `<property>` is previously defined `OrbitalFrame`.
which must all be subtypes of `OrbitalFrames.OrbitalFrame`.
By using the shortened syntax, `@frame <name>`, your new `abstract` 
type `<name>` will only have `OrbitalFrames.OrbitalFrame` as
a `supertype`.

**Examples**
```julia
@frame AbstractFrame # second argument defaults to the top-level `OrbitalFrame`
@frame EarthCenteredInertial is BodycentricInertial
@frame SunEarthSynodic is BarycentricInertial
```
"""
macro frame(name)
    quote
        abstract type $name <: $OrbitalFrame end
    end
end
macro frame(name, _is, super)
    @assert :($_is) == :is "Incorrect Usage! Use `is` as the second argument. For example: `@frame ECI is BodycentricInertialFrame`."
    quote
        abstract type $name <: $super end
    end
end
