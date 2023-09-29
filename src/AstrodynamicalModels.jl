"""
Provides astrodynamical models as `AstrodynamicalModels.ODESystems`.
Check out the `ModelingToolkit` docs to learn how to use these
systems for orbit propagation with `DifferentialEquations`, or
see `GeneralAstrodynamics` for some convenient orbit propagation
wrappers.

# Extended help

## License
$(LICENSE)

## Exports
$(EXPORTS)

## Imports
$(IMPORTS)
"""
module AstrodynamicalModels

# Export every model
export R2BSystem, CR3BSystem, NBSystem, PlanarEntrySystem, AttitudeSystem

# Export every `ODEFunction`
export R2BFunction, CR3BFunction, NBFunction, PlanarEntryFunction, AttitudeFunction

# Export every array type
export CartesianState, R2BParameters, CR3BParameters, AttitudeState, AttitudeParameters, PlanarEntryState, PlanarEntryParameters

# Export every orbit type
export R2BOrbit, CR3BOrbit

using Symbolics
using SciMLBase
using Memoize
using LinearAlgebra
using ModelingToolkit
using StaticArrays

using DocStringExtensions
@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(SIGNATURES)

                                         $(DOCSTRING)
                                         """

@template (TYPES, CONSTANTS) = """
                               $(TYPEDEF)

                               $(DOCSTRING)
                               """

"""
An abstract supertype for all astrodynamical state vectors.
"""
abstract type AstrodynamicalState{F,N} <: FieldVector{N,F} end

"""
An abstract supertype for all astrodynamical parameter vectors.
"""
abstract type AstrodynamicalParameters{F,N} <: FieldVector{N,F} end

Base.similar(state::AstrodynamicalState) = typeof(state)(undef)
Base.similar(parameters::AstrodynamicalParameters) = typeof(parameters)(undef)

"""
A mutable vector, with labels, for 6DOF Cartesian states.
"""
mutable struct CartesianState{F} <: AstrodynamicalState{F,6}
    x::F
    y::F
    z::F
    ẋ::F
    ẏ::F
    ż::F

    CartesianState{F}(::UndefInitializer) where {F} = new{F}()
    CartesianState(::UndefInitializer) = CartesianState{Float64}(undef)

    CartesianState{F}(x, y, z, ẋ, ẏ, ż) where {F} = new{F}(convert(F, x), convert(F, y), convert(F, z), convert(F, ẋ), convert(F, ẏ), convert(F, ż))
    CartesianState(x, y, z, ẋ, ẏ, ż) = new{promote_type(typeof(x), typeof(y), typeof(z), typeof(ẋ), typeof(ẏ), typeof(ż))}(x, y, z, ẋ, ẏ, ż)
    CartesianState{F}(; x=zero(F), y=zero(F), z=zero(F), ẋ=zero(F), ẏ=zero(F), ż=zero(F)) where {F} = CartesianState{F}(x, y, z, ẋ, ẏ, ż)
    CartesianState(; x=0.0, y=0.0, z=0.0, ẋ=0.0, ẏ=0.0, ż=0.0) = CartesianState(x, y, z, ẋ, ẏ, ż)
    CartesianState{F}(values::NamedTuple) where {F} =
        let (; x, y, z, ẋ, ẏ, ż) = merge((; x=zero(F), y=zero(F), z=zero(F), ẋ=zero(F), ẏ=zero(F), ż=zero(F)), values)
            CartesianState{F}(x, y, z, ẋ, ẏ, ż)
        end
    CartesianState(values::NamedTuple) =
        let F = Float64, (; x, y, z, ẋ, ẏ, ż) = merge((; x=zero(F), y=zero(F), z=zero(F), ẋ=zero(F), ẏ=zero(F), ż=zero(F)), values)
            CartesianState(x, y, z, ẋ, ẏ, ż)
        end
end

"""
A mutable matrix, with labels, for a 6DOF Cartesian state transition matrix.
"""
mutable struct CartesianSTM{F} <: FieldMatrix{6,6,F}
    xx::F
    xy::F
    xz::F
    xẋ::F
    xẏ::F
    xż::F
    yx::F
    yy::F
    yz::F
    yẋ::F
    yẏ::F
    yż::F
    zx::F
    zy::F
    zz::F
    zẋ::F
    zẏ::F
    zż::F
    ẋx::F
    ẋy::F
    ẋz::F
    ẋẋ::F
    ẋẏ::F
    ẋż::F
    ẏx::F
    ẏy::F
    ẏz::F
    ẏẋ::F
    ẏẏ::F
    ẏż::F
    żx::F
    ży::F
    żz::F
    żẋ::F
    żẏ::F
    żż::F

    CartesianSTM{F}(::UndefInitializer) where {F} = new{F}()
    CartesianSTM(::UndefInitializer) = CartesianSTM{Float64}(undef)

    CartesianSTM{F}(xx, xy, xz, xẋ, xẏ, xż, yx, yy, yz, yẋ, yẏ, yż, zx, zy, zz, zẋ, zẏ, zż, ẋx, ẋy, ẋz, ẋẋ, ẋẏ, ẋż, ẏx, ẏy, ẏz, ẏẋ, ẏẏ, ẏż, żx, ży, żz, żẋ, żẏ, żż) where {F} = new{F}(xx, xy, xz, xẋ, xẏ, xż, yx, yy, yz, yẋ, yẏ, yż, zx, zy, zz, zẋ, zẏ, zż, ẋx, ẋy, ẋz, ẋẋ, ẋẏ, ẋż, ẏx, ẏy, ẏz, ẏẋ, ẏẏ, ẏż, żx, ży, żz, żẋ, żẏ, żż)
    CartesianSTM(xx, xy, xz, xẋ, xẏ, xż, yx, yy, yz, yẋ, yẏ, yż, zx, zy, zz, zẋ, zẏ, zż, ẋx, ẋy, ẋz, ẋẋ, ẋẏ, ẋż, ẏx, ẏy, ẏz, ẏẋ, ẏẏ, ẏż, żx, ży, żz, żẋ, żẏ, żż) = new{promote_type(typeof(xx), typeof(xy), typeof(xz), typeof(xẋ), typeof(xẏ), typeof(xż), typeof(yx), typeof(yy), typeof(yz), typeof(yẋ), typeof(yẏ), typeof(yż), typeof(zx), typeof(zy), typeof(zz), typeof(zẋ), typeof(zẏ), typeof(zż), typeof(ẋx), typeof(ẋy), typeof(ẋz), typeof(ẋẋ), typeof(ẋẏ), typeof(ẋż), typeof(ẏx), typeof(ẏy), typeof(ẏz), typeof(ẏẋ), typeof(ẏẏ), typeof(ẏż), typeof(żx), typeof(ży), typeof(żz), typeof(żẋ), typeof(żẏ), typeof(żż))}(xx, xy, xz, xẋ, xẏ, xż, yx, yy, yz, yẋ, yẏ, yż, zx, zy, zz, zẋ, zẏ, zż, ẋx, ẋy, ẋz, ẋẋ, ẋẏ, ẋż, ẏx, ẏy, ẏz, ẏẋ, ẏẏ, ẏż, żx, ży, żz, żẋ, żẏ, żż)
    CartesianSTM{F}(; xx=one(F), xy=zero(F), xz=zero(F), xẋ=zero(F), xẏ=zero(F), xż=zero(F), yx=zero(F), yy=one(F), yz=zero(F), yẋ=zero(F), yẏ=zero(F), yż=zero(F), zx=zero(F), zy=zero(F), zz=one(F), zẋ=zero(F), zẏ=zero(F), zż=zero(F), ẋx=zero(F), ẋy=zero(F), ẋz=zero(F), ẋẋ=one(F), ẋẏ=zero(F), ẋż=zero(F), ẏx=zero(F), ẏy=zero(F), ẏz=zero(F), ẏẋ=zero(F), ẏẏ=one(F), ẏż=zero(F), żx=zero(F), ży=zero(F), żz=zero(F), żẋ=zero(F), żẏ=zero(F), żż=one(F)) where {F} = CartesianSTM{F}(xx, xy, xz, xẋ, xẏ, xż, yx, yy, yz, yẋ, yẏ, yż, zx, zy, zz, zẋ, zẏ, zż, ẋx, ẋy, ẋz, ẋẋ, ẋẏ, ẋż, ẏx, ẏy, ẏz, ẏẋ, ẏẏ, ẏż, żx, ży, żz, żẋ, żẏ, żż)
    CartesianSTM(; xx=1.0, xy=0.0, xz=0.0, xẋ=0.0, xẏ=0.0, xż=0.0, yx=0.0, yy=1.0, yz=0.0, yẋ=0.0, yẏ=0.0, yż=0.0, zx=0.0, zy=0.0, zz=1.0, zẋ=0.0, zẏ=0.0, zż=0.0, ẋx=0.0, ẋy=0.0, ẋz=0.0, ẋẋ=1.0, ẋẏ=0.0, ẋż=0.0, ẏx=0.0, ẏy=0.0, ẏz=0.0, ẏẋ=0.0, ẏẏ=1.0, ẏż=0.0, żx=0.0, ży=0.0, żz=0.0, żẋ=0.0, żẏ=0.0, żż=1.0) = CartesianSTM(xx, xy, xz, xẋ, xẏ, xż, yx, yy, yz, yẋ, yẏ, yż, zx, zy, zz, zẋ, zẏ, zż, ẋx, ẋy, ẋz, ẋẋ, ẋẏ, ẋż, ẏx, ẏy, ẏz, ẏẋ, ẏẏ, ẏż, żx, ży, żz, żẋ, żẏ, żż)
    CartesianSTM{F}(values::NamedTuple) where {F} =
        let (; xx, xy, xz, xẋ, xẏ, xż, yx, yy, yz, yẋ, yẏ, yż, zx, zy, zz, zẋ, zẏ, zż, ẋx, ẋy, ẋz, ẋẋ, ẋẏ, ẋż, ẏx, ẏy, ẏz, ẏẋ, ẏẏ, ẏż, żx, ży, żz, żẋ, żẏ, żż) = merge((; xx=one(F), xy=zero(F), xz=zero(F), xẋ=zero(F), xẏ=zero(F), xż=zero(F), yx=zero(F), yy=one(F), yz=zero(F), yẋ=zero(F), yẏ=zero(F), yż=zero(F), zx=zero(F), zy=zero(F), zz=one(F), zẋ=zero(F), zẏ=zero(F), zż=zero(F), ẋx=zero(F), ẋy=zero(F), ẋz=zero(F), ẋẋ=one(F), ẋẏ=zero(F), ẋż=zero(F), ẏx=zero(F), ẏy=zero(F), ẏz=zero(F), ẏẋ=zero(F), ẏẏ=one(F), ẏż=zero(F), żx=zero(F), ży=zero(F), żz=zero(F), żẋ=zero(F), żẏ=zero(F), żż=one(F)), values)
            CartesianSTM{F}(xx, xy, xz, xẋ, xẏ, xż, yx, yy, yz, yẋ, yẏ, yż, zx, zy, zz, zẋ, zẏ, zż, ẋx, ẋy, ẋz, ẋẋ, ẋẏ, ẋż, ẏx, ẏy, ẏz, ẏẋ, ẏẏ, ẏż, żx, ży, żz, żẋ, żẏ, żż)
        end
    CartesianSTM(values::NamedTuple) =
        let F = Float64, (; xx, xy, xz, xẋ, xẏ, xż, yx, yy, yz, yẋ, yẏ, yż, zx, zy, zz, zẋ, zẏ, zż, ẋx, ẋy, ẋz, ẋẋ, ẋẏ, ẋż, ẏx, ẏy, ẏz, ẏẋ, ẏẏ, ẏż, żx, ży, żz, żẋ, żẏ, żż) = merge((; xx=one(F), xy=zero(F), xz=zero(F), xẋ=zero(F), xẏ=zero(F), xż=zero(F), yx=zero(F), yy=one(F), yz=zero(F), yẋ=zero(F), yẏ=zero(F), yż=zero(F), zx=zero(F), zy=zero(F), zz=one(F), zẋ=zero(F), zẏ=zero(F), zż=zero(F), ẋx=zero(F), ẋy=zero(F), ẋz=zero(F), ẋẋ=one(F), ẋẏ=zero(F), ẋż=zero(F), ẏx=zero(F), ẏy=zero(F), ẏz=zero(F), ẏẋ=zero(F), ẏẏ=one(F), ẏż=zero(F), żx=zero(F), ży=zero(F), żz=zero(F), żẋ=zero(F), żẏ=zero(F), żż=one(F)), values)
            CartesianSTM(xx, xy, xz, xẋ, xẏ, xż, yx, yy, yz, yẋ, yẏ, yż, zx, zy, zz, zẋ, zẏ, zż, ẋx, ẋy, ẋz, ẋẋ, ẋẏ, ẋż, ẏx, ẏy, ẏz, ẏẋ, ẏẏ, ẏż, żx, ży, żz, żẋ, żẏ, żż)
        end
end

"""
An abstract supertype for all orbits. 

# Extended Help

To support the `AstrodynamicalOrbit` interface, you must implement the following methods.

1. `AstrodynamicalModels.states(orbit)`
2. `AstrodynamicalModels.parameters(orbit)`

"""
abstract type AstrodynamicalOrbit{U,P} end

"""
A full representation of an orbit, including a numerical state, and the parameters of the system.
"""
struct Orbit{U<:AbstractVector,P<:AbstractVector} <: AstrodynamicalOrbit{U,P}
    state::U
    parameters::P

    Orbit(state, parameters) = new{typeof(state),typeof(parameters)}(state, parameters)
    Orbit(; state, parameters) = Orbit(state, parameters)
    Orbit(orbit::Orbit; state=orbit.state, parameters=orbit.parameters) = Orbit(state, parameters)
end

"""
Return the state vector for an `Orbit`.
"""
state(orbit::Orbit) = orbit.state

"""
Return the parameter vector for an `Orbit`.
"""
parameters(orbit::Orbit) = orbit.parameters

include("R2BP.jl")
include("CR3BP.jl")
include("NBP.jl")
include("Entry.jl")
include("Attitude.jl")

end # module
