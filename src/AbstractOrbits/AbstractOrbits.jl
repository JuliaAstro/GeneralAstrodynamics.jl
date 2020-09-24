"""
    AbstractOrbits

Contains abstractions for describing orbital states and bodies.
Implementations are provided in TwoBody, and NBody.
"""
module AbstractOrbits

using Reexport
@reexport using Unitful, UnitfulAstro

export AbstractOrbit, Body, earth, sun

"""
    AbstractOrbit

Abstract type for all orbital states.
"""
abstract type AbstractOrbit end

"""
    Body(μ)

Type representing large bodies in space. Current only Earth 
and the Sun are supported. All bodies are treated as point 
masses. 

"""
struct Body
    μ::Quantity
end

earth = Body(upreferred(1.0u"GMearth"))
sun   = Body(upreferred(1.0u"GMsun"))

end