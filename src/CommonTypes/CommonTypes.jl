"""
    CommonTypes

Contains abstractions for describing orbital states and bodies.
Implementations are provided in TwoBody, and NBody.
"""
module CommonTypes

using Reexport
using DocStringExtensions

@reexport using Unitful, UnitfulAstro

export AbstractBody, OrbitalSystem, PropagationResult
export SOURCECODE

struct SourceCode <: DocStringExtensions.Abbreviation end
const SOURCECODE = SourceCode()

function DocStringExtensions.format(abbrv::SourceCode, buf, doc)
    if include_source_in_docstring
        println("put source in docstring")
        file = doc.data[:path]
        if isfile(file)
            lines = Base.Iterators.drop(eachline(file), doc.data[:linenumber] - 1)
            text = join(lines, '\n')
            _, from = Meta.parse(text, 1; greedy=false)
            _, to = Meta.parse(text, from)
            println(buf, "```julia")
            println(buf, rstrip(text[from:to]))
            println(buf, "```")
        end
    else
        println("no source in docstring")
    end
    return nothing
end

include_source_in_docstring = false
include_sourcecode(b::Bool) = include_source_in_docstring = b

@template DEFAULT = 
"""
$(SIGNATURES)
$(DOCSTRING)
$(SOURCECODE)
"""

""" 
    AbstractBody

Abstract type for bodies in space: both `CelestialBody`s (in
`TwoBody.jl`), and `Body`s (in `NBody.jl`).
"""
abstract type AbstractBody end

"""
    AbstractSystem

Abstract type describing all states in select Astrodynamics problems.
"""
abstract type OrbitalSystem end

"""
    PropagationResult

Abstract type describing a collection of states resulting from 
"""
abstract type PropagationResult end

end
