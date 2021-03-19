# 
# Provides common DocStringExtensions abbreviations
# for UnitfulAstrodynamics.jl sub-modules. This file
# is NOT meant to be included by anything else! It 
# includes `export` line(s), which will cause an 
# error you this file is included in any non-module
# context.
# 

using DocStringExtensions

struct TerseMethods <: DocStringExtensions.Abbreviation end
const TERSEMETHODS = TerseMethods()

function DocStringExtensions.format(::TerseMethods, buf, doc)
    local binding = doc.data[:binding]
    local typesig = doc.data[:typesig]
    local modname = doc.data[:module]
    local func = Docs.resolve(binding)
    local groups = DocStringExtensions.methodgroups(func, typesig, modname; exact = false)
    if !isempty(groups)
        println(buf)
        println(buf, "```julia")
        for group in groups
            for method in group
                DocStringExtensions.printmethod(buf, binding, func, method)
                println(buf)
            end
        end
        println(buf, "```\n")
        println(buf)
    end
    return nothing
end

struct SourceCode <: DocStringExtensions.Abbreviation end
const SOURCECODE = SourceCode()

function DocStringExtensions.format(abbrv::SourceCode, buf, doc)

    if include_source_in_docstring
        println("Adding source code!")
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
        println("Skipping sourceode.")
    end

    return nothing
    
end

include_source_in_docstring = false
include_sourcecode(b::Bool) = include_source_in_docstring = b

@template (FUNCTIONS, METHODS, MACROS) =
"""
$(METHODLIST)
$(DOCSTRING)
"""
