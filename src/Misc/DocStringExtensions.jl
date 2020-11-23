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

            for (i, method) in enumerate(group)
                N = length(DocStringExtensions.arguments(method))
                # return a list of tuples that represent type signatures
                tuples = DocStringExtensions.find_tuples(typesig)
                # The following will find the tuple that matches the number of arguments in the function
                # ideally we would check that the method signature matches the Tuple{...} signature
                # but that is not straightforward because of how expressive Julia can be
                if Sys.iswindows()
                    t = tuples[findlast(t -> t isa DataType && t <: Tuple && length(t.types) == N, tuples)]
                else
                    t = tuples[findfirst(t -> t isa DataType && t <: Tuple && length(t.types) == N, tuples)]
                end
                DocStringExtensions.printmethod(buf, binding, func, method, t)
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
$(TERSEMETHODS)
$(DOCSTRING)
"""
